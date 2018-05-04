
include("Quasi1DEuler.jl")
include("parameters.jl")
using Quasi1DEuler

shocks = true
max_degree = 2
direc = "results/transonic"
n_elems = [10, 20, 40, 80]
if area_star < 1.0 - 1e-8
    shocks = false
    direc = "results/subsonic"
    max_degree = 3
end
if shocks
    j1ex = -0.38482425832590694
    j2ex = -0.10202215753949603
    n_elems = [10, 15, 20, 30, 40]
else
    j1ex = -0.35194635479522557
    j2ex = -0.11000142657405404
end

if !isdir("results")
    mkdir("results")
end

if !isdir(direc)
    mkdir(direc)
end

# n_elems = [40]

for p = 1 : max_degree
# for p = 1 : 1
    degree = p
    fname1 = string(direc, "/p$(p)_J1.dat")
    fname2 = string(direc, "/p$(p)_J2.dat")
    fout1 = open(fname1, "w")
    fout2 = open(fname2, "w")
    for numelem in n_elems
        vis0 = 10.0/(numelem*degree)
        s0 = -(1.0 + 4.0*log10(degree))
        println("p = $p, numelem = $numelem")
        #---------------------------------------------------
        # solve in coarse space
        solver, area, q, jac = setup_for_implicit_solve(degree, numelem, shocks)
        if shocks
            solveHomotopy!(solver, area, q, jac, maxouter=maxouter, maxinner=10, 
                           tol=tol, inner_tol=1e-4, alpha_max=0.02, 
                           phi_targ=2.0, display=display)
        else
            solveNewton!(solver, area, q, jac, numiter=maxouter, tol=tol)
        end

        J1 = calcIntegratedSource(solver, area, q)
        J2 = calcWeightedSource(solver, area, q)
        J1_error = abs(J1 - j1ex)
        J2_error = abs(J2 - j2ex)

        #---------------------------------------------------
        # solve in fine space
        solver_f, area_f, q_f, jac_f = setup_for_implicit_solve(degree+1, numelem, shocks)
        q_f_proj = zeros(q_f)
        interpSolution!(solver, solver_f, q, q_f_proj)
        # q_f[:] = q_f_proj[:]
        if shocks
            solveHomotopy!(solver_f, area_f, q_f, jac_f, maxouter=maxouter, maxinner=10, 
                           tol=tol, inner_tol=1e-4, alpha_max=0.02, 
                           phi_targ=2.0, display=display)
        else
            solveNewton!(solver_f, area_f, q_f, jac_f, numiter=maxouter, tol=tol)
        end

        # take into account J_H(qH) - J_h(qh)
        J1_h = calcIntegratedSource(solver_f, area_f, q_f_proj)
        J2_h = calcWeightedSource(solver_f, area_f, q_f_proj)
        #---------------------------------------------------
        # adjoint solve in fine space, functional correction
        adj_f = zeros(q_f)
        dJdq_f = zeros(q_f)
        res_f = zeros(q_f)
        calcStateJacobian(solver_f, area_f, q_f, jac_f)
        calcIntegratedSourcedJdq!(solver_f, area_f, q_f, dJdq_f)
        solveAdjoint!(solver_f, area_f, q_f, jac_f, dJdq_f, adj_f)
        calcWeakResidual!(solver_f, area_f, q_f_proj, res_f)
        dJ1 = -dot(adj_f, res_f)
        corrected_J1 = J1_h - dJ1
        fout_elem_error = open(direc * "/elem_J1_error_p$(p)_numelem$(numelem).dat", "w")
        for el = 1 : solver.numelem
            adj_f_elem = view(adj_f, :, :, el)
            res_f_elem = view(res_f, :, :, el)
            adj_vec = reshape(adj_f_elem, (length(adj_f_elem)))
            res_vec = reshape(res_f_elem, (length(res_f_elem)))
            elem_error = abs(dot(res_vec, adj_vec))
            println(fout_elem_error, el, " ", 0.5*(solver.x[1,1,el] + solver.x[1,2,el]), " ", elem_error)
        end
        close(fout_elem_error)

        calcWeightedSourcedJdq!(solver_f, area_f, q_f, dJdq_f)
        solveAdjoint!(solver_f, area_f, q_f, jac_f, dJdq_f, adj_f)
        calcWeakResidual!(solver_f, area_f, q_f_proj, res_f)
        dJ2 = -dot(adj_f, res_f)
        fout_elem_error = open(direc * "/elem_J2_error_p$(p)_numelem$(numelem).dat", "w")
        for el = 1 : solver.numelem
            adj_f_elem = view(adj_f, :, :, el)
            res_f_elem = view(res_f, :, :, el)
            adj_vec = reshape(adj_f_elem, (length(adj_f_elem)))
            res_vec = reshape(res_f_elem, (length(res_f_elem)))
            elem_error = abs(dot(res_vec, adj_vec))
            mid_x = abs(0.5*(solver.x[1,1,el] + solver.x[1,2,el]))
            println(fout_elem_error, el, " ", mid_x, " ", elem_error)
        end
        close(fout_elem_error)
        corrected_J2 = J2_h - dJ2
        corrected_J1_error = abs(corrected_J1 - j1ex)
        corrected_J2_error = abs(corrected_J2 - j2ex)

        println(fout1, 1.0/numelem, " ", 
                J1_error, " ", 
                abs(J1 - corrected_J1), " ", 
                corrected_J1_error)
        println(fout2, 1.0/numelem, " ", 
                J2_error, " ", 
                abs(J2 - corrected_J2), " ", 
                corrected_J2_error)
        flush(fout1)
        flush(fout2)
    end
    close(fout1)
    close(fout2)
end
