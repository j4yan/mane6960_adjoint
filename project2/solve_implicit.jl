
include("assign2_example.jl")
include("jacobian.jl")
include("adjoint.jl")

function implicit_solve(degree::Integer, numelem::Integer)
    
    solver, area, q, jac = setup_for_implicit_solve(degree, numelem)

    newton_iterate(solver, area, q, jac)

    postprocessing(solver, area, q)

    return solver, area, q, jac
end

function newton_iterate{T,Tarea,Tres}(solver::EulerSolver{T},
                                      area::AbstractArray{Tarea,2},
                                      q::AbstractArray{Tres,3},
                                      jac::SparseMatrixCSC,
                                      print_res=true)
    if print_res
        println("\nNewton iterate...")
    end
    linSolverType = "direct" # "direct" or "gmres"
    krylov_tol = 1.e-13
    res_tol = 1e-13
    ilu_tol = 1.e-4
    iter_max = 100

    res = zeros(q, Float64)
    q_vec   = reshape(q, length(q))
    res_vec = reshape(res, length(res))

    for iter_newton = 1 : iter_max
        # residual
        calcWeakResidual!(solver, area, q, res)
        l2norm = norm(res_vec)
        if print_res
            println("iter = ", iter_newton-1, ", res_l2norm = ", l2norm)
        end

        if(maximum(l2norm) < res_tol)
            break
        end

        # jacobian
        calcJac(solver, area, q, jac)

        # solve for update
        if linSolverType == "direct"
            # factorize jacobian
            jac_f = factorize(jac)
            dq_vec = Array{Float64}(length(q))
            dq_vec[:] = jac_f \ res_vec
        elseif linSolverType == "gmres"
            # restart GMRES with ILU
            pcILU = crout_ilu(sys.jacMat, τ=ilu_tol)
            dq_vec, hist = gmres(sys.jacMat, res_vec, Pl=pcILU,
                                 restart=50, maxiter=1000, tol=krylov_tol,
                                 log=true)
            printKrylovConvergenceHistory(hist)
        else
            error("Unknown linear solver $linSolverType")
        end

        # update q
        for i = 1 : length(q_vec)
            q_vec[i] -= dq_vec[i]
        end
    end

    return nothing
end


"""
Postprocessing: calculate solution error, J1 and J2
"""
function postprocessing{T,Tarea,Tres}(solver::EulerSolver{T},
                                      area::AbstractArray{Tarea,2},
                                      q::AbstractArray{Tres,3})
    println("\nPostprocessing...")
    # evaluate and print the solution error
    area_left = 2.0
    area_right = 1.5
    area_mid = 1.0
    area_star = 0.8
    γ = Quasi1DEuler.gamma
    temp_stag = 300.0
    press_stag = 100000.0
    rho, rho_u, e = getFlowExact(area_star, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(γ*press/rho_ref)
    qexact = calcExactSolution(area_star, rho_ref, a_ref, area)
    q_float = zeros(q, Float64)
    q_float[:] = real(q)
    L2_err, max_err = calcSolutionError(solver, q_float, qexact)
    println("L2 solution error  = ",L2_err)
    println("max solution error = ",max_err)
    # evaluate and print the functional values
    j1ex = -0.35194635479522557
    j2ex = 0.6896586332699256
    J1 = calcIntegratedSource(solver, area, q)
    J2 = calcOutletPressure(solver, area, q)
    println("J_1: Integrated Source = ", real(J1))
    println("J_2: Outlet Pressure   = ", real(J2))
    J1_error = abs(real(J1) - j1ex)
    J2_error = abs(real(J2) - j2ex)
    println("J1_error = ", J1_error)
    println("J2_error = ", J2_error)

    return L2_err, max_err, J1_error, J2_error
end

"""
Preprocessing before implicit solve, including 
    specifying constant, area, 
    constructing an EulerSolver object,
    allocate q, res, jac
"""
function setup_for_implicit_solve(degree::Integer, numelem::Integer)
    println("\npreprocessing (constructiong solver, area, q, jac)...")
    area_left = 2.0
    area_right = 1.5
    area_mid = 1.0
    area_star = 0.8
    γ = Quasi1DEuler.gamma
    temp_stag = 300.0
    press_stag = 100000.0

    # set the inlet and outlet flow states (for the boundary conditions)
    rho, rho_u, e = getFlowExact(area_star, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(γ*press/rho_ref)
    bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    rho, rho_u, e = getFlowExact(area_star, area_right, temp_stag, press_stag)
    bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    #--------------------------------------------------------------------------------
    # construct the solver object
    # degree = 1   # degree of the polynomial basis
    # numelem = 10 # number of elements
    solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out)

    #--------------------------------------------------------------------------------
    # define the nozzle-area array using the function nozzleArea() defined above
    area = zeros(solver.sbp.numnodes, solver.numelem)
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            area[i,k] = nozzleArea(solver.x[1,i,k])
        end
    end

    #--------------------------------------------------------------------------------
    # initialize the discrete solution to the inlet state (i.e. q = constant)
    Tsln = Complex128
    q = zeros(Tsln, 3, solver.sbp.numnodes, numelem)
    q[1,:,:] = bc_in[1]
    q[2,:,:] = bc_in[2]
    q[3,:,:] = bc_in[3]

    jac = createJac(q)

    return solver, area, q, jac
end
"""
calculate jacobian matrix using coloring complex step method.
"""
function calcJac{T,Tarea,Tres}(solver::EulerSolver{T},
                               area::AbstractArray{Tarea,2},
                               q::AbstractArray{Tres,3},
                               jac::SparseMatrixCSC,
                               nu=0.0)
    fill!(jac.nzval, 0.0)
    n_elem = size(q, 3)
    n_dof_per_elem = size(q, 2)
    n_field = size(q, 1)

    # color and group elements.
    n_color = 3
    flag = [1,2,0]
    group_size = zeros(Int32, n_color)
    group_size[:] = div(n_elem, n_color)
    for i = 1 : n_elem % n_color
        group_size[i] += 1
    end
    groups = zeros(Int32, maximum(group_size), n_color)

    for color = 1 : n_color
        cnt = 1
        for el = 1 : n_elem
            if el % n_color != flag[color]
                continue
            end
            groups[cnt, color] = el
            cnt += 1
        end
    end

    blk_size = n_field * n_dof_per_elem
    eps = 1e-15
    res = zeros(q, Complex128)

    # start the complex step
    for color = 1 : n_color
        for dof = 1 : n_dof_per_elem
            for field = 1 : n_field
                for i = 1 : group_size[color]
                    el = groups[i, color]
                    q[field, dof, el] += Complex128(0., eps)
                end

                calcWeakResidual!(solver, area, q, res, nu=nu)

                for i = 1 : group_size[color]
                    el = groups[i, color]
                    q[field, dof, el] -= Complex128(0., eps)
                end

                for i = 1 : group_size[color]
                    el = groups[i, color]
                    # starting and ending row of nonzero residual
                    # due to disturb in dofs in element \p el
                    row0 = (el - 2) * blk_size + 1
                    row1 = row0 + blk_size*3 - 1
                    if el == 1
                        row0 = 1
                    end
                    if el == n_elem
                        row1 = length(q)
                    end

                    # local jacobian, only one column since we just 
                    # disturbed one filed at a single node
                    dRdq = zeros(Float64, row1 - row0 + 1)

                    for row = row0 : row1
                        loc = row - row0 + 1
                        dRdq[loc] = imag(res[row]) / eps
                    end

                    # assemble
                    offset_i = row0 - 1
                    col = (el - 1) * blk_size + (dof - 1) * n_field + field
                    add_matrix(jac, dRdq, offset_i, col)
                end
            end
        end
    end

    # ********* begin debug **********
    # verifyJac(solver, area, q, jac)
    # ********** end debug ***********
    return nothing
end

function verifyJac{T,Tarea,Tres}(solver::EulerSolver{T},
                                 area::AbstractArray{Tarea,2},
                                 q::AbstractArray{Tres,3},
                                 jac::SparseMatrixCSC)

    eps = 1e-8
    res0 = zeros(q, Float64)
    res1 = zeros(q, Float64)
    calcWeakResidual!(solver, area, q, res0)

    full_jac = zeros(length(q), length(q))
    for j = 1 : length(q)
        q[j] += eps
        calcWeakResidual!(solver, area, q, res1)
        q[j] -= eps
        for i = 1 : length(q)
            full_jac[i, j] = (res1[i] - res0[i]) / eps
        end
    end

    diff = maximum(full_jac - full(jac))
    println("max_diff = $diff")
    if diff > 1e-6
        exit("jacobian matrix is probaly incorrect!")
    end

    return nothing
end

function add_matrix(J::SparseMatrixCSC,
                    Jsub::AbstractArray{Float64, 1},
                    offset_i::Integer,
                    col::Integer)
    for i = 1 : length(Jsub)
        row = i + offset_i

        found = false
        for k = J.colptr[col] : J.colptr[col+1] - 1
            if J.rowval[k] == row
                J.nzval[k] += Jsub[i]
                found = true
                break
            end
        end
        if(!found)
            println(row, ", ", col)
        end
        @assert(found)
    end
    return nothing
end


"""
print out convergence history
"""
function printKrylovConvergenceHistory(hist)
    println("Krylov iteration:")
    println("  isConverged = $(hist.isconverged)")
    println("  n_iters     = $(hist.iters)")
    n_iters = hist.iters
    for i = 1 : n_iters
        println("  iter = $i, res = $(hist[:resnorm][i])")
    end
    return nothing
end

