include("solve_implicit.jl")
using Quasi1DEuler

if !isdir("results")
    mkdir("results")
end

if !isdir("results/total_deriv")
    mkdir("results/total_deriv")
end

numelems = [20, 40, 80]
degrees = [1, 2, 3, 4]
# numelems = [20]
# degrees = [1,2,3,4]
debug = false

# J1
functional = calcIntegratedSource
# J2
functional = calcOutletPressure

for p in degrees
    for numelem in numelems
        fout = open("results/total_deriv/p$(p)_$(numelem).dat", "w")

        solver, area, q, jac = setup_for_implicit_solve(p, numelem)
        newton_iterate(solver, area, q, jac)
        DJDA = calcDJDA(solver, area, q, jac, functional)

        # sort nodes
        x_s = zeros(solver.x)  
        sortVec!(solver, solver.x, x_s)
        # sort vector to be printed out
        DJDA_s = zeros(DJDA)
        numnodes = solver.sbp.numnodes
        n_elem = size(q, 3)
        DJDA_2D = reshape(DJDA, (numnodes, n_elem))
        DJDA_s_2D = reshape(DJDA_s, (numnodes, n_elem))
        sortVec!(solver, DJDA_2D, DJDA_s_2D)
        for i = 1 : length(x_s)
            println(fout, x_s[i], " ", DJDA_s[i])
        end
        close(fout)

        # ********* begin debug **********
        if !debug
            continue
        end
        solver, area, q, jac = setup_for_implicit_solve(p, numelem)
        print_res = false
        newton_iterate(solver, area, q, jac, print_res)
        J0 = functional(solver, area, q)
        eps = 1e-5
        idL = 1
        idR = 2

        for el = 2 : size(q, 3) - 1
            # 
            for id = 3 : 5
                if id > p + 1
                    continue
                end

                area[id, el] -= eps
                newton_iterate(solver, area, q, jac, print_res)
                Jm = functional(solver, area, q)

                area[id, el] += 2*eps
                newton_iterate(solver, area, q, jac, print_res)
                Jp = functional(solver, area, q)

                area[id, el] -= eps
                idx = (el -1) * size(area, 1) + id
                deriv = (Jp - Jm) / (2*eps)
                # println("idx = $idx, deriv = $deriv")
                diff = abs(DJDA[idx] + deriv) / abs(DJDA[idx])
                println("diff = $diff")
            end

            # 
            id = 1
            area[id, el] -= eps
            area[idR, el-1] -= eps
            newton_iterate(solver, area, q, jac, print_res)
            Jm = functional(solver, area, q)

            area[id, el] += 2*eps
            area[idR, el-1] += 2*eps
            newton_iterate(solver, area, q, jac, print_res)
            Jp = functional(solver, area, q)

            area[id, el] -= eps
            area[idR, el-1] -= eps
            idx = (el -1) * size(area, 1) + id
            deriv = (Jp - Jm) / (2*eps)
            # println("idx = $idx, deriv = $deriv")
            diff = abs(DJDA[idx] + deriv) / abs(DJDA[idx])
            println("diff = $diff")

            # 
            id = 2
            area[id, el] -= eps
            area[idL, el+1] -= eps
            newton_iterate(solver, area, q, jac, print_res)
            Jm = functional(solver, area, q)

            area[id, el] += 2*eps
            area[idL, el+1] += 2*eps
            newton_iterate(solver, area, q, jac, print_res)
            Jp = functional(solver, area, q)

            area[id, el] -= eps
            area[idL, el+1] -= eps
            idx = (el -1) * size(area, 1) + id
            deriv = (Jp - Jm) / (2*eps)
            # println("idx = $idx, deriv = $deriv")
            diff = abs(DJDA[idx] + deriv) / abs(DJDA[idx])
            println("diff = $diff")
        end
        # ********** end debug ***********
    end
end

