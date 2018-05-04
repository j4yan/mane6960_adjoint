"""
Solve the adjoint equation. Assumes the jacobian is available.
"""
function adjointSolve{T, Tarea, Tsln}(solver::EulerSolver{T},
                                      area::AbstractArray{Tarea,2},
                                      q::AbstractArray{Tsln, 3},
                                      jac::SparseMatrixCSC,
                                      functional::Function)
    println("\nAdjoint solving...")
    println("  computing dJdq...")
    adj_sln = zeros(q, Float64)
    adj_vec = reshape(adj_sln, length(adj_sln))
    dJdq = calcDJdq(solver, area, q, functional)
    println("  recalculating jacobian matrix...")
    calcJac(solver, area, q, jac)
    println("  transposing jacobian matrix...")
    transposeJacobian(jac)

    linSolverType = "direct"
    krylov_tol = 1.e-13
    ilu_tol = 1.e-3

    println("  solving linear system...")
    if linSolverType == "direct"
        jac_f = factorize(jac)
        adj_vec[:] = jac_f \ dJdq
    elseif linSolverType == "gmres"
        # restart GMRES with ILU
        prILU = crout_ilu(jac, τ=ilu_tol)
        adj_vec, hist = gmres(jac, dJdq, Pl=prILU, restart=50, 
                              maxiter=1000, tol=krylov_tol, log=true)
        printKrylovConvergenceHistory(hist)
    else
        error("Unknown linear solver $linSolverType")
    end

    transposeJacobian(jac)
    return adj_sln
end

"""
Transpose jac. Given that jac is symmetric,
we don't need additional intermediate variables
"""
function transposeJacobian(jac::SparseMatrixCSC)
    for col = 1 : jac.m
        idx0 = jac.colptr[col]
        idx1 = jac.colptr[col + 1] - 1
        for idx = idx0 : idx1
            row = jac.rowval[idx] 
            if row <= col
                continue
            end
            j0 = jac.colptr[row]
            j1 = jac.colptr[row + 1] - 1
            for j = j0 : j1
                row1 = jac.rowval[j]
                if row1 == col
                    tmp = jac.nzval[idx]
                    jac.nzval[idx] = jac.nzval[j]
                    jac.nzval[j] = tmp
                end
            end
        end
    end
    return nothing
end

function calcDJdq{T, Tarea, Tsln}(solver::EulerSolver{T},
                                  area::AbstractArray{Tarea,2},
                                  q::AbstractArray{Tsln, 3},
                                  functional::Function)
    dJdq = zeros(Float64, length(q))
    eps = 1e-15
    
    for i = 1 : length(q)
        q[i] += Complex128(0, eps)
        J = functional(solver, area, q)
        q[i] -= Complex128(0, eps)
        dJdq[i] = imag(J) / eps
    end

    return dJdq
end

"""
Calculates the partial derivative of J wrt area.
"""
function calcdJdA{T, Tarea, Tsln}(solver::EulerSolver{T},
                                  area::AbstractArray{Tarea,2},
                                  q::AbstractArray{Tsln,3},
                                  functional::Function)
    n_node = size(area, 1)
    n_elem = size(area, 2)
    dJdA = zeros(area, Float64)
    eps = 1e-15
    idL = 1
    idR = 1
    xmin = solver.x[1, 1, 1]
    xmax = solver.x[1, 1, 1]
    for n = 2 : size(solver.x, 2)
        if xmin > solver.x[1,n,1]
            xmin = solver.x[1,n,1]
            idL = n
        end
        if xmax < solver.x[1,n,1]
            xmax = solver.x[1,n,1]
            idR = n
        end
    end
    # idL = 1
    # idR = 2

    area_cplx = zeros(area, Complex128)
    for j = 1 : length(area)
        area_cplx[j] = area[j]
    end

    # only loop over interior nodes
    for el = 1 : n_elem
        for n = 1 : n_node
            if n == idL || n == idR
                continue
            end
            area_cplx[n, el] += Complex128(0, eps)
            J = functional(solver, area_cplx, q)
            dJdA[n, el] = imag(J) / eps
            area_cplx[n, el] -= Complex128(0, eps)
        end
    end

    # interface nodes
    for el = 1 : n_elem - 1
        area_cplx[idR, el] += Complex128(0, eps)
        area_cplx[idL, el+1] = area_cplx[idR, el]
        J = functional(solver, area_cplx, q)
        dJdA[idR, el] = imag(J) / eps
        dJdA[idL, el+1] = dJdA[idR, el]
        area_cplx[idR, el] -= Complex128(0, eps)
        area_cplx[idL, el+1] = area_cplx[idR, el]
    end

    # first node
    n = idL
    el = 1
    area_cplx[n, el] += Complex128(0, eps)
    J = functional(solver, area_cplx, q)
    dJdA[n, el] = imag(J) / eps
    area_cplx[n, el] -= Complex128(0, eps)

    # last node
    n = idR
    el = n_elem
    area_cplx[n, el] += Complex128(0, eps)
    J = functional(solver, area_cplx, q)
    dJdA[n, el] = imag(J) / eps
    area_cplx[n, el] -= Complex128(0, eps)

    dJdA_vec = reshape(dJdA, (n_node * n_elem))
    return dJdA_vec
end

"""
Calculates the partial derivative of J wrt area.
"""
function calcdRdA{T, Tarea, Tsln}(solver::EulerSolver{T},
                                  area::AbstractArray{Tarea,2},
                                  q::AbstractArray{Tsln,3})
    res = zeros(q, Complex128)
    dRdA = zeros(Float64, length(res), length(area))
    eps = 1e-15
    idL = 1
    idR = 1
    xmin = solver.x[1, 1, 1]
    xmax = solver.x[1, 1, 1]
    for n = 2 : size(solver.x, 2)
        if xmin > solver.x[1,n,1]
            xmin = solver.x[1,n,1]
            idL = n
        end
        if xmax < solver.x[1,n,1]
            xmax = solver.x[1,n,1]
            idR = n
        end
    end

    area_cplx = zeros(area, Complex128)
    for j = 1 : length(area)
        area_cplx[j] = area[j]
    end

    n_node = size(area, 1)
    n_elem = size(area, 2)

    # only loop over interior nodes
    for el = 1 : n_elem
        for n = 1 : n_node
            if n == idL || n == idR
                continue
            end
            area_cplx[n, el] += Complex128(0, eps)
            calcWeakResidual!(solver, area_cplx, q, res)
            area_cplx[n, el] -= Complex128(0, eps)
            idx = (el -1) * n_node + n
            for i = 1 : length(res)
                dRdA[i,idx] = imag(res[i]) / eps
            end
        end
    end

    # interface nodes
    for el = 1 : n_elem - 1
        area_cplx[idR, el] += Complex128(0, eps)
        area_cplx[idL, el+1] = area_cplx[idR, el]
        calcWeakResidual!(solver, area_cplx, q, res)
        area_cplx[idR, el] -= Complex128(0, eps)
        area_cplx[idL, el+1] = area_cplx[n, el]
        idx1 = (el -1) * n_node + idR
        idx2 = el * n_node + idL
        for i = 1 : length(res)
            dRdA[i,idx1] = imag(res[i]) / eps
            dRdA[i,idx2] = dRdA[i, idx1]
        end
    end

    # first node
    n = idL
    el = 1
    idx = (el - 1) * n_node + n
    area_cplx[n, el] += Complex128(0, eps)
    calcWeakResidual!(solver, area_cplx, q, res)
    area_cplx[n, el] -= Complex128(0, eps)
    for i = 1 : length(res)
        dRdA[i,idx] = imag(res[i]) / eps
    end
    # last node
    n = idR
    el = n_elem
    idx = (el - 1) * n_node + n
    area_cplx[n, el] += Complex128(0, eps)
    calcWeakResidual!(solver, area_cplx, q, res)
    area_cplx[n, el] -= Complex128(0, eps)
    for i = 1 : length(res)
        dRdA[i,idx] = imag(res[i]) / eps
    end

    return dRdA
end
"""
This function compute the total derivative of J wrt A 
using adjoint-based approch. It assumes that we've already
solved the primal problem.
"""
function calcDJDA{T, Tarea, Tsln}(solver::EulerSolver{T},
                                  area::AbstractArray{Tarea, 2},
                                  q::AbstractArray{Tsln, 3},
                                  jac::SparseMatrixCSC,
                                  functional::Function)

    adj_sln = adjointSolve(solver, area, q, jac, functional)
    adj = reshape(adj_sln, length(adj_sln))
    dJdA = calcdJdA(solver, area, q, functional)
    dRdA = calcdRdA(solver, area, q)

    DJDA = dJdA + dRdA'*adj
    return DJDA
end
