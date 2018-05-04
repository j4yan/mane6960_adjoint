export setup_for_implicit_solve, implicit_solve
export calcStateJacobian
export transposeJacobian
export postprocessing, calcExactSolution
"""
Solves 1D Euler equation using Newton iteration, 
only works in absent of shocks
"""
function implicit_solve(degree::Integer, numelem::Integer)
    
    solver, area, q, jac = setup_for_implicit_solve(degree, numelem)

    newton_iterate(solver, area, q, jac)

    postprocessing(solver, area, q)

    return solver, area, q, jac
end

"""
setup system and variables used in implicit solve
"""
function setup_for_implicit_solve(degree::Integer, 
                                  numelem::Integer,
                                  shocks::Bool)
    println("\npreprocessing (constructiong solver, area, q, jac)...")
    # area_left = 2.0
    # area_right = 1.5
    # area_mid = 1.0
    # area_star = 0.8
    # γ = Quasi1DEuler.gamma
    # temp_stag = 300.0
    # press_stag = 100000.0
    vis0 = 10.0/(numelem*degree)
    s0 = -(1.0 + 4.0*log10(degree))

    # set the inlet and outlet flow states (for the boundary conditions)
    rho, rho_u, e = getFlowExact(area_star, 0.0, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(γ*press/rho_ref)
    bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    rho, rho_u, e = getFlowExact(area_star, 1.0, area_right, temp_stag, 
                                 press_stag, subsonic=!shocks, 
                                 area_shock=nozzleArea(xshock))
    bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    #--------------------------------------------------------------------------------
    # construct the solver object
    # degree = 1   # degree of the polynomial basis
    # numelem = 10 # number of elements
    solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out, shocks=shocks, 
                                  vis0=vis0, s0=s0, kappa=kappa)

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
    Tsln = Float64
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
function calcStateJacobian{T,Tarea,Tres}(solver::EulerSolver{T},
                                         area::AbstractArray{Tarea,2},
                                         qfloat::AbstractArray{Tres,3},
                                         jac::SparseMatrixCSC;
                                         nu::Float64=0.0)
    q = convert(Array{Complex128,3}, qfloat)
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
    eps = 1e-25
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

"""
Verify using finite difference if the SparseMatrixCSC 
jacobian matrix is computed correctly
"""
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

"""
Add a dense submatrix to a SparseMatrixCSC matrix
"""
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
Create Jacobian matrix of type of SparseMatrixCSC
"""
function createJac{Tsol}(q::AbstractArray{Tsol, 3})
    n_field = size(q, 1)
    n_dof_per_elem = size(q, 2)
    n_elem = size(q, 3)

    m = Int64(length(q))
    n = m
    blk_size = n_field * n_dof_per_elem
    nnz = blk_size * blk_size * (n_elem * 3 - 2)
    colptr = zeros(Int64, m+1)
    rowval = zeros(Int64, nnz)
    nzval = zeros(Float64, nnz)

    # special treatment for 1st and last element
    colptr[1] = 1
    colptr[2 : blk_size+1] = blk_size * 2
    colptr[m-blk_size+1 : m+1] = blk_size * 2
    for el = 2 : n_elem - 1
        for j = 1 : blk_size
            col = (el - 1) * blk_size + j
            colptr[col+1] = blk_size * 3
        end
    end

    for j = 1 : m
        colptr[j+1] += colptr[j]
    end

    @assert(colptr[m+1] == nnz + 1)

    # special treatment for 1st and last element
    for j = 1 : blk_size    # loop over local columns
        # first element
        row_offset = 0
        col_offset = 0
        col = col_offset + j
        for i = 1 : 2*blk_size # loop over rows
            loc = colptr[col] + i - 1
            rowval[loc] = i + row_offset
        end
        # last element
        row_offset = m - 2*blk_size
        col_offset = m - blk_size
        col = col_offset + j
        for i = 1 : 2*blk_size    # loop over rows 
            loc = colptr[col] + i - 1
            rowval[loc] = i + row_offset
        end
    end

    for el = 2 : n_elem - 1
        col0 = (el - 1) * blk_size + 1 
        col1 = col0 + blk_size - 1
        row0 = col0 - blk_size
        row1 = row0 + 3*blk_size - 1
        for col = col0 : col1
            for i = 1 : blk_size * 3
                rowval[colptr[col] + i - 1] = row0 + i -1
            end
        end
    end
    
    return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

"""
Postprocessing: calculate solution error, J1 and J2
"""
function postprocessing{T,Tarea,Tres}(solver::EulerSolver{T},
                                      area::AbstractArray{Tarea,2},
                                      q::AbstractArray{Tres,3})
    println("\nPostprocessing...")
    # evaluate and print the solution error
    rho, rho_u, e = getFlowExact(area_star, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(γ*press/rho_ref)
    qexact = calcExactSolution(area_star, rho_ref, a_ref, 
                               solver.x, area, subsonic=false)
    q_float = zeros(q, Float64)
    q_float[:] = real(q)
    L2_err, max_err = calcSolutionError(solver, q_float, qexact)
    println("L2 solution error  = ",L2_err)
    println("max solution error = ",max_err)
    # evaluate and print the functional values
    if solver.shocks
        j1ex = -0.38482425832590694
        j2ex = -0.10202215753949603
    else
        j1ex = -0.35194635479522557
        j2ex = -0.11000142657405404
    end
    J1 = calcIntegratedSource(solver, area, q)
    J2 = calcWeightedSource(solver, area, q)
    println("J_1: Integrated Source = ", real(J1))
    println("J_2: Outlet Pressure   = ", real(J2))
    J1_error = abs(real(J1) - j1ex)
    J2_error = abs(real(J2) - j2ex)
    println("J1_error = ", J1_error)
    println("J2_error = ", J2_error)

    return L2_err, max_err, J1_error, J2_error
end

"""
Transpose jac. Given that the sparsity is symmetric,
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
