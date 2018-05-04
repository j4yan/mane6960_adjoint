# methods related to solving the quasi-1D Euler equations explicitly

"""
    applyInverseMass!(solver, res)

Multiplies 'res' by the inverse mass matrix
"""
function applyInverseMass!{T,Tres}(solver::EulerSolver{T},
                                   res::AbstractArray{Tres,3})
  h = 1.0/solver.numelem
  for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
      fac = 1.0./(solver.sbp.w[i]*solver.jac[i,k])
      for j = 1:3
        res[j,i,k] *= fac
      end
    end
  end  
end

"""
    stepRK4(solver, area, q, dt)

Apply one step of RK4 to a vector based on the problem defined in solver
"""
function stepRK4!{T,Tarea,Tsol}(solver::EulerSolver{T},
                                area::AbstractArray{Tarea,2},
                                q::AbstractArray{Tsol,3}, dt::T)
  @assert(dt > 0.0)

  # 1st stage
  k1 = zeros(q)
  calcWeakResidual!(solver, area, q, k1)
  applyInverseMass!(solver, k1)
  scale!(k1, -1.0)
  q1 = deepcopy(q)
  q1 += (0.5*dt).*k1

  # 2nd stage
  k2 = zeros(q)
  calcWeakResidual!(solver, area, q1, k2)
  applyInverseMass!(solver, k2)
  scale!(k2, -1.0)
  q2 = deepcopy(q)
  q2 += (0.5*dt).*k2

  # 3rd stage
  k3 = zeros(q)
  calcWeakResidual!(solver, area, q2, k3)
  applyInverseMass!(solver, k3)
  scale!(k3, -1.0)
  q3 = deepcopy(q)
  q3 += dt.*k3

  # 4th stage
  k4 = zeros(q)
  calcWeakResidual!(solver, area, q3, k4)
  applyInverseMass!(solver, k4)
  scale!(k4, -1.0)

  # update solution and compute norm of residual (at beginning of step)
  res_norm = 0.0
  for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
      fac = solver.sbp.w[i]*solver.jac[i,k]
      for j = 1:3
        q[j,i,k] += (dt/6.0)*(k1[j,i,k] + 2.0*k2[j,i,k]
                              + 2.0*k3[j,i,k] + k4[j,i,k])
        res_norm += k1[j,i,k]*fac*k1[j,i,k]
      end
    end
  end
  return sqrt(res_norm)
end

"""
    solveRK4(solver, area, q)

Solve the quasi-1D Euler equations defined by `solver` by time-marching with the
classical 4th-order Runge-Kutta method.
"""
function solveRK4!{T,Tarea,Tsol}(solver::EulerSolver{T},
                                 area::AbstractArray{Tarea,2},
                                 q::AbstractArray{Tsol,3};
                                 cfl::T=1.0,
                                 maxiter::Int=10000,
                                 tol::Float64=1e-10,
                                 display::Bool=false)
  # use left boundary state to estimate the maximum allowable time step
  h = 0.5*minimum(solver.sbp.w)/solver.numelem
  p = Quasi1DEuler.calcPressure(solver.bc_in)
  a = sqrt(gamma*p/solver.bc_in[1])
  u = solver.bc_in[2]/solver.bc_in[1]
  dt = cfl*h/(abs(u)+a)
  if display
    println("Using time-step size dt = ",dt)
  end
  # time march until iterations are exceeded or residual norm is satisfied
  local res_norm0
  for n = 0:maxiter
    res_norm = stepRK4!(solver, area, q, dt)
    if n == 0
      res_norm0 = res_norm
      if display
        println("initial residual norm = ",res_norm0)
      end
    else
      if display
        println("iter ",n,": relative residual norm = ",res_norm/res_norm0)
      end
      if res_norm < tol*res_norm0
        break
      end
    end
  end    
end
