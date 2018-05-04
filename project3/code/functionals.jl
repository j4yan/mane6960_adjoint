# methods related to computing various functionals

"""
    calcOutletPressure(solver, area, q)

Returns the pressure at the outlet evaluated using the numerical solution.
"""
function calcOutletPressure{T,Tarea,Tsol}(solver::EulerSolver{T},
                                          area::AbstractArray{Tarea,2},
                                          q::AbstractArray{Tsol,3})
  J = calcPressure(q[:,2,solver.numelem])
  return J
end

"""
    calcIntegratedSource(solver, area, q)

Returns the integrated momentum source.
"""
function calcIntegratedSource{T,Tarea,Tsol}(solver::EulerSolver{T},
                                            area::AbstractArray{Tarea,2},
                                            q::AbstractArray{Tsol,3})
  J = zero(Tarea)*zero(Tsol)*zero(T)
  dAdx = zeros(Tarea, (solver.sbp.numnodes) )
  pdAdx = zeros(typeof(J), (solver.sbp.numnodes))
  HpdAdx = zeros(pdAdx)
  for k = 1:solver.numelem
    # get the derivative of area with respect to x
    fill!(dAdx, zero(Tarea))
    differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
    for i = 1:solver.sbp.numnodes
      # scale by pressure
      pdAdx[i] = calcPressure(view(q,:,i,k))*dAdx[i]
    end
    fill!(HpdAdx, zero(typeof(J)))
    # scale by quadrature weights and sum
    volumeIntegrateElement!(solver.sbp, pdAdx, HpdAdx)
    J += sum(HpdAdx)
  end
  return J
end

"""
    calcWeightededSource(solver, area, q)

Similar to `calcIntegratedSource`, but weights the source by a localized
Guassian bump kernel.  The bump is centered at `xloc` and has a width
deteremined by `sigma`.
"""
function calcWeightedSource{T,Tarea,Tsol}(solver::EulerSolver{T},
                                          area::AbstractArray{Tarea,2},
                                          q::AbstractArray{Tsol,3};
                                          xloc::Float64=0.25,
                                          sigma::Float64=0.05)
  J = zero(Tarea)*zero(Tsol)*zero(T)
  dAdx = zeros(Tarea, (solver.sbp.numnodes) )
  pdAdx = zeros(typeof(J), (solver.sbp.numnodes))
  HpdAdx = zeros(pdAdx)
  for k = 1:solver.numelem
    # get the derivative of area with respect to x
    fill!(dAdx, zero(Tarea))
    differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
    for i = 1:solver.sbp.numnodes
      # scale by pressure and kernel
      x = solver.x[1,i,k]
      fac = exp(-((x-xloc)./(sigma))^2)
      pdAdx[i] = fac*calcPressure(view(q,:,i,k))*dAdx[i]      
    end
    fill!(HpdAdx, zero(typeof(J)))
    # scale by quadrature weights and sum
    volumeIntegrateElement!(solver.sbp, pdAdx, HpdAdx)
    J += sum(HpdAdx)
  end
  return J
end
