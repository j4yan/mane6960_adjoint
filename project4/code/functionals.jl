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

"""
    calcTimeIntegratedObjective!(solver, q, savefile, time, nsteps[, sensor_x,
                                 sensor_sig, press_targ])

A time- and space-integrated functional; the pressure difference, `press(q) -
press_targ`, is squared and weighted by a squared-exponential that is centered
at `sensor_x` and has a width of `sensor_sig`.  The data required to integrate
the solution is in `savefile`, which must be consistent in dimension with `q`
and the number of time steps.
"""
function calcTimeIntegratedObjective!{T,Tsol}(solver::EulerSolver{T},
                                              q::AbstractArray{Tsol,3},
                                              savefile::AbstractString,
                                              time::T, nsteps::Int;
                                              sensor_x::T=0.85, sensor_sig::T=0.05,
                                              press_targ::T=0.2)
  dt = time/nsteps
  # open the solution file
  fsol = open(savefile, "r")
  integrand = zeros(solver.sbp.numnodes)
  scaled_int = zeros(integrand)
  fun = 0.0
  # loop over the time steps
  for step = 0:nsteps
    # read solution
    read!(fsol, q)
    fun_step = 0.0
    for k = 1:solver.numelem
      for i = 1:solver.sbp.numnodes
        dx = (solver.x[1,i,k] - sensor_x)/sensor_sig
        kern = exp(-dx*dx)
        dpress = calcPressure(view(q,:,i,k)) - press_targ
        integrand[i] = 0.5*dpress*dpress*kern*solver.jac[i,k]
      end
      fill!(scaled_int, zero(Tsol))
      volumeIntegrateElement!(solver.sbp, integrand, scaled_int)
      fun_step += sum(scaled_int)
    end
    # get trapezoid quadrature weight and then add contribution
    fac = dt
    if step == 0 || step == nsteps
      fac *= 0.5
    end
    fun += fac*fun_step
  end
  close(fsol)
  return fun  
end
