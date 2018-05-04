# methods related to computing the residual

"""
    calcPressure(q)

Return the pressure based on the given state
"""
function calcPressure{Tsol}(q::AbstractArray{Tsol,1})
  return (gamma-1.)*(q[3] - 0.5*(q[2]*q[2])/q[1])
end

"""
    calcEntropy(q)

Returns the entropy based on the given state
"""
function calcEntropy{Tsol}(q::AbstractArray{Tsol,1})
  p = calcPressure(q)
  rho = q[1]
  return -log(p/rho^gamma)
end

"""
    calcEulerFlux!(area, q, flux)

Compute the quasi-1D-Euler flux and store in `flux`
"""
function calcEulerFlux!{Tarea,Tsol,Tflux}(area::Tarea, q::AbstractArray{Tsol,1},
                                          flux::AbstractArray{Tflux,1})
  rho = q[1]
  u = q[2]/rho
  e = q[3]
  press = calcPressure(q)
  flux[1] = rho*u*area
  flux[2] = (rho*u*u + press)*area
  flux[3] = u*(e + press)*area
end

"""
    calcRoeFlux!(area, q, flux)

Compute the Roe numerical flux based on the "left" and "right" states and areas,
and store in `flux`
"""
function calcRoeFlux!{Tarea,Tsol,Tflux}(areaL::Tarea, areaR::Tarea,
                                        qL::AbstractArray{Tsol,1},
                                        qR::AbstractArray{Tsol,1},
                                        flux::AbstractArray{Tflux,1})
  # Declaring constants 
  d1_0 = one(Tflux)
  d0_0 = zero(Tflux)
  d0_5 = 0.5*one(Tflux)
  sat_Vn = 0.025*one(Tflux)
  sat_Vl = 0.025*one(Tflux)

  # find the Roe-average state
  dA = 0.5*(areaL + areaR)
  
  fac = d1_0/qL[1]
  uL = qL[2]*fac
  phi = d0_5*uL*uL
  HL = gamma*qL[3]*fac - gami*phi
  
  fac = d1_0/qR[1]
  uR = qR[2]*fac
  phi = d0_5*uR*uR
  HR = gamma*qR[3]*fac - gami*phi

  sqL = sqrt(qL[1]); sqR = sqrt(qR[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*u*u
  a = sqrt(gami*(H - phi))
  Un = u*dA

  # compute the wave speeds
  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = cabs(Un) + dA*a
  lambda1 = d0_5*(max(cabs(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = d0_5*(max(cabs(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = d0_5*(max(cabs(lambda3),sat_Vl *rhoA) - lambda3)

  dq1 = qL[1] - qR[1] 
  dq2 = qL[2] - qR[2]
  dq3 = qL[3] - qR[3]

  calcEulerFlux!(dA, qL, flux)

  # diagonal matrix multiply
  flux[1] += lambda3*dq1
  flux[2] += lambda3*dq2
  flux[3] += lambda3*dq3

  # get E1*dq
  E1dq = zeros(Tflux, 3)
  E1dq[1] = phi*dq1 - u*dq2 + dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*H
  
  # get E2*dq
  E2dq = zeros(Tflux, 3)
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + dq2*dA
  E2dq[3] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*dA  

  # add to flux
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  for n = 1:3
    flux[n] += tmp1*(tmp2*E1dq[n] + tmp3*E2dq[n])
  end  
  
  # get E3*dq
  E1dq[1] = -Un*dq1 + dA*dq2
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*H

  # get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 + dq3
  E2dq[3] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*dA

  # add to flux
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  for n = 1:3
    flux[n] += tmp1*(E1dq[n] + gami*E2dq[n])
  end
end

"""
    calcWeakResidual!(solver)

Return the weak-form of the quasi-1d-Euler residual.
"""
function calcWeakResidual!{T,Tarea,Tsol,Tres}(solver::EulerSolver{T},
                                              area::AbstractArray{Tarea,2},
                                              q::AbstractArray{Tsol,3},
                                              res::AbstractArray{Tres,3};
                                              nu::T=0.0)
  @assert( size(q,1) == size(res,1) == 3 )
  @assert( size(q,2) == size(res,2) == size(area,1) == solver.sbp.numnodes )
  @assert( size(q,3) == size(res,3) == size(area,2) == solver.numelem )
  fill!(res, zero(Tres))

  # add the volume contributions to the residual
  flux = zeros(Tres, (3,solver.sbp.numnodes) )
  dAdx = zeros(Tarea, (solver.sbp.numnodes) )
  for k = 1:solver.numelem
    # add the volume contribution -F*dv/dx
    for i = 1:solver.sbp.numnodes
      calcEulerFlux!(area[i,k], view(q,:,i,k), view(flux,:,i))
    end
    weakDifferentiateElement!(solver.sbp, 1, flux, view(res,:,:,k),
                              SummationByParts.Subtract(), true)  
    # add the volume contribution p*dA/dx to the momentum
    fill!(dAdx, zero(Tarea))
    differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
    fill!(flux, zero(Tres))
    for i = 1:solver.sbp.numnodes
      flux[2,i] = -calcPressure(view(q,:,i,k))*dAdx[i]
    end
    volumeIntegrateElement!(solver.sbp, flux, view(res,:,:,k))
  end

  # add the interior-face penalty terms
  qfaceL = zeros(Tsol, (3,1))
  qfaceR = zeros(qfaceL)
  areaL = zeros(Tarea, (1))
  areaR = zeros(Tarea, (1))
  flux_face = zeros(Tres, (3,1))
  for (findex, face) in enumerate(solver.ifaces)
    # interpolate the solution to the face
    interiorFaceInterpolate!(solver.sbpface, face,
                             view(q,:,:,face.elementL),
                             view(q,:,:,face.elementR),
                             qfaceL, qfaceR)
    # interpolate the area to the face
    interiorFaceInterpolate!(solver.sbpface, face,
                             view(area,:,face.elementL),
                             view(area,:,face.elementR),
                             areaL, areaR)
    calcRoeFlux!(areaL[1], areaR[1], view(qfaceL,:,1), view(qfaceR,:,1),
                 view(flux_face,:,1))
    interiorFaceIntegrate!(solver.sbpface, face, flux_face,
                           view(res,:,:,face.elementL),
                           view(res,:,:,face.elementR))
  end

  # add the boundary-condition penalties; here, the internal state is stored in
  # qfaceR and the boundary state is stored in qfaceL, even if the state is
  # actually the "left" state.
  for (findex, face) in enumerate(solver.bfaces)
    # interpolate the solution to the face
    boundaryFaceInterpolate!(solver.sbpface, face.face,
                             view(q,:,:,face.element), qfaceR)
    # interpolate the area to the face
    boundaryFaceInterpolate!(solver.sbpface, face.face,
                             view(area,:,face.element), areaR)
    if solver.sbpface.normal[face.face] < 0.0
      # the element is at the left-side (inlet) boundary
      qfaceL[:,1] = solver.bc_in[:]
      calcRoeFlux!(areaR[1], areaR[1], view(qfaceL,:,1), view(qfaceR,:,1),
                   view(flux_face,:,1))
      scale!(flux_face, -1.0)
    else
      # the element is at the right-side (outlet) boundary
      qfaceL[:,1] = solver.bc_out[:]
      calcRoeFlux!(areaR[1], areaR[1], view(qfaceR,:,1), view(qfaceL,:,1),
                   view(flux_face,:,1))
    end
    boundaryFaceIntegrate!(solver.sbpface, face.face, flux_face,
                           view(res,:,:,face.element))
  end

  if solver.shocks || nu > 0.0
    # add artificial viscosity if necessary
    scale!(res, (1.0-nu))
    addLaplacianArtificialViscosity!(solver, area, q, res, nu=nu)
  end
end

"""
    calcWeakResidual!(solver, src, q, res)

Return the weak-form of the 1D-Euler residual for unsteady flows with a source.
"""
function calcWeakResidual!{T,Tsrc,Tsol,Tres}(solver::EulerSolver{T}, src::Tsrc,
                                             q::AbstractArray{Tsol,3},
                                             res::AbstractArray{Tres,3};
                                             zerofill::Bool=true)
  @assert( size(q,1) == size(res,1) == 3 )
  @assert( size(q,2) == size(res,2) == solver.sbp.numnodes )
  @assert( size(q,3) == size(res,3) == solver.numelem )
  if zerofill
    fill!(res, zero(Tres))
  end

  # add the volume contributions to the residual
  flux = zeros(Tres, (3,solver.sbp.numnodes) )
  for k = 1:solver.numelem
    # add the volume contribution -F*dv/dx
    for i = 1:solver.sbp.numnodes
      calcEulerFlux!(1.0, view(q,:,i,k), view(flux,:,i))
    end
    weakDifferentiateElement!(solver.sbp, 1, flux, view(res,:,:,k),
                              SummationByParts.Subtract(), true)
  end

  # add the interior-face penalty terms
  qfaceL = zeros(Tsol, (3,1))
  qfaceR = zeros(qfaceL)
  flux_face = zeros(Tres, (3,1))
  for (findex, face) in enumerate(solver.ifaces)
    # interpolate the solution to the face
    interiorFaceInterpolate!(solver.sbpface, face,
                             view(q,:,:,face.elementL),
                             view(q,:,:,face.elementR),
                             qfaceL, qfaceR)
    calcRoeFlux!(1.0, 1.0, view(qfaceL,:,1), view(qfaceR,:,1),
                 view(flux_face,:,1))
    interiorFaceIntegrate!(solver.sbpface, face, flux_face,
                           view(res,:,:,face.elementL),
                           view(res,:,:,face.elementR))
  end

  # add the boundary-condition penalties; here, the internal state is stored in
  # qfaceR and the boundary state is stored in qfaceL, even if the state is
  # actually the "left" state.
  for (findex, face) in enumerate(solver.bfaces)
    # interpolate the solution to the face
    boundaryFaceInterpolate!(solver.sbpface, face.face,
                             view(q,:,:,face.element), qfaceR)
    if solver.sbpface.normal[face.face] < 0.0
      # the element is at the left-side (inlet) boundary
      qfaceL[:,1] = solver.bc_in[:]
      calcRoeFlux!(1.0, 1.0, view(qfaceL,:,1), view(qfaceR,:,1),
                   view(flux_face,:,1))
      scale!(flux_face, -1.0)
    else
      # the element is at the right-side (outlet) boundary
      qfaceL[:,1] = solver.bc_out[:]
      calcRoeFlux!(1.0, 1.0, view(qfaceR,:,1), view(qfaceL,:,1),
                   view(flux_face,:,1))
    end
    boundaryFaceIntegrate!(solver.sbpface, face.face, flux_face,
                           view(res,:,:,face.element))
  end

  # add the source term
  src_elem = zeros(Tsrc, (3,solver.sbp.numnodes) )
  for k = 1:solver.numelem    
    for i = 1:solver.sbp.numnodes
      dx = solver.x[1,i,k] - solver.src_x
      src_elem[2,i] = src*exp(-dx*dx/solver.src_sig2)*solver.jac[i,k]
    end
    volumeIntegrateElement!(solver.sbp, src_elem, view(res,:,:,k),
                            SummationByParts.Subtract())
  end  
end
