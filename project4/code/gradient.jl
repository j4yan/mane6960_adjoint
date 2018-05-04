# methods related to compute the derivative with respect to parameters

"""
    calcOutletPressuredJda!(solver, area, q, dJda)

Returns the derivative of the outlet pressure with respect to area (its zero).
"""
function calcOutletPressuredJda!{T,Tarea,Tsol,Tfun}(solver::EulerSolver{T},
                                                    area::AbstractArray{Tarea,2},
                                                    q::AbstractArray{Tsol,3},
                                                    dJda::AbstractArray{Tfun,2})
  # the outlet pressure functional has no explicit dependence on the area
  fill!(dJda, zero(Tfun))
end

"""
    calcIntegratedSourcedJda(solver, area, q)

Returns the derivative of the integrated source with respect to area
"""
function calcIntegratedSourcedJda!{T,Tarea,Tsol,Tfun}(solver::EulerSolver{T},
                                                      area::AbstractArray{Tarea,2},
                                                      q::AbstractArray{Tsol,3},
                                                      dJda::AbstractArray{Tfun,2})
  fill!(dJda, zero(Tfun))
  J_bar = one(Tfun)  
  dAdx_bar = zeros(Tarea, (solver.sbp.numnodes) )
  pdAdx_bar = zeros(Tfun, (solver.sbp.numnodes))
  HpdAdx_bar = zeros(pdAdx_bar)
  for k = 1:solver.numelem
    # J += sum(HpdAdx)    
    HpdAdx_bar[:] = J_bar
    # volumeIntegrateElement!(solver.sbp, pdAdx, HpdAdx)
    fill!(pdAdx_bar, zero(Tfun))
    volumeIntegrateElement_rev!(solver.sbp, pdAdx_bar, HpdAdx_bar)
    fill!(dAdx_bar, zero(Tarea))
    for i = 1:solver.sbp.numnodes
      # pdAdx[i] = calcPressure(view(q,:,i,k))*dAdx[i]
      dAdx_bar[i] += calcPressure(view(q,:,i,k))*pdAdx_bar[i]
    end
    # differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
    differentiateElement_rev!(solver.sbp, 1, view(dJda,:,k), dAdx_bar)
  end
end

"""
    calcWeakResidual_rev!(solver)

Return the reverse-mode derivative of the residual with respect to area.
"""
function calcWeakResidual_rev!{T,Tarea,Tsol,Tres}(solver::EulerSolver{T},
                                                  area_bar::AbstractArray{Tarea},
                                                  q::AbstractArray{Tsol,3},
                                                  res_bar::AbstractArray{Tres,3})
  @assert( size(q,1) == size(res_bar,1) == 3 )
  @assert( size(q,2) == size(res_bar,2) == size(area_bar,1)
           == solver.sbp.numnodes )
  @assert( size(q,3) == size(res_bar,3) == size(area_bar,2) == solver.numelem )
  fill!(area_bar, zero(Tres))
  
  # add the volume contributions to the derivative
  flux = zeros(Tsol, (3))
  flux_bar = zeros(Tres, (3,solver.sbp.numnodes) )
  dAdx_bar = zeros(Tres, (solver.sbp.numnodes) )
  for k = 1:solver.numelem
    # volumeIntegrateElement!(solver.sbp, flux, view(res,:,:,k))
    fill!(flux_bar, zero(Tres))
    volumeIntegrateElement_rev!(solver.sbp, flux_bar, view(res_bar,:,:,k))
    fill!(dAdx_bar, zero(Tres))
    for i = 1:solver.sbp.numnodes
      # flux[2,i] = -calcPressure(view(q,:,i,k))*dAdx[i]
      dAdx_bar[i] -= calcPressure(view(q,:,i,k))*flux_bar[2,i]
    end
    # differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
    differentiateElement_rev!(solver.sbp, 1, view(area_bar,:,k), dAdx_bar)

    # weakDifferentiateElement!(solver.sbp, 1, flux, view(res,:,:,k),
    #                           SummationByParts.Subtract(), true)
    fill!(flux_bar, zero(Tres))
    weakDifferentiateElement_rev!(solver.sbp, 1, flux_bar, view(res_bar,:,:,k),
                                  SummationByParts.Subtract(), true)
    for i = 1:solver.sbp.numnodes
      # Note: Euler flux depends linearly on area
      # calcEulerFlux!(area[i,k], view(q,:,i,k), view(flux,:,i))
      calcEulerFlux!(1.0, view(q,:,i,k), flux)
      area_bar[i,k] += dot(vec(flux_bar[:,i]), flux)    
    end
  end

  # add interior face contributions
  qfaceL = zeros(Tsol, (3,1))
  qfaceR = zeros(qfaceL)
  areaL_bar = zeros(Tarea, (1))
  areaR_bar = zeros(Tarea, (1))
  flux_face = zeros(Tres, (3,1))
  flux_face_bar = zeros(Tres, (3,1))
  for (findex, face) in enumerate(solver.ifaces)
    # interpolate the solution to the face; this is needed for reverse sweep
    interiorFaceInterpolate!(solver.sbpface, face,
                             view(q,:,:,face.elementL),
                             view(q,:,:,face.elementR),
                             qfaceL, qfaceR)
    
    # interiorFaceIntegrate!(solver.sbpface, face, flux_face,
    #                       view(res,:,:,face.elementL),
    #                       view(res,:,:,face.elementR))
    fill!(flux_face_bar, zero(Tres))
    interiorFaceIntegrate_rev!(solver.sbpface, face, flux_face_bar,
                               view(res_bar,:,:,face.elementL),
                               view(res_bar,:,:,face.elementR))
    # As with the Euler flux, the Roe flux depends linearly on area
    # calcRoeFlux!(areaL[1], areaR[1], view(qfaceL,:,1), view(qfaceR,:,1),
    #              view(flux_face,:,1))
    calcRoeFlux!(1.0, 1.0, view(qfaceL,:,1), view(qfaceR,:,1),
                 view(flux_face,:,1))
    areaL_bar[1] = 0.5*dot(vec(flux_face),vec(flux_face_bar))
    areaR_bar[1] = areaL_bar[1]    
    # interiorFaceInterpolate!(solver.sbpface, face,
    #                          view(area,:,face.elementL),
    #                          view(area,:,face.elementR),
    #                          areaL, areaR)
    interiorFaceInterpolate_rev!(solver.sbpface, face,
                                 view(area_bar,:,face.elementL),
                                 view(area_bar,:,face.elementR),
                                 areaL_bar, areaR_bar)
  end

  # add boundary face contributions
  for (findex, face) in enumerate(solver.bfaces)
    # interpolate the solution to the face; needed for reverse sweep
    boundaryFaceInterpolate!(solver.sbpface, face.face,
                             view(q,:,:,face.element), qfaceR)
    
    # boundaryFaceIntegrate!(solver.sbpface, face.face, flux_face,
    #                        view(res,:,:,face.element))
    fill!(flux_face_bar, zero(Tres))
    boundaryFaceIntegrate_rev!(solver.sbpface, face.face, flux_face_bar,
                               view(res_bar,:,:,face.element))
    if solver.sbpface.normal[face.face] < 0.0
      # the element is at the left-side (inlet) boundary      
      qfaceL[:,1] = solver.bc_in[:]
      # scale!(flux_face, -1.0)      
      scale!(flux_face_bar, -1.0)
      # calcRoeFlux!(areaR[1], areaR[1], view(qfaceL,:,1), view(qfaceR,:,1),
      #              view(flux_face,:,1))
      calcRoeFlux!(1.0, 1.0, view(qfaceL,:,1), view(qfaceR,:,1),
                   view(flux_face,:,1))
      areaR_bar[1] = dot(vec(flux_face),vec(flux_face_bar))
    else
      # the element is at the right-side (outlet) boundary
      qfaceL[:,1] = solver.bc_out[:]
      # calcRoeFlux!(areaR[1], areaR[1], view(qfaceR,:,1), view(qfaceL,:,1),
      #              view(flux_face,:,1))
      calcRoeFlux!(1.0, 1.0, view(qfaceR,:,1), view(qfaceL,:,1),
                   view(flux_face,:,1))
      areaR_bar[1] = dot(vec(flux_face),vec(flux_face_bar))
    end
    # boundaryFaceInterpolate!(solver.sbpface, face.face,
    #                          view(area,:,face.element), areaR)
    boundaryFaceInterpolate_rev!(solver.sbpface, face.face,
                                 view(area_bar,:,face.element), areaR_bar)
  end
end

