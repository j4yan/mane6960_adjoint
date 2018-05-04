"""
    sortVec!(solver, u, usort)

Sorts the vector q sequentially by x, so it can be plotted 
"""
function sortVec!{T}(solver::EulerSolver{T}, u::AbstractArray{T,2},
                     usort::AbstractArray{T,2})
  idx = sortperm(vec(solver.x[1,:,1]))
  tmp = zeros(T, (solver.sbp.numnodes))
  for k = 1:solver.numelem
    tmp[:] = u[idx,k]
    usort[:,k] = tmp[:]
  end
end

function sortVec!{T}(solver::EulerSolver{T}, u::AbstractArray{T,3},
                     usort::AbstractArray{T,3})
  idx = sortperm(vec(solver.x[1,:,1]))
  tmp = zeros(T, (size(u,1), solver.sbp.numnodes))
  for k = 1:solver.numelem
    tmp[:,:] = u[:,idx,k]
    usort[:,:,k] = tmp[:,:]
  end
end

"""
    interpSolution!(solver_c, solver_f, u_c, u_f)

Interpolates coarse solution `u_c` onto fine solution `u_f`.  The solver objects
`solver_c` and `solver_f` refer to the coarse- and fine-spaces for p enrichment,
and they must have the same number of elements.
"""
function interpSolution!{T,Tsol}(solver_c::EulerSolver{T},
                                 solver_f::EulerSolver{T},
                                 u_c::AbstractArray{Tsol,3},
                                 u_f::AbstractArray{Tsol,3})
  @assert( size(u_c,2) == solver_c.sbp.numnodes &&
           size(u_c,3) == solver_c.numelem )
  @assert( size(u_f,2) == solver_f.sbp.numnodes &&
           size(u_f,3) == solver_f.numelem )
  @assert( solver_c.numelem == solver_f.numelem )
  x_f = calcnodes(solver_f.sbp, solver_f.sbp.vtx)
  R = SummationByParts.buildinterpolation(solver_c.sbp, x_f)
  for k = 1:solver_c.numelem
    for j = 1:solver_c.sbp.numnodes
      for i = 1:solver_f.sbp.numnodes
        for l = 1:3
          u_f[l,i,k] += R[i,j]*u_c[l,j,k]
        end
      end
    end
  end
end

"""
    calcInnerProd(solver, u, v)

Returns the (approximate) integral inner product between `u` and `v`
"""
function calcInnerProd{T}(solver::EulerSolver{T}, u::AbstractArray{T,3},
                          v::AbstractArray{T,3})
  prod = 0.0
  h = 1.0/solver.numelem
  for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
      fac = solver.sbp.w[i]*h
      for j = 1:3
        prod += u[j,i,k]*v[j,i,k]*fac
      end
    end
  end
  return prod
end

"""
    getSolutionError(solver, q, qexact)

Return the (approximate) L2 and max solution error in the three solution
components, ρ, ρu, and e based on the given exact solution `qexact`
"""
function calcSolutionError{T}(solver::EulerSolver{T}, 
                              q::AbstractArray{Float64,3},
                              qexact::AbstractArray{Float64,3})
  error = zeros(3,solver.sbp.numnodes)
  max_err = zeros(3)
  L2_err = zeros(3)
  h = 1.0/solver.numelem
  for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
      error[1,i] = abs(qexact[1,i,k] - q[1,i,k])
      error[2,i] = abs(qexact[2,i,k] - q[2,i,k])
      error[3,i] = abs(qexact[3,i,k] - q[3,i,k])
    end
    max_err = max.(maximum(error,2), max_err)
    error[:,:] .*= error[:,:]
    volumeIntegrateElement!(solver.sbp, error, error)
    L2_err[:] += h*sum(error,2)
  end
  return sqrt.(L2_err), max_err
end

export nozzleArea, calcExactSolution


include("parameters.jl")
"""
    nozzleArea(x)

Returns the nozzle cross-sectional area at a given location `x`.
"""
function nozzleArea(x::T) where T
  a = area_left
  b = 4.0*area_mid - 5.0*area_left + area_right
  c = -4.0*(area_right - 2.0*area_left + area_mid)
  d = 4.0*(area_right - area_left)
  return a + x*(b + x*(c + x*d))
end

"""
    calcExactSolution(area_star, rho_ref, a_ref, x, area[, subsonic=true,
                      area_shock=nozzleArea(0.7)])

Returns an array containing the exact solution; in contrast with
Quasi1DEuler.getFlowExact, this function nondimensionalizes the state.  Use this
function if you want to compute the solution error.
"""
function calcExactSolution(area_star::Float64, rho_ref::Float64,
                           a_ref::Float64, x::Array{Float64,3},
                           area::Array{Float64,2};
                           subsonic::Bool=true, 
                           area_shock::Float64=nozzleArea(xshock))
  numnodes = size(area,1)
  numelem = size(area,2)
  q = zeros(3,numnodes,numelem)
  for k = 1:numelem
    for i = 1:numnodes
      q[1,i,k], q[2,i,k], q[3,i,k] = getFlowExact(area_star, x[1,i,k], area[i,k],
                                                  temp_stag, press_stag,
                                                  subsonic=subsonic,
                                                  area_shock=area_shock)
      q[1,i,k] /= rho_ref
      q[2,i,k] /= (a_ref*rho_ref)
      q[3,i,k] /= (rho_ref*a_ref*a_ref)
    end
  end
  return q
end

