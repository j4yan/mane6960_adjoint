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
