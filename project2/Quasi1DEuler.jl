module Quasi1DEuler

using SummationByParts
using Roots

import Base.isless

# define some constants needed for the PDE
const gamma = 1.4
const gami = gamma - 1.0
const R_gas = 287

export EulerSolver
export getFlowExact
export calcWeakResidual!
export solveRK4!
export sortVec!
export calcOutletPressure
export calcIntegratedSource
export calcSolutionError

# the following functions are needed for the complex-step
function cabs(x::Number)
  if real(x) >= 0.0
    return x
  else
    return -x
  end
end
function isless(a::Number, b::Number)
  if real(a) < real(b)
    return true
  else
    return false
  end
end

"""
    EulerSolver

Stores information for a quasi-1D-Euler solver
"""
type EulerSolver{T<:Number}
  """number of elements"""
  numelem::Int
  """summation-by-parts operator type"""
  sbp::LineSegSBP{T}
  """summation-by-parts face type"""
  sbpface::LineSegFace{T}
  """boundary values of (rho, rho*u, e) at the inlet of the nozzle"""
  bc_in::Array{T,1}
  """boundary values of (rho, rho*u, e) at the outlet of the nozzle"""
  bc_out::Array{T,1}
  """node locations for each element"""
  x::Array{T,3}
  """determinant of the mapping Jacobian"""
  jac::Array{T,2}
  """list of boundary faces"""
  bfaces::Array{Boundary}
  """list of element interfaces"""
  ifaces::Array{Interface} 
  
  """inner constructor"""
  function EulerSolver{T}(degree::Int, numelem::Int,
                          bc_in::AbstractArray{T,1},
                          bc_out::AbstractArray{T,1}) where {T<:Number}
    @assert( degree > 0 )
    @assert( numelem > 0 )
    # define the operators needed for FEM operations
    sbp = getLineSegSBPLobbato(degree=degree, Tsbp=T)
    sbpface = getLineSegFace(degree, sbp.cub, sbp.vtx)
    x = zeros(T, (1, sbp.numnodes, numelem) )
    jac = zeros(T, (sbp.numnodes, numelem))
    # a uniform mesh is used
    h = 1.0/numelem
    for k = 1:numelem
      vtx = reshape([(k-1)*h; k*h], (2,1))
      x[1,:,k] = calcnodes(sbp, vtx)
      jac[:,k] = 0.5*h
    end
    # set-up the element interfaces and boundaries
    ifaces = Array{Interface}(numelem-1)
    for k = 1:numelem-1
      ifaces[k] = Interface(k,k+1,2,1,1)
    end
    bfaces = Array{Boundary}(2)
    bfaces[1] = Boundary(1,1)
    bfaces[2] = Boundary(numelem,2)    
    new(numelem, sbp, sbpface, bc_in, bc_out, x, jac, bfaces, ifaces)
  end
end

include("utils.jl")
include("exact_nozzle.jl")
include("residual.jl")
include("solve_explicit.jl")
include("functionals.jl")

end
