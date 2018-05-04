
"""
    getMachRelation(area, area_star)

Return a function that provides the Mach relation for the given area ratio
"""
function getMachRelation{T}(area_star::T, area::T)
  function machRelation(M::T)
    return (1.0./M)*(((2.0/(gamma+1.0))*(1.0 + ((gamma-1.0)/2.0)*M*M))^
                     ((gamma+1.0)/(2.0*(gamma-1.0)))) - area/area_star
  end
  return machRelation
end

"""
    getFlowExact(area_star, area, temp_stag, press_stag[, subsonic=true])

Return the exact conservative variables based on the Area-Mach relation
"""
function getFlowExact{T}(area_star::T, area::T, temp_stag::T, press_stag::T,
                         subsonic::Bool=true)
  local Mach::T
  if abs(area_star - area) < eps(T)
    Mach = one(T)
  else
    machRelation = getMachRelation(area_star, area)
    bracket = zeros(2)
    if subsonic
      # subsonic branch
      bracket = [0.001; 1.0]
    else
      # supersonic branch
      bracket = [1.0; 4.0]
    end
    Mach = fzero(machRelation, bracket)
  end
  temp = 1.0./(1.0 + ((gamma - 1.0)/2.0)*(Mach*Mach))
  press = press_stag*(temp^(gamma/(gamma-1.0)))
  temp *= temp_stag
  ρ = press/(R_gas*temp)
  u = Mach*sqrt(gamma*press/ρ)
  ρu = ρ*u
  e = 0.5*ρ*u*u + press/(gamma-1.0)
  return ρ, ρu, e
end

"""
    getFlowExact(area_star, x, area, temp_stag, press_stag[, subsonic=true])

Return the exact conservative variables based on the Area-Mach relation at a set
of points. `area_shock` denotes the nozzle area at which the shock is located.

"""
function getFlowExact{T}(area_star::T, x::AbstractArray{T,1},
                         area::AbstractArray{T,1},
                         temp_stag::T, press_stag::T; subsonic::Bool=true,
                         area_shock::T=area_star)
  num_nodes = length(area)
  q = zeros(3, num_nodes)
  if subsonic
    for i = 1:num_nodes
      ρ, ρu, e = getFlowExact(area_star, area[i], temp_stag, press_stag, 
                              true)
      q[:,i] = [ρ; ρu; e]
    end
  else
    # find the Mach number just upstream of the shock
    ρ, ρu, e = getFlowExact(area_star, area_shock, temp_stag, press_stag, 
                            false)
    press = (gamma-1)*(e - 0.5*ρu*ρu/ρ)
    Ms = (ρu/ρ)/sqrt(gamma*press/ρ)
    # normal shock relations
    Mafter = sqrt((1 + 0.5*(gamma-1)*Ms*Ms)/(gamma*Ms*Ms - (gamma-1)*0.5))
    press_stag_after = press*(((Ms*(gamma+1))^2/(4*gamma*Ms^2 - 2*(gamma-1)))^
                              (gamma/(gamma-1))*
                              (1 - gamma + 2gamma*Ms*Ms)/(gamma+1))
    area_star_after = area_shock/sqrt((1/Mafter^2).*
                                      ((2/(gamma+1)).*(1 + 0.5*(gamma-1).*
                                                       (Mafter.^2))).^
                                      ((gamma+1)/(gamma-1)))
    for i = 1:num_nodes
      local ρ::Float64
      local ρu::Float64
      local e::Float64
      if x[i] <= 0.5
        ρ, ρu, e = getFlowExact(area_star, area[i], temp_stag, press_stag, 
                                true)
      elseif x[i] > 0.5 && area[i] < area_shock
        ρ, ρu, e = getFlowExact(area_star, area[i], temp_stag, press_stag, 
                                false)
      else
        ρ, ρu, e = getFlowExact(area_star_after, area[i], temp_stag,
                                press_stag_after, true)
      end
      q[:,i] = [ρ; ρu; e]
    end
  end
  return q
end

