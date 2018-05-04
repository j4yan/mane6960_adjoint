
"""
Parameters are moved here!
"""
#
# parameters : geometric and gas properties
#
area_left = 2.0
area_right = 1.5
area_mid = 1.0
area_star = 1.0
γ = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0
const xshock = 0.7

#
# parameters: discretization and iteration 
#
degree = 2
numelem = 40
maxouter = 1000
tol = 1e-12
kappa = 0.5
display = true

