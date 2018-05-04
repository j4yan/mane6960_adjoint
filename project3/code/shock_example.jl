# This is a simple test that solves the "standard" quasi-1d nozzle flow problem
# To run this script, open Julia and then type include("example.jl")

include("Quasi1DEuler.jl")
include("parameters.jl")
using Quasi1DEuler
# using PyPlot

shocks = true
if area_star < 1.0 - 1e-8
    shocks = false
end

solver, area, q, jac = setup_for_implicit_solve(degree, numelem, shocks)
if shocks
    solveHomotopy!(solver, area, q, jac, maxouter=maxouter, maxinner=10, 
                   tol=tol, inner_tol=1e-4, alpha_max=0.02, 
                   phi_targ=5.0, display=display)
else
    solveNewton!(solver, area, q, jac, numiter=maxouter, tol=tol)
end

postprocessing(solver, area, q)
exit()
