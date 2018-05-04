# This is a simple test that solves the 1d-Euler flow problem
# To run this script, open Julia and then type include("assign4_example.jl")

include("Quasi1DEuler.jl")
using Quasi1DEuler
# using PyPlot

# set some constants and parameters
γ = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0
src_sigma = 0.05

implicit = false  # if true, use CN implicit time-marching method
implicit = true  # if true, use CN implicit time-marching method
degree = 10       # degree of the polynomial basis
numelem = 10      # number of elements
final_time = 2.0  # evolve to time = final_time
# final_time = 1.0  # evolve to time = final_time
CFL = 1.0         # CFL number sets the step size, but we cannot go much beyond 1.0

function initCondition!{T,Tsol}(x::AbstractArray{T,3}, q::AbstractArray{Tsol,3})
    xc = 0.25
    sig = src_sigma
    fill!(q, 1.0)
    for k = 1:size(x,3)
        for i = 1:size(x,2)
            q[2,i,k] -= 0.01*exp(-((x[1,i,k]-xc)/sig)^2)
        end
    end
end

#--------------------------------------------------------------------------------
# set the inlet and outlet flow states (for the boundary conditions)
bc_in = [1.0; 1.0; 1.0]
bc_out = bc_in

#--------------------------------------------------------------------------------
# construct the solver object
solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out, src_x=0.5, src_sig=0.05)

#--------------------------------------------------------------------------------
# initialize the discrete solution using initCondition function above
q = zeros(3, solver.sbp.numnodes, numelem)
initCondition!(solver.x, q)

#--------------------------------------------------------------------------------
# determine the number of time steps and initialize the source term
N = Quasi1DEuler.calcNumSteps(solver, final_time, cfl=CFL)
src = zeros(N+1)

#--------------------------------------------------------------------------------
# evolve the state using RK4 or Crank-Nicholson
display = true   # if true, display information at each iteration
save = true      # if true, the solution (for each time) is saved to file
filename = "solsave.dat"
if implicit
    unsteadySolveCN!(solver, src, q, final_time, N, display=display, save=save,
                     savefile=filename, tol=1e-8, abstol=1e-10)
else
    unsteadySolveRK4!(solver, src, q, final_time, N, display=display, save=save,
                      savefile=filename)
end

#--------------------------------------------------------------------------------
# evaluate the functional
fun = calcTimeIntegratedObjective!(solver, q, filename, final_time, N)
println("functional value = ",fun)

#--------------------------------------------------------------------------------
# plot the final pressure; the nodes on each element are not ordered
# sequentially, so the coordinates and solution have to be sorted first for
# plotting.
if false
    # NOTE: if you want to see a non-constant pressure field, you should set
    # final_time < 1.0
    press = zeros(solver.sbp.numnodes, solver.numelem)
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            press[i,k] = Quasi1DEuler.calcPressure(view(q,:,i,k))
        end
    end
    x_s = zeros(solver.x)  
    sortVec!(solver, solver.x, x_s)
    press_s = zeros(press)
    sortVec!(solver, press, press_s)
    # PyPlot.plot(vec(x_s), vec(press_s), "-k")
end

#-----------------------------------------------------------
# below this line is added by Jianfeng Yan
#-----------------------------------------------------------

#----------------------------------------------------------------------------------
psi_all = unsteadyAdjointSolveCN(solver, final_time, N, savefile=filename)
fout = open("psi_x=1.dat", "w")
right_end_idx = 1
xmax = solver.x[1,1,1]
for i = 2 : size(solver.x, 2)
    if xmax < solver.x[1, i, 1]
        right_end_idx = i
    end
end
psi_x1  = view(psi_all, :, right_end_idx, size(psi_all, 3), :)

for i = 0 : N
    println(fout, i*final_time/N, " ", psi_x1[1, i+1], " ", psi_x1[2, i+1], " ", psi_x1[3, i+1])
end
close(fout)
# read all primal solution
q_all = zeros(size(psi_all))
nsteps = size(q_all, 4) - 1
fsave = open(filename, "r")
for i = 0 : nsteps
    q = zeros(Float64, 3, solver.sbp.numnodes, solver.numelem)
    read!(fsave, q)
    q_all[:,:,:, i+1] = q[:,:,:]
end

for i = 0 : N
    time = i * final_time / N
    filename = ".dat"
    if i == 0
        filename = "psi_0.0.dat"
    elseif i == N
        filename = "psi_2.0.dat"
    elseif abs(time - 0.25) < 1e-8
        filename = "psi_0.25.dat"
    elseif abs(time - 0.5) < 1e-8
        filename = "psi_0.5.dat"
    else
        continue
    end
    fout = open(filename, "w")
    time_step = i+1
    psi = view(psi_all, :, :, :, time_step)
    psi_rho = view(psi, 1, :, :)
    psi_rhou = view(psi, 2, :, :)
    psi_rhoE = view(psi, 3, :, :)

    # the coordinates and solution have to be sorted first for plotting.
    x_s = zeros(solver.x)  
    sortVec!(solver, solver.x, x_s)
    psi_rho_s = zeros(psi_rho)
    sortVec!(solver, psi_rho, psi_rho_s)
    psi_rhou_s = zeros(psi_rhou)
    sortVec!(solver, psi_rhou, psi_rhou_s)
    psi_rhoE_s = zeros(psi_rhoE)
    sortVec!(solver, psi_rhoE, psi_rhoE_s)

    for j = 1 : length(x_s)
        println(fout, x_s[j], " ", psi_rho_s[j], " ", psi_rhou_s[j], " ", psi_rhoE_s[j])
    end
    close(fout)
end

#------------------------------------------------------------------------
# compute total derivative of the functional with respect to the source,
# using adjoint-based approach
dt = final_time / N
DJDs = cmptTotalDerivOfFunc(solver, dt, src, q_all, psi_all)
println("outputting DJDs...")
fout = open("DJDs_adjoint.dat", "w")
for i = 0 : N
    println(fout, "i = $i, ", DJDs[i+1])
end
close(fout)

#------------------------------------------------------------------
# using FD to compute dfunc_ds
FD_verification = true
FD_verification = false
if FD_verification
    println("compute DJDs using FD...")
    eps_pert = 1e-3
    fout = open("DJDs_FD.dat", "w")
    dfunc_ds = zeros(size(src))

    # for i = 0 : div(N, 50)
    for i = 200 : 205
        src[i+1] += eps_pert 
        q = zeros(3, solver.sbp.numnodes, numelem)
        initCondition!(solver.x, q)
        if implicit
            unsteadySolveCN!(solver, src, q, final_time, N, display=false, save=save,
                             savefile=filename, tol=1e-8, abstol=1e-10)
        else
            unsteadySolveRK4!(solver, src, q, final_time, N, display=false, save=save,
                              savefile=filename)
        end

        fun_purt = calcTimeIntegratedObjective!(solver, q, filename, final_time, N)
        src[i+1] -= eps_pert
        dfunc_ds[i+1] = (fun_purt - fun) / eps_pert
        println(fout, "i = $i, ", dfunc_ds[i+1], ", ", DJDs[i+1], ", ", DJDs[i+1] / dfunc_ds[i+1])
    end
    close(fout)
end
