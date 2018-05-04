include("solve_implicit.jl")
using Quasi1DEuler

if !isdir("results")
    mkdir("results")
end
if !isdir("results/adjoints")
    mkdir("results/adjoints")
end
numelems = [100]
# J1
functional = calcIntegratedSource
# J2
# functional = calcOutletPressure
for p = 1 : 4
    for numelem in numelems
        solver, area, q, jac = implicit_solve(p, numelem)
        adj_sln = adjointSolve(solver, area, q, jac, functional)
        # sort nodes
        x_s = zeros(solver.x)  
        sortVec!(solver, solver.x, x_s)
        # sort vector to be printed out
        adj_s = zeros(adj_sln)
        sortVec!(solver, adj_sln, adj_s)

        fout = open("results/adjoints/p$p.dat", "w")
        n_node = size(q, 2) * size(q, 3)
        adj = reshape(adj_s, (3, n_node))
        for i = 1 : n_node
            println(fout, x_s[i], " ", adj[1,i], " ", adj[2,i], " ", adj[3,i])
        end
        close(fout)
    end
end

