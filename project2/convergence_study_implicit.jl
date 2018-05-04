# include("assign2_example.jl")
include("solve_implicit.jl")

j1ex = -0.35194635479522557
j2ex = 0.6896586332699256
n_elems = [5, 10, 20, 40, 80]

if !isdir("results")
    mkdir("results")
end

if !isdir("results/functional")
    mkdir("results/functional")
end

for p = 1 : 4
    fout = open("results/functional/p$p.dat", "w")
    for n_elem in n_elems
        solver, area, q = implicit_solve(p, n_elem)

        J1 = calcIntegratedSource(solver, area, q)
        J2 = calcOutletPressure(solver, area, q)
        J1 = abs(J1 - j1ex)
        J2 = abs(J2 - j2ex)
        println(fout, 1.0/n_elem, " ", J1, " ", J2)
        flush(fout)
    end
    close(fout)
end
