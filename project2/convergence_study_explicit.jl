include("assign2_example.jl")

j1ex = -0.35194635479522557
j2ex = 0.6896586332699256
n_elems = [5, 10, 20, 40, 80]
for p = 1 : 4
    fout = open("j1j2_p$p.dat", "w")
    for n_elem in n_elems
        J1, J2 = explicit_solve(p, n_elem)
        J1 = abs(J1 - j1ex)
        J2 = abs(J2 - j2ex)
        println(fout, 1.0/n_elem, " ", J1, " ", J2)
        flush(fout)
    end
    close(fout)
end

# for p = 1 : 4
#     n_elem = 5
#     J1, J2 = explicit_solve(p, n_elem)
#     J1 = abs(J1 - j1ex)
#     J2 = abs(J2 - j2ex)
#     println(1.0/n_elem, " ", J1, " ", J2)
#     println(1.0/n_elem, " ", J1, " ", J2)
#     println(1.0/n_elem, " ", J1, " ", J2)
#     println(1.0/n_elem, " ", J1, " ", J2)
# end
