include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/Potentials.jl")

using .FQH_states
using .Potentials
using LinearAlgebra
using SparseArrays

#Ne = 10
#No = 2*Ne-1

function main()
    println("------ Potential pins on the sphere ------")
    println()
    println("Input root file name: ")
    fname = readline()
    roots = readwf(fname)
    jack_list = Vector{FQH_state}()

    println("How many potential pins?")
    npins = parse(Int, readline())

    println("Input location θ ϕ of each potential pins (in units of π)")
    θ_list = Float64[]
    ϕ_list = Float64[]
    for i in 1:npins
        println("    Pin #$i:")
        locread = readline()
        θ,ϕ = map(x->parse(Float64, x), split(locread))
        push!(θ_list, θ*π)
        push!(ϕ_list, ϕ*π)
    end


    for root in roots.basis
        rootstring = prod(string.((Int.(root))))
        print("\r$(rootstring)\t")
        jack = sphere_normalize(readwf("jacks/J_$(rootstring)"))
        push!(jack_list, jack)
    end

    Ne = count(jack_list[1].basis[1])
    No = length(jack_list[1].basis[1])
    println("\n\n$(Ne) electrons and $(No) orbitals\n")

    all_basis, all_coef = collate_many_vectors(jack_list; separate_out=true, collumn_vector=true)

    println("Orthonormalizing basis using QR decomposition")
    @time all_coef_ortho = transpose(Matrix(qr(all_coef).Q))

    dim = length(all_basis)
    mat = spzeros(Complex{Float64},(dim,dim))

    println("------")
    
    #R = 3.5

    #pos1 = R
    #pos2 = conj(pos1)
    #pos3 = -pos1
    #pos4 = -pos2

    #mat = diracdelta_matrix(all_basis, [pos1,pos2,pos3,pos4])
    @time for i in 1:npins
        mat += sphere_point_matrix(all_basis, θ_list[i], ϕ_list[i])
    end

    mat += 5 * sparse(I, dim, dim)
    #@time mat = sphere_point_matrix(all_basis,θ_list, ϕ_list,1.0,5.0)

    @time mat_MR = conj.(all_coef_ortho) * mat * transpose(all_coef_ortho)


    ED = eigen(mat_MR)


    println("Eigenvalues = ")
    if length(ED.values) > 3
        display(ED.values[1:5])
        println("\t.\n\t.\n\t.")
    else
        display(ED.values)
    end


    vecs = ED.vectors[:,sortperm(abs.(ED.values))] # Sort by increasing abs(E)


    gs_coef = transpose(all_coef_ortho) * vecs[:,1]

    gs = FQH_state(all_basis, gs_coef)
    printwf(gs;fname="$(fname)_gs_$(npins)pins_0")

    gs2_coef = transpose(all_coef_ortho) * vecs[:,2]
    gs2 = FQH_state(all_basis, gs2_coef)
    printwf(gs2;fname="$(fname)_gs_$(npins)pins_1")


    open("$(fname)_gs_$(npins)pins_eigen", "w") do f
        for val in ED.values
            write(f,"$(real(val))\n")
        end
    end
end

@time main()