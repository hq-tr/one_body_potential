## In this version, the basis states in FQH_state variable is stored as BitVector instead of Vector{String}


# -------------- MODULE FOR VECTOR SPACE OPERATIONS ON FQH STATES
module FQH_states

using LinearAlgebra
import Base.+, Base.*
import LinearAlgebra.⋅
import Base.display

include("Misc.jl")
using .MiscRoutine


abstract type  AbstractFQH_state end

struct FQH_state <: AbstractFQH_state
    basis::Vector{BitVector}
    coef::Vector{Number}
end


mutable struct FQH_state_mutable <: AbstractFQH_state
    basis::Vector{BitVector}
    coef::Vector{Number}
end



disk_coef(M,m) = sqfactorial(M) * sqfactorial(m)
# single-particle coefficient on the disk (Landau gauge)

sphere_coef(S,m) =   sqfactorial(S-m)/sqfactorial(S+m+1, 2S+1)

# BASIC HANDLING OF A STATE
dim(vec::AbstractFQH_state) = length(vec.basis)

countorbital(vec::AbstractFQH_state) = length(vec.basis[1])

wfnorm(vec::AbstractFQH_state) = norm(vec.coef)

wfnormalize(vec::FQH_state) = FQH_state(vec.basis, vec.coef/wfnorm(vec))

state2dict(state::AbstractFQH_state) = Dict(state.basis[i] => state.coef[i] for i in 1:dim(state))

function sphere_normalize(vec::FQH_state)
    S = (countorbital(vec)-1.)/2.
    multiplier = [prod([sphere_coef(S,S-m) for m in bin2dex(config)]) for config in vec.basis]
    new_coef = vec.coef.*multiplier
    new_coef /= norm(new_coef)
    return FQH_state(vec.basis, new_coef)
end

function disk_normalize(vec::FQH_state)
    S = (countorbital(vec)-1.)/2.
    multiplier = [prod([sqfactorial(m)/sqrt(2^m) for m in bin2dex(config)]) for config in vec.basis]
    new_coef = vec.coef.*multiplier
    new_coef /= norm(new_coef)
    return FQH_state(vec.basis, new_coef)
end

function wfnormalize!(vec::FQH_state_mutable)
    vec.coef /= wfnorm(vec)
end

function sphere_normalize!(vec::FQH_state_mutable)
    S = (countorbital(vec)-1.)/2.
    multiplier = [prod([sphere_coef(S,S-m) for m in bin2dex(config)]) for config in vec.basis]
    vec.coef .*= multiplier
    vec.coef /= norm(vec.coef)
    return
end

function disk_normalize!(vec::FQH_state_mutable)
    multiplier = [prod([sqfactorial(m)/sqrt(2^m) for m in bin2dex(config)]) for config in vec.basis]
    vec.coef .*= multiplier
    vec.coef /= norm(new_coef)
    return
end




# ------------I/O
function printwf(state::AbstractFQH_state; fname = "")
    D = dim(state)
    if length(fname) == 0
        println(D)
        for i in 1:D
            println(prod(string.((Int.(state.basis[i])))))
            println(replace("$(state.coef[i])", "im"=>"j", " "=>""))
        end
    else
        open(fname,"w") do f
            write(f,"$D\n")
            for i in 1:D
                writebasis = prod(string.((Int.(state.basis[i]))))
                writecoef  = replace("$(state.coef[i])", "im"=>"j", " "=>"")
                write(f,"$writebasis\n$writecoef\n")
            end 
        end

    end
end

function display(vec::AbstractFQH_state) printwf(vec) end

function readwf(fname::String; mutable=false)
        f = open(fname)
        content = readlines(f)
        dim = parse(Int64, content[1])
        basis = [BitVector(map(x->parse(Bool, x), split(y,""))) for y in content[2:2:end]]
        #println(basis[1])
        println("The dimension is $dim")
        try
            global co = [parse(Float64, x) for x in content[3:2:end]]
        catch ArgumentError
            println("Reading coefficients as complex numbers")
            global co = [parse(Complex{Float64}, x) for x in content[3:2:end]]
        finally
            close(f)
            #println("Success")
        end
        if mutable
            return FQH_state_mutable(basis,co)
        else
            return FQH_state(basis, co)
        end
end


#------------- COLLATE VECTORS WITH DIFFERENT BASIS
# (Useful for subsequent operations)
function collate_vector(vec1::FQH_state, vec2::FQH_state)
    new_basis = collect(Set(vec1.basis) ∪ Set(vec2.basis))
    dict1 = state2dict(vec1)
    dict2 = state2dict(vec2)
    new_coef1 = map(x -> (try dict1[x] catch KeyError 0 end), new_basis)
    new_coef2 = map(x -> (try dict2[x] catch KeyError 0 end), new_basis)
    return FQH_state(new_basis, new_coef1), FQH_state(new_basis, new_coef2)
end


function collate_many_vectors(vectors::Vector{FQH_state};separate_out=false, collumn_vector=false)
    if separate_out
        N = length(vectors)
        new_basis = collect(Set(Iterators.flatten(x.basis for x in vectors)))
        dim = length(new_basis)
        println("Total dimension is $(dim)")
        if collumn_vector
            new_coef = zeros((dim,N))
        else
            new_coef = zeros((N,dim))
        end

        for i in 1:N
            vec = state2dict(vectors[i])
            if collumn_vector
                new_coef[:,i] = map(x -> (try vec[x] catch KeyError 0 end), new_basis)
            else
                new_coef[i,:] = map(x -> (try vec[x] catch KeyError 0 end), new_basis)
            end
        end
        newvectors = (new_basis, new_coef)
    else
        new_basis = collect(Set(Iterators.flatten(x.basis for x in vectors)))
        newvectors = map(vec -> FQH_state(new_basis,map(x -> (try state2dict(vec)[x] catch KeyError 0 end), new_basis)), vectors)
    end
    return newvectors
end

# ------------ VECTOR SPACE OPERATIONS

function Base.:+(vec1::FQH_state, vec2::FQH_state) # Addition
    v1,v2 = collate_vector(vec1, vec2)
    return FQH_state(v1.basis, v1.coef+v2.coef)
end

function Base.:*(vec1::FQH_state, multiplier::Number) # Scalar multiplication
    return FQH_state(vec1.basis, multiplier * vec1.coef)
end

Base.:*(multiplier::Number, vec1::FQH_state) = vec1*multiplier

function overlap(vec1::FQH_state, vec2::FQH_state) # Inner product
    v1, v2 = collate_vector(vec1, vec2)
    return v1.coef⋅v2.coef
end

LinearAlgebra.:⋅(vec1::FQH_state, vec2::FQH_state) = overlap(vec1,vec2)


# -------- Calculate density 
include("Density.jl")
using .ParticleDensity

function get_density_disk(vec::FQH_state, x::Real, y::Real)
    z = 0.5(x + y*im)
    single_particle = [single_particle_state_disk(z,m) for m in 0:countorbital(vec)]
    D = dim(vec)
    element(i,j) = 2^(i!=j) * conj(vec.coef[i]) * vec.coef[j] * density_element_gen(vec.basis[i], vec.basis[j], single_particle)
    density = abs(sum(sum(element(i,j) for i in j:D) for j in 1:D))
    #density = sum(abs(vec.coef[i])^2 * density_element_gen(vec.basis[i], vec.basis[i], single_particle) for i in 1:D)
    #density += 2*sum(real(vec.coef[j]*sum(conj(vec.coef[i]) * density_element_gen(vec.basis[i], vec.basis[j], single_particle) for i in (j+1):D)) for j in 1:(D-1))
    return density
end

function get_density_disk(vec::FQH_state, x::Array{T} where T<: Number, y::Array{T} where T <: Number)
    z = 0.5*(x .+ y*im)
    single_particle = [single_particle_state_disk.(z,m) for m in 0:countorbital(vec)]
    D = dim(vec)

    den = zeros(size(z))
    for i in 1:D
        for j in i:D 
            coef = 2^(i!=j) * conj(vec.coef[i]) * vec.coef[j]
            density_element_gen!(den, coef, vec.basis[i], vec.basis[j], single_particle)
        end
    end
    return den
end

export AbstractFQH_state, FQH_state, FQH_state_mutable, readwf, printwf, sphere_normalize, disk_normalize, wfnormalize, dim, get_density_disk, overlap, +, *, ⋅, collate_many_vectors

end # ----- END MODULE