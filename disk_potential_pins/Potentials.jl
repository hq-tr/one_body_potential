module Potentials
include("Misc.jl")
include("FQH_state_v2.jl")
using .FQH_states
using .MiscRoutine
using LinearAlgebra
using SparseArrays
using Arpack



function diracdelta_element!(mat::SparseMatrixCSC{ComplexF64, Int64}, i::Int, j::Int, Lam::BitVector,Mu::BitVector, pos::Number)
    R = abs(pos)
    θ = angle(pos) + π
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        term = (-1)^(a+b) * (-1)^(Mu_b-Lam_a) * π * R^(Lam_a+Mu_b) * exp(-R^2/2) * exp(im*(Mu_b-Lam_a)*θ)/((√2)^(Lam_a+Mu_b)*sqfactorial(Lam_a)*sqfactorial(Mu_b))
        mat[i,j] += term
        if i!=j mat[j,i] += conj(term) end
    elseif count_difference == 0
        #println(bin2dex(Lam))
        for m in bin2dex(Lam)
            mat[i,j] += π*R^(2m) * exp(-R^2/2)/(2^m*factorial(big(m)))
        end
    end    
end

function diracdelta_matrix(basis_list::Vector{BitVector}, pos::Number,shift=0.)
    dim = length(basis_list)
    mat = spzeros(Complex{Float64},(dim,dim))
    for i in 1:dim
        print("\r$i\t")
        for j in i:dim
            diracdelta_element!(mat, i, j, basis_list[i], basis_list[j], pos)
        end
    end
    if shift!=0 mat += shift * sparse(I, dim, dim) end
    return mat
end 

function diracdelta_matrix(basis_list::Vector{BitVector}, pos::Vector{T} where T<:Number,shift=0.)
    dim = length(basis_list)
    mat = spzeros(Complex{Float64},(dim,dim))
    for i in 1:dim
        print("\r$i\t")
        for j in i:dim
            for po in pos
                diracdelta_element!(mat, i, j, basis_list[i], basis_list[j], po)
            end
        end
    end
    if shift!=0 mat += shift * sparse(I, dim, dim) end
    return mat
end 

function diracdelta_groundstate(basis_list::Vector{BitVector}, pos::Number, shift=0.;fname="")
    mat = diracdelta_matrix(basis_list, pos, shift)
    E, V = eigs(mat, nev=5, sigma=0)
    println("Energy eigenvalues:")
    println(E)
    gs = FQH_state(basis_list, V[:,1])
    if !isempty(fname) printwf(gs, fname=fname) end
    return gs    
end

function diracdelta_groundstate(basis_list::Vector{BitVector}, pos::Vector{T} where T<:Number, shift=0.;fname="groundstate")
    # More than one trap
    dim = length(basis_list)
    mat = spzeros(Complex{Float64}, (dim,dim))
    for p in pos
        mat += diracdelta_matrix(basis_list, p)
    end
    if shift!=0
        mat += shift * sparse(I,dim,dim)
    end
    E, V = eigs(mat, nev=5, sigma=0)
    println("Energy eigenvalues:")
    println(E)
    if !isempty(fname) printwf(gs; fname=fname) end
    return gs    
end

export diracdelta_matrix, diracdelta_groundstate

end