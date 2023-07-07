# ------------------ MODULE FOR CALCULATING THE ELECTRON DENSITY
module ParticleDensity

include("Misc.jl")
using .MiscRoutine

single_particle_state_disk(z::Number,m::Integer) = z.^m * exp.(-abs.(z)^2) * sqrt(2^m / (2π)) / sqfactorial(m)

single_particle_state_disk(z::Vector{T} where T <: Number,m::Integer) = z.^m * exp.(-abs.(z)^2) * sqrt(2^m / (2π)) / sqfactorial(m)

single_particle_state_sphere(θ::Number, φ::Number, S::Number, m::Number) = cos(θ/2)^(S+m) * sin(θ/2)^(S-m) * exp(m*im*φ) / sphere_coef(S,m)

function density_element_disk(Lam::BitVector,Mu::BitVector, x::Real,y::Real)
    z = 0.5(x+y*im)
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        return (-1)^(a+b) * conj(single_particle_state_disk(z, Lam_a)) * single_particle_state_disk(z,Mu_b)
    elseif count_difference == 0
        #println("Same partition")
        #println(bin2dex(Lam))
        return sum([abs(single_particle_state_disk(z,m))^2 for m in bin2dex(Lam)])
    else
        return 0
    end    
end

function density_element_sphere(Lam::BitVector,Mu::BitVector, θ::Real, φ::Real)
    No = length(Lam)
    S  = (No-1)/2.
    println(S)
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        return (-1)^(a+b) * conj(single_particle_state_sphere(θ, φ, S, S-Lam_a)) * single_particle_state_sphere(θ,φ,S,S-Mu_b)
    elseif count_difference == 0
        #println("Same partition")
        #println(bin2dex(Lam))
        return sum([abs(single_particle_state_sphere(θ,φ,S,S-m))^2 for m in bin2dex(Lam)])
    else
        return 0
    end    
end

function density_element_gen(Lam::BitVector,Mu::BitVector, single_particle_function::Vector{T} where T <: Number)
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        return (-1)^(a+b) * conj(single_particle_function[Lam_a+1]) * single_particle_function[Mu_b+1]
    elseif count_difference == 0
        #println("Same partition")
        #println(bin2dex(Lam))
        return sum([abs2(single_particle_function[m+1]) for m in bin2dex(Lam)])
    else
        return 0
    end    
end

function density_element_gen(Lam::BitVector,Mu::BitVector, single_particle_function::Vector{Matrix{ComplexF64}})
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        return (-1)^(a+b) .* conj.(single_particle_function[Lam_a+1]) .* single_particle_function[Mu_b+1]
    elseif count_difference == 0
        #println("Same partition")
        #println(bin2dex(Lam))
        return sum([abs2.(single_particle_function[m+1]) for m in bin2dex(Lam)])
    else
        return zeros(size(single_particle_function[1]))
    end    
end

function density_element_gen!(den::Matrix{Float64}, coef::Number, Lam::BitVector,Mu::BitVector, single_particle_function::Vector{Matrix{ComplexF64}})
    check_difference = Lam .⊻ Mu
    count_difference = count(check_difference)
    if count_difference == 2
        Lam_a = bin2dex(check_difference.*Lam)[1]
        Mu_b  = bin2dex(check_difference.*Mu)[1]
        a = count(Lam[1:Lam_a])
        b = count(Mu[1:Mu_b])
        den .+= real.(coef * (-1)^(a+b) * conj.(single_particle_function[Lam_a+1]) .* single_particle_function[Mu_b+1])
    elseif count_difference == 0
        #println(bin2dex(Lam))
        den .+= real.(coef*sum([abs2.(single_particle_function[m+1]) for m in bin2dex(Lam)]))
    end    
end

export density_element_disk, density_element_sphere, density_element_gen, single_particle_state_disk, single_particle_state_sphere, density_element_gen!


end
