module MiscRoutine

# Basis format conversions
bin2dex(config::BitVector) = [i-1 for i in 1:length(config) if config[i]]

dex2bin(state::Vector{Int}, N_orb::Int) = BitVector([i-1 in state for i in 1:N_orb]) 


# Miscellaneous functions for LLL physics

sqfactorial(N) = prod(map(sqrt, 1:N)) 
# square root of N!, avoid overflow for large N and more efficient than sqrt(factorial(big(N)))
# (overflow starts at N=21)

sqfactorial(n,N) = prod(map(sqrt, n:N))

export bin2dex, sqfactorial, dex2bin
end