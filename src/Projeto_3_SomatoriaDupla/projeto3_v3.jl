using LinearAlgebra, Random

"""
   (4) `zeros(N,N)` --> `Array{eltype(β)}(undef, N,N)`
   (5) expand the dot product
   (6) change `exp(im*x)` by `cis(x)`
"""
@views function scattering_v3(β, n̂, r)
    N = length(β)
    intensity = ComplexF64(0)
    
    βₙₘ = Array{eltype(β)}(undef, N,N) # (4)
    for n=1:N
        for m=1:N
            if n≠m
                βₙₘ[n,m] = conj(β[n])*β[m]
            end
        end
    end

    for n=1:N
        rₙ = r[:, n]
        for m=1:N
            if n≠m
                rₙₘ =  rₙ - r[:, m]
                dot_n_r = n̂[1]*rₙₘ[1] + n̂[2]*rₙₘ[2] + n̂[3]*rₙₘ[3] # (5)

                intensity += βₙₘ[n,m]*cis( dot_n_r ) # (6)
            end
        end
    end
    return intensity
end

N = 3000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N)

@time I_v3 = scattering_v3(β, n̂, r);

I_v1 ≈ I_v3
