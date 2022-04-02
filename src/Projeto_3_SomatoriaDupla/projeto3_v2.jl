using LinearAlgebra, MKL, Random

"""
   (1) Put @views to access all collumns without producing SubArray
   (2) Compute `conj(β[n])*β[m]` once
   (3) Avoid accessing position rₙ inside loop m (r_n = r[:, n])
"""
@views function scattering_v2(β, n̂, r)
    N = length(β)
    intensity = ComplexF64(0)
    
    βₙₘ = zeros(ComplexF64, N, N)
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
                rₘ = r[:, m]
                intensity += βₙₘ[n,m]*exp( im*(n̂⋅(rₙ-rₘ)) )
            end
        end
    end
    return intensity
end

N = 3000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N) # julia é 'collumn-major'

@time I_v2 = scattering_v2(β, n̂, r);

I_v1 ≈ I_v2