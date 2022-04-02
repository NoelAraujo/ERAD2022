using LinearAlgebra, MKL, Random
"""
    Intesidade = ∑ₙ∑ₘ conj(βₙ)βₘ exp( im*n̂⋅(rₙ - rₘ) )

- ∑ₙ∑ₘ :  n,m ∈ [1..N], with n≠m
- N = number of particles
- βⱼ = atom expected value of σ⁻ (j ∈ [1..N])
- n̂ = versor of the sensor ( n = r_sensor./norm(r_sensor) )
- rⱼ = atom position (j ∈ [1..N])
"""
function scattering_v1(β, n̂, r)
    N = length(β)
    intensity = ComplexF64(0)

    for n=1:N
        for m=1:N
            if n≠m
                rₙ, rₘ = r[:, n], r[:, m]
                βₙ, βₘ = β[n], β[m]
                intensity += conj(βₙ)*βₘ*exp( im*( n̂⋅(rₙ-rₘ) ))
                # warning
                # im*(c⋅(a-b))
                # (im*c)⋅(a-b)
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

@time I_v1 = scattering_v1(β, n̂, r);
