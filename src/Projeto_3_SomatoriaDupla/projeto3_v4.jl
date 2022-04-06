using LinearAlgebra, Random

"""
   (7) remove the `if n≠m` condition from `βₙₘ` and `rₙₘ`, and only retain in the intensity
   (8) stored in a matrix the vectors of `rₙ - rₘ`
   (9) `count` variable avoids me to think about indices positions inside the matrix
"""
@views function scattering_v4(β, n̂, r)
    N = length(β)

    βₙₘ = Array{eltype(β)}(undef, N,N)
    for n=1:N
        for m=1:N # (7)
            βₙₘ[n,m] = conj(β[n])*β[m]
        end
    end

    rₙₘ = Array{eltype(r)}(undef, 3, N^2) # (8)
    count = 1
    for n=1:N
        r_n = r[:,n]
        for m=1:N # (7)
            # por tentativa e erro, eu vi que foi mais rápido
            # acessar cada posição  do Array
            rₙₘ[1,count] = r_n[1] - r[1,m]
            rₙₘ[2,count] = r_n[2] - r[2,m]
            rₙₘ[3,count] = r_n[3] - r[3,m]
            count += 1 # (9)
        end
    end

    intensity = zero(eltype(β))
    count = 1    
    for n=1:N
        for m=1:N
            if n≠m  # (7)
                dot_n_r = n̂[1]*rₙₘ[1, count] + n̂[2]*rₙₘ[2, count] + n̂[3]*rₙₘ[3, count]  # (9)
                intensity += βₙₘ[n,m]*cis( dot_n_r )
            end
            count += 1
        end
    end
    return intensity
end

N = 3000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N)

@time I_v4 = scattering_v4(β, n̂, r);

I_v1 ≈ I_v4
