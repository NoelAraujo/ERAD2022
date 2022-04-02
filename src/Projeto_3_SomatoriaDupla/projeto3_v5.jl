using LinearAlgebra, Random, MKL

"""
    (9) using complex number relation: 2real(z) = z + conj(z).  Therefore, I don't need to compute 
        `conj(z)` terms, which are the lower diagonal of βₙₘ and rₙₘ. 
    (9.1) I return: 2real(intensity)
        
    (10) To avoid the `if n≠m`, I begin with upper diagonal  `m = (n+1):N`. Example:
        for n=1:N
            for m = (n+1):N  # <--- here is the change
                ...
            end
        end
    (10.1) The consequence is that I have to specify the exact number of terms in upper matrix
"""
@views function scattering_v5(β, n̂, r)
    N = length(β)
    number_configurations = ((N^2)÷2 - N÷2)
    
    βₙₘ = Array{eltype(β)}(undef, number_configurations)
    count = 1
    for n=1:N
        for m=(n+1):N            
            βₙₘ[count] = conj(β[n])*β[m]
            count += 1
        end
    end

    rₙₘ = Array{eltype(r)}(undef, 3, number_configurations)
    count = 1
    for n=1:N
        r_n = r[:,n]
        for m=(n+1):N
            rₙₘ[1,count] = r_n[1] - r[1,m]
            rₙₘ[2,count] = r_n[2] - r[2,m]
            rₙₘ[3,count] = r_n[3] - r[3,m]
            count += 1
        end
    end

    intensity = zero(eltype(β))
    count = 1    
    for n=1:N
        for m=(n+1):N
            dot_n_r = n̂[1]*rₙₘ[1, count] + n̂[2]*rₙₘ[2, count] + n̂[3]*rₙₘ[3, count]
            intensity += βₙₘ[count]*cis( dot_n_r )
            count += 1
        end
    end
    return 2real(intensity)
end

N = 3000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N) # julia é 'collumn-major'

@time I_v5 = scattering_v5(β, n̂, r);

I_v1 ≈ I_v5
