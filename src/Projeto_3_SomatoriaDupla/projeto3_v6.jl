using LinearAlgebra, MKL, Random, ThreadsX


"""
    (12) use multi-threading of `ThreadsX.mapreduce` to transform double loop into a sum-like expression
    (13) add `@inbounds` to tell the compiler to disable out of memory checking
"""
@views function scattering_v6(β, n̂, r)
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
        @inbounds r_n = r[:,n]
        for m=(n+1):N
            rₙₘ[1,count] = r_n[1] - r[1,m]
            rₙₘ[2,count] = r_n[2] - r[2,m]
            rₙₘ[3,count] = r_n[3] - r[3,m]
            count += 1
        end
    end
     
    intensity = ThreadsX.mapreduce(+, 1:number_configurations) do count # (12)
                    @inbounds dot_n_r = n̂[1]*rₙₘ[1, count] + n̂[2]*rₙₘ[2, count] + n̂[3]*rₙₘ[3, count] # (13)
                    @inbounds βₙₘ[count]*cis( dot_n_r ) # (13)
                end
    return 2real(intensity)
end

N = 3000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N) # julia é 'collumn-major'

@time I_v6 = scattering_v6(β, n̂, r);

I_v1 ≈ I_v6
