# Apenas as comparações

N = 5000

Random.seed!(2022)
β = rand(ComplexF64, N)
sensor = rand(3); n̂ = sensor./norm(sensor)
r = rand(3, N) # julia é 'collumn-major'

print("v1"); @time I_v1 = scattering_v1(β, n̂, r);
print("v2"); @time I_v2 = scattering_v2(β, n̂, r); @assert I_v1 ≈ I_v2
print("v3"); @time I_v3 = scattering_v3(β, n̂, r); @assert I_v1 ≈ I_v3
print("v4"); @time I_v4 = scattering_v4(β, n̂, r); @assert I_v1 ≈ I_v4
print("v5"); @time I_v5 = scattering_v5(β, n̂, r); @assert I_v1 ≈ I_v5
print("v6"); @time I_v6 = scattering_v6(β, n̂, r); @assert I_v1 ≈ I_v6

# @profview scattering_v6(β, n̂, r)