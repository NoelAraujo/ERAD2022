using DifferentialEquations, Random, LinearAlgebra, BenchmarkTools

import Pkg; Pkg.add("CUDA")
using CUDA

function get_G(N)
    r = rand(N, 3)
    r_ij = reshape( [norm(r[i,:] - r[j,:]) for i=1:N, j=1:N], N, N)
    
    G = -0.5*cis.(-r_ij)./(im.*r_ij)
    G[diagind(r_ij)] .= zero(eltype(G))

    return G
end

function estados_atomicos_v4(du, u, p, t)
    G, Ωₙ, G_βₙ, Wₙ, Δ, N, Γ, _temp1, _temp2 = p

    βₙ = u[1:N]
    zₙ = u[N+1:end]
    
    mul!(G_βₙ, G, βₙ) # == G_βₙ = G*βₙ
    Wₙ     .=  Ωₙ/2 + im.*G_βₙ
    _temp1 .=  (im*Δ - Γ/2) * βₙ.+ im.*Wₙ.*zₙ
    _temp2 .=  -Γ*(1 .+ zₙ) .- 4 * imag.(βₙ.*conj.(Wₙ))
    du[:]  .=  vcat(_temp1, _temp2)
    
    return nothing
end

const Γ = 1
N = 100

Random.seed!(2022)

G      = get_G(N) |> CuArray
Ωₙ     = rand(ComplexF64, N) |> CuArray
Δ      = rand()
Wₙ     = similar(Ωₙ)
G_βₙ   = similar(Ωₙ)
_temp1 = similar(Ωₙ)
_temp2 = similar(Ωₙ)
p      = G, Ωₙ, G_βₙ, Wₙ, Δ, N, Γ, _temp1, _temp2

u0         = [zeros(ComplexF64, N); -ones(ComplexF64, N)] |> CuArray
tmin, tmax = 0.0, 10.0
prob_v4    = ODEProblem( estados_atomicos_v4, u0, (tmin, tmax), p);

@time solution_v4 = solve(prob_v4, Tsit5(), saveat = 0:1:10);
solução_v1 .≈ Array(solution_v4.u[end])

@benchmark solve(prob_v4, Tsit5(), saveat = 0:1:10)