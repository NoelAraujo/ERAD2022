using DifferentialEquations, Random, LinearAlgebra, BenchmarkTools
using TensorOperations # uso tensores

function get_G(N)
    r = rand(N, 3)
    r_ij = reshape( [norm(r[i,:] - r[j,:]) for i=1:N, j=1:N], N, N)
    
    G = -0.5*cis.(-r_ij)./(im.*r_ij)
    G[diagind(r_ij)] .= zero(eltype(G))

    return G
end

function estados_atomicos_v3(du, u, p, t)
    G, Ωₙ, Δ, N, Γ, _temp1 = p

    βₙ = @views u[1:N]
    zₙ = @views u[N+1:end]
    
    # Removo a somatória e converto em  multiplicação matricial FORA do loop
    mul!(_temp1, G, βₙ) # _temp1 = G*βₙ
    for j=1:N
        termo1 = (im*Δ - Γ/2)*βₙ[j]
        termo2 = 0.5*im*Ωₙ[j]*zₙ[j]
        termo3 = _temp1[j]*zₙ[j]
        du[j] = termo1 + termo2  - termo3
    end
    
    # Deveria criar outra variavel temporária:
    # _temp2 .= (G*βₙ .- diag(G).*βₙ) 

    # Mas a diagonal já é nula, então não preciso fazer nada,
    # pois _temp2 seria equivalente a _temp1
    for j=1:N
        termo1 = im*conj(Ωₙ[j])*βₙ[j] - im*Ωₙ[j]*conj(βₙ[j])
        termo2 = Γ*(1 + zₙ[j])
        termo3 = _temp1[j]*conj(βₙ[j])
        du[N+j] = termo1 - termo2  + (2/Γ)*(  termo3 +  conj.(termo3)  )
    end 
    
    return nothing
end

const Γ = 1
N = 100

Random.seed!(2022)


G = get_G(N)
Ωₙ = rand(ComplexF64, N)
Δ = rand()
_temp1 = similar(Ωₙ)
p = G, Ωₙ, Δ, N, Γ, _temp1

u0 = [zeros(ComplexF64, N); -ones(ComplexF64, N)]
tmin, tmax = 0.0, 10.0
prob_v3 = ODEProblem( estados_atomicos_v3, u0, (tmin, tmax), p);

@time solution_v3 = solve(prob_v3, Tsit5(), saveat = 0:1:10);
solução_v1 .≈ solution_v3.u[end]

@benchmark solve(prob_v3, Tsit5(), saveat = 0:1:10)


#=
 Conclusão : Não tem mais espaço para optimização 
=#
@profview solve(prob_v3, Tsit5())