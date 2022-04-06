using DifferentialEquations, Random, LinearAlgebra, BenchmarkTools

function get_G(N)
    r = rand(N, 3)
    r_ij = reshape( [norm(r[i,:] - r[j,:]) for i=1:N, j=1:N], N, N)
    
    G = -0.5*cis.(-r_ij)./(im.*r_ij)
    G[diagind(r_ij)] .= zero(eltype(G))

    return G
end
"""
    - posso colocar @views na criação da função
"""
@views function estados_atomicos_v3(du, u, p, t)
    G, Ω, Δ, N, Γ, _temp1 = p

    β = u[1:N]
    z = u[N+1:end]
    
    # Removo a somatória e converto em  multiplicação matricial FORA do loop
    mul!(_temp1, G, β) # _temp1 = G*βₙ
    for j=1:N
        termo1 = (im*Δ - Γ/2)*β[j]
        termo2 = 0.5*im*Ω[j]*z[j]
        termo3 = _temp1[j]*z[j]
        du[j] = termo1 + termo2  - termo3
    end
    
    # Deveria criar outra variavel temporária:
    # _temp2 .= (G*βₙ .- diag(G).*βₙ) 

    # Mas a diagonal já é nula, então não preciso fazer nada,
    # pois '_temp2' já é equivalente a '_temp1'
    for j=1:N
        termo1 = im*conj(Ω[j])*β[j] - im*Ω[j]*conj(β[j])
        termo2 = Γ*(1 + z[j])
        termo3 = _temp1[j]*conj(β[j]) # apenas utilizo a variavel _temp1 já calculada
        du[N+j] = termo1 - termo2  + (2/Γ)*(  termo3 +  conj.(termo3)  )
    end 
    
    return nothing
end

const Γ = 1
N = 100

Random.seed!(2022)
G = get_G(N)
T = eltype(G)

Ω = rand(T, N)
Δ = rand()
_temp1 = similar(Ω)
p = G, Ω, Δ, N, Γ, _temp1

u0 = [zeros(T, N); -ones(T, N)]
tmin, tmax = 0.0, 10.0
prob_v3 = ODEProblem( estados_atomicos_v3, u0, (tmin, tmax), p);

@time solution_v3 = solve(prob_v3, Tsit5(), saveat = 0:1:10);
solução_v1 .≈ solution_v3.u[end]

@benchmark solve(prob_v3, Tsit5(), saveat = 0:1:10)


#=
 Conclusão : Não tem mais espaço para optimização 
=#
@profview solve(prob_v3, Tsit5())