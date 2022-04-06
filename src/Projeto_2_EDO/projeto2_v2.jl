using DifferentialEquations, Random, LinearAlgebra, BenchmarkTools

function get_G(N)
    r = rand(N, 3)
    r_ij = reshape( [norm(r[i,:] - r[j,:]) for i=1:N, j=1:N], N, N)
    
    G = -0.5*cis.(-r_ij)./(im.*r_ij)
    G[diagind(r_ij)] .= zero(eltype(G))

    return G
end

"""
    - @views
    - multiplicação matricial com 'mul!' para salvar respostas na variável temporaria
"""
function estados_atomicos_v2(du, u, p, t)
    G, Ω, Δ, N, Γ, _temp1 = p

    # @views fará um acesso do Vetor sem criar uma cópia temporária
    β = @views u[1:N]
    z = @views u[N+1:end]
    
    # Removo a somatória e converto em  multiplicação matricial FORA do loop
    mul!(_temp1, G, β) # _temp1 = G*βₙ
    for j=1:N
        termo1 = (im*Δ - Γ/2)*β[j]
        termo2 = 0.5*im*Ω[j]*z[j]
        termo3 = _temp1[j]*z[j]
        du[j] = termo1 + termo2  - termo3
    end
    

    for j=1:N
        termo1 = im*conj(Ω[j])*β[j] - im*Ω[j]*conj(β[j])
        termo2 = Γ*(1 + z[j])
        termo3 = sum(G[j,m]*β[m]*conj(β[j]) for m =1:N if m≠j)
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
prob_v2 = ODEProblem( estados_atomicos_v2, u0, (tmin, tmax), p);


@time solution_v2 = solve(prob_v2, Tsit5(), saveat = 0:1:10);
solução_v1 .≈ solution_v2.u[end]


@benchmark solve(prob_v2, Tsit5(), saveat = 0:1:10)


#=
 Sugestões : 
 - Melhorar a segunda Somatória
=#
@profview solve(prob_v2, Tsit5());