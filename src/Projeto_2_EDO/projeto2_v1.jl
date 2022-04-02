using DifferentialEquations, Random, LinearAlgebra, BenchmarkTools

function get_G(N)
    r = rand(N, 3)
    r_ij = reshape( [norm(r[i,:] - r[j,:]) for i=1:N, j=1:N], N, N)
    
    G = -0.5*cis.(-r_ij)./(im.*r_ij)
    G[diagind(r_ij)] .= zero(eltype(G))

    return G
end

function estados_atomicos(du, u, p, t)
    G, Ωₙ, Δ, N, Γ = p

    βₙ = u[1:N]
    zₙ = u[N+1:end]

    for j=1:N
        termo1 = (im*Δ - Γ/2)*βₙ[j]
        termo2 = 0.5*im*Ωₙ[j]*zₙ[j]
        termo3 = sum( G[j,m]*βₙ[m] for m =1:N if j ≠ m  )*zₙ[j]
        du[j] = termo1 + termo2  - termo3
    end
    for j=1:N
        termo1 = im*conj(Ωₙ[j])*βₙ[j] - im*Ωₙ[j]*conj(βₙ[j])
        termo2 = Γ*(1 + zₙ[j])
        termo3 = sum(G[j,m]*βₙ[m]*conj(βₙ[j]) for m =1:N if m≠j)
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

p = G, Ωₙ, Δ, N, Γ

u0 = [zeros(ComplexF64, N); -ones(ComplexF64, N)]
tmin, tmax = 0.0, 10.0
prob = ODEProblem( estados_atomicos, u0, (tmin, tmax), p)

solution = solve(prob, Tsit5(), saveat = 0:1:10);
solução_v1 = solution.u[end];

@benchmark solve(prob, Tsit5(), saveat = 0:1:10)


#=
 Sugestões : 
 - Dá para fazer uma multiplicação matriz-vetor no primeiro loop
 - @views
=#
@profview solve(prob, Tsit5(), saveat = 0:1:10)