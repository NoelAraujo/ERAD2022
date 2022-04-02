# -------- GS_process.py --------
using LinearAlgebra
using QuadGK
using Memoize
using Symbolics

@variables x

@memoize poly(order) =  x^order

@memoize function l2_inner_product(func1, func2, params)
    density, a, b, norm_factor = params
    
    dotProduct = build_function( func1*func2*density/norm_factor, [x], expression = Val{false} )
    proj, err = quadgk(dotProduct, a,b; rtol=1e-12, atol=1e-12)
    
    return proj
end

@memoize multiply(coefficient, func) = coefficient*func
@memoize difference(func1, func2) = func1 - func2

@memoize function orthonormal(order, parameters)    
    if(order == 0)
        norm = sqrt(l2_inner_product(poly(0), poly(0), parameters))
    
        return poly(0)/norm
    else
        temp_func = poly(order)
        
        for deg in 1:order
            ref_func = orthonormal(deg - 1, parameters)
            proj_deg = l2_inner_product(poly(order), ref_func, parameters)
            temp_func = difference(temp_func, multiply(proj_deg, ref_func))
        end
        norm = sqrt(l2_inner_product(temp_func, temp_func, parameters))
        normed_func = temp_func/norm
        
        return normed_func
    end
end

# -------- run.py --------
using Plots

density = 1 # antes era função, agora é constante
lower_bound = -1
upper_bound = 1
density_normalization = 1
parameters = density, lower_bound, upper_bound, density_normalization

deg = 10
orth_1 = build_function( orthonormal(deg, parameters), [x], expression = Val{false} )

plot(lower_bound:0.01:upper_bound, orth_1)

empty!(memoize_cache(poly))
empty!(memoize_cache(l2_inner_product))
empty!(memoize_cache(multiply))
empty!(memoize_cache(difference))
empty!(memoize_cache(orthonormal))

for deg in range(1,20)
    @time orthonormal(deg, parameters)    #orthonormal function
end

# -------- profiler --------
#=
    Agora a construção da Função a ser integrada demora mais que o calculo numérico 
=#
empty!(memoize_cache(poly))
empty!(memoize_cache(l2_inner_product))
empty!(memoize_cache(multiply))
empty!(memoize_cache(difference))
empty!(memoize_cache(orthonormal))

@profview orthonormal(15, parameters)