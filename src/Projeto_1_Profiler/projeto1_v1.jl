# -------- GS_process.py --------
using LinearAlgebra
using QuadGK

poly(order) = x->x^order

function l2_inner_product(func1, func2, params...)
    density, a, b, norm_factor = params
    
    proj, err = quadgk(a,b; rtol=1e-12, atol=1e-12) do x
        func1(x)*func2(x)*density(x)/norm_factor
    end
    
    return proj
end

multiply(coefficient, func) = x->coefficient*func(x)
difference(func1, func2) = x->(func1(x) - func2(x))

function orthonormal(order, parameters)    
    if(order == 0)
        norm = sqrt(l2_inner_product(poly(0), poly(0), parameters...))
        
        return x->poly(0)(x)/norm
    
    else
        temp_func = x->poly(order)(x)
        
        for deg in 1:order
            ref_func = orthonormal(deg - 1, parameters)
            proj_deg = l2_inner_product(poly(order), ref_func, parameters...)
            temp_func = difference(temp_func, multiply(proj_deg, ref_func))
        end
        norm = sqrt(l2_inner_product(temp_func, temp_func, parameters...))
        normed_func = x->temp_func(x)/norm
        
        return normed_func
    end
end

# -------- run.py --------
using Plots

density = x->1
lower_bound = -1
upper_bound = 1
density_normalization = 1
parameters = density, lower_bound, upper_bound, density_normalization

plot(lower_bound:0.01:upper_bound, orthonormal(13, parameters))

for deg in range(1,14)
    @time orthonormal(deg, parameters)    #orthonormal function
end

# -------- profiler --------
#=
    Muitas chamadas recursivas de funções
=#
@profview orthonormal(13, parameters)

