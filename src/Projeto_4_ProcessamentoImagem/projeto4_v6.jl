using Images
using ParallelStencil
# @everywhere begin
#     import Pkg; Pkg.add("ParallelStencil")
# end

const USE_GPU = false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2); # 2 == dimensão (trabalhamos com matrizes 2D)
else
    @init_parallel_stencil(Threads, Float64, 2); # como não é todo mundo que tem CUDA, usa as Threads
end

"""
    Criar CUDA Kernel. Mas para quem não tem CUDA disponivel, pode usar ParallelStencil

        https://github.com/omlins/ParallelStencil.jl
"""
@parallel_indices (x,y) function extrairBordas_v6!(originalImage, newImage, threshold)
    if (    x >= 2  && x <= size(originalImage,1) &&
            y >= 2  && y <= size(originalImage,2) )
        
        currentPixel = originalImage[x,y]
        leftPixel    = originalImage[x-1,y]
        bottomPixel  = originalImage[x,y-1]

        if abs(currentPixel-leftPixel) ≥ threshold || abs(currentPixel-bottomPixel) ≥ threshold
            newImage[x-1,y-1] = 0.0
        end
    end
    return
end

originalImage = load("src/Projeto_4_ProcessamentoImagem/image2.png")
nx,ny = size(originalImage)

imageInput = @zeros(nx,ny)
#  'Data.Array' = backend-agnostic initialization function from ParallelStencil
imageInput .= Data.Array(  (red.(originalImage) + green.(originalImage) + blue.(originalImage))/3  )

imageOutput = @ones(nx,ny) # quando eu for colocar em escala de cinza, '1' será o branco

threshold = 0.02
@time @parallel extrairBordas_v6!(imageInput, imageOutput, threshold)
[RGB.(Array(originalImage)); RGB.(Array(imageOutput))]


using BenchmarkTools
@benchmark @parallel extrairBordas_v6!(imageInput, imageOutput, threshold)