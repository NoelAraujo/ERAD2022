"""
    Criar CUDA Kernel. Mas para quem não tem CUDA disponivel, pode usar ParallelStencil

        https://github.com/omlins/ParallelStencil.jl
"""

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

@parallel_indices (x,y) function extrairBordas_v6!(originalImage, newImage, threshold)
    if (    x >= 2  && x <= size(originalImage,1) &&
            y >= 2  && y <= size(originalImage,2) )
        
        currentPixel = originalImage[x,y]
        leftPixel    = originalImage[x-1,y]
        bottomPixel  = originalImage[x,y-1]

        oldAverage    = (red(currentPixel)+ green(currentPixel)+ blue(currentPixel))/3
        leftAverage   = (red(leftPixel)   + green(leftPixel)   + blue(leftPixel))   /3
        bottomAverage = (red(bottomPixel) + green(bottomPixel) + blue(bottomPixel)) /3

        if abs(oldAverage-leftAverage) ≥ threshold || abs(oldAverage-bottomAverage) ≥ threshold
            newImage[x-1,y-1] = RGBA(0,0,0.,1)
        else
            newImage[x-1,y-1] = RGBA(1,1,1.,1)
        end
    end
    return
end

originalImage = load("ERAD2022/src/Projeto_4_ProcessamentoImagem/image2.png")
newImage = deepcopy(originalImage)
threshold = 0.02

nx,ny = size(originalImage)
BLOCKX = 512
BLOCKY = 512

GRIDX = nx÷BLOCKX + 1
GRIDY = ny÷BLOCKY + 1

cuthreads = (BLOCKX, BLOCKX, 1)
cublocks = (GRIDX, GRIDY, 1)


@time @parallel cublocks cuthreads extrairBordas_v6!(originalImage, newImage, threshold)
[originalImage; newImage]

using BenchmarkTools
@benchmark @parallel cublocks cuthreads extrairBordas_v6!(originalImage, newImage, threshold)