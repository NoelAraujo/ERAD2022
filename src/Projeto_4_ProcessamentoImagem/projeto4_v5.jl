"""
    Usar Processos e Não Threads
        Tenho que dizer qual a região ( Δx e Δy ) cada processo deve trabalhar
        Por algum motivo estranho, tenho que inverter os indices x e y
            pois 'height, width = size(originalImage)' não esta retornando o que esperava
"""

using Distributed
addprocs(4, exeflags = `--threads 2`)
# @everywhere begin 
#     import Pkg; Pkg.add("Images")
#     using Images
# end
@everywhere using SharedArrays
@everywhere using Images

@everywhere function extrairBordas_v5!(originalImage, newImage, threshold, cor, y_min_max, x_min_max)
    @sync for x = x_min_max
        Threads.@threads for y = y_min_max
            currentPixel = originalImage[x,y]
            leftPixel    = originalImage[x-1,y]
            bottomPixel  = originalImage[x,y-1]

            currentAverage= (red(currentPixel)+ green(currentPixel)+ blue(currentPixel))/3.0
            leftAverage   = (red(leftPixel)   + green(leftPixel)   + blue(leftPixel))   /3.0
            bottomAverage = (red(bottomPixel) + green(bottomPixel) + blue(bottomPixel)) /3.0

            if abs(currentAverage-leftAverage) ≥ threshold || abs(currentAverage-bottomAverage) ≥ threshold
                newImage[x,y] = cor[ myid()-1 ]
            else
                newImage[x,y] = RGBA(1,1,1.,1)
            end
        end
    end
    return nothing
end


originalImage = load("ERAD2022/src/Projeto_4_ProcessamentoImagem/image2.png")
# Preciso de uma variável compartilhada entre todos os processos locais: SharedArrays
_height, _width = size(originalImage)
originalImage_Shared  = SharedArray{eltype(originalImage)}(_height, _width);
originalImage_Shared .= originalImage;
newImage_Shared = deepcopy(originalImage_Shared);

threshold = 0.02

# preciso definir os pedaços da imagem 
x_min_max, y_min_max = 2:_width, 2:_height
chunkSize = length(x_min_max) ÷ (nworkers())
x_chunks = Base.Iterators.partition(x_min_max, chunkSize)
collect(x_chunks)

pmap(x_chunks) do pedaço_x
    extrairBordas_v5!(originalImage_Shared, newImage_Shared, threshold, cor, pedaço_x, y_min_max)
end
[originalImage_Shared; newImage_Shared]


using BenchmarkTools
@benchmark pmap(x_chunks) do x
    extrairBordas_v5!(originalImage_Shared, newImage_Shared, threshold, cor, x, y_min_max)
end


