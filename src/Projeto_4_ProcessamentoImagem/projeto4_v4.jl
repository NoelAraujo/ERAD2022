using Images
# Algumas vezes o Caminho Relativo não identifica arquivos ou pastas.
# Então precisa pegar o caminho atual com o `pwd()` - mesmo pwd do shell do linux
include(pwd()*"/src/Projeto_4_ProcessamentoImagem/cores.jl")

"""
    Multithreading Nativo
        Cada thread tem uma cor diferente
"""
function extrairBordas_v4!(originalImage, newImage, threshold, cor)
    height, width = size(originalImage)
    @sync for x = 2:height # '@sync' wait to finish
        Threads.@threads for y = 2:width
            currentPixel = originalImage[x,y]
            leftPixel    = originalImage[x-1,y]
            bottomPixel  = originalImage[x,y-1]

            currentAverage= (red(currentPixel)+ green(currentPixel)+ blue(currentPixel))/3.0
            leftAverage   = (red(leftPixel)   + green(leftPixel)   + blue(leftPixel))   /3.0
            bottomAverage = (red(bottomPixel) + green(bottomPixel) + blue(bottomPixel)) /3.0

            if abs(currentAverage-leftAverage) ≥ threshold || abs(currentAverage-bottomAverage) ≥ threshold
                ## Tenho 8 threads mas quero usar apenas 4 cores (para deixar mais limpo o resultado)
                ## então eu uso a função 'mod' para limitar as threads de [0..3]. 
                ## Como não existe thread '0', eu faço '+1', assim fica [1..4]
                newImage[x,y] = cor[ mod(Threads.threadid(), 4) + 1 ]
            else
                newImage[x,y] = RGBA(1,1,1.,1)
            end
        end
    end
    return nothing
end


originalImage = load("src/Projeto_4_ProcessamentoImagem/image2.png")
newImage = deepcopy(originalImage)
threshold = 0.02

Threads.nthreads() |> println
## ---> ver arquivo 'cor.jl' <---
@time extrairBordas_v4!(originalImage, newImage, threshold, cor)
[originalImage; newImage]

using BenchmarkTools
@benchmark extrairBordas_v4!(originalImage, newImage, threshold, cor)


