"""
    Apenas cria uma função com o Kernel da computação
"""

using Images

function extrairBordas_v2!(originalImage, newImage, threshold)
    height, width = size(originalImage)
    for x = 2:height
        for y = 2:width
            currentPixel = originalImage[x,y]
            leftPixel    = originalImage[x-1,y]
            bottomPixel  = originalImage[x,y-1]

            oldAverage    = (red(currentPixel)+ green(currentPixel)+ blue(currentPixel))/3
            leftAverage   = (red(leftPixel)   + green(leftPixel)   + blue(leftPixel))   /3
            bottomAverage = (red(bottomPixel) + green(bottomPixel) + blue(bottomPixel)) /3

            if abs(oldAverage-leftAverage) ≥ threshold || abs(oldAverage-bottomAverage) ≥ threshold
                newImage[x,y] = RGBA(0,0,0.,1)
            else
                newImage[x,y] = RGBA(1,1,1.,1)
            end
        end
    end
    return nothing
end

originalImage = load("ERAD2022/src/Projeto_4_ProcessamentoImagem/image2.png")
newImage = deepcopy(originalImage)
threshold = 0.02

@time extrairBordas_v2!(originalImage, newImage, threshold)
[originalImage; newImage]

using BenchmarkTools
@benchmark extrairBordas_v2!(originalImage, newImage, threshold)

# analisar onde posso melhorar com profiler
@profview extrairBordas_v2!(originalImage, newImage, threshold)
# problema: conversão de tipos Float64, para que ">=" não faça conversão