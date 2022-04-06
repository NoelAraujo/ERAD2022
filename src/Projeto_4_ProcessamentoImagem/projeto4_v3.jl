using Images
"""
    Faço uma divisão explicita por Float64 para o calculo das médias
"""
function extrairBordas_v3!(originalImage, newImage, threshold)
    height, width = size(originalImage)
    for x = 2:height
        for y = 2:width
            currentPixel = originalImage[x,y]
            leftPixel    = originalImage[x-1,y]
            bottomPixel  = originalImage[x,y-1]

            currentAverage= (red(currentPixel)+ green(currentPixel)+ blue(currentPixel))/3.0
            leftAverage   = (red(leftPixel)   + green(leftPixel)   + blue(leftPixel))   /3.0
            bottomAverage = (red(bottomPixel) + green(bottomPixel) + blue(bottomPixel)) /3.0

            if abs(currentAverage-leftAverage) ≥ threshold || abs(currentAverage-bottomAverage) ≥ threshold
                newImage[x,y] = RGBA(0,0,0.,1)
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

@time extrairBordas_v3!(originalImage, newImage, threshold)
[originalImage; newImage]

using BenchmarkTools
@benchmark extrairBordas_v3!(originalImage, newImage, threshold)

## não tem mais o que fazer para o código sequencial
@profview extrairBordas_v3!(originalImage, newImage, threshold)

