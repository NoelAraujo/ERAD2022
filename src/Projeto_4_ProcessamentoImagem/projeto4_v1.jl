using Images

originalImage = load("ERAD2022/src/Projeto_4_ProcessamentoImagem/image2.png")
newImage = deepcopy(originalImage)
height, width = size(originalImage)

threshold = 0.02

@time for x = 2:height
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

[originalImage; newImage]