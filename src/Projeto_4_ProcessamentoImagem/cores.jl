
function convert_hex_rhb(hex_color)
    r = parse(Int, hex_color[1:2], base = 16) / 255
    g = parse(Int, hex_color[3:4], base = 16) / 255
    b = parse(Int, hex_color[5:6], base = 16) / 255
    return RGBA(r, g, b)
end

# hexColors = ["3B707D", "FFB85D", "FF5412", "030A04"]
hexColors = ["90ADC6", "FF4500", "FAD02C", "333652"]

cor = convert_hex_rhb.(hexColors)