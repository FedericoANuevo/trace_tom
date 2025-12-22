pro record_gif,dir,filename

image24 = TVRD(True=1)
image2d = Color_Quan(image24, 1, r, g, b)
write_GIF, dir+filename, image2d, r, g, b

return
end
