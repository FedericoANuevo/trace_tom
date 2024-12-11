PRO ps1,archivo
set_plot,'ps'
device,filename=archivo,bits_per_pixel=8,/color,/encapsulated
return
end
