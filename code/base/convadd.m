function outadd=convadd(inadd)

outadd=dec2hex(hex2dec('FFFFFFFFFFFF')-hex2dec(fliplr(inadd)));
while length(outadd)<12, outadd=['0' outadd]; end