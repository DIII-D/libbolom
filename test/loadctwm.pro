pro loadctwm
s = 36
r = bytscl(hanning(s*2+2))
rampup = r[0:s]
rampdn = r[s+1:s*2+1]
r = bytscl(bytscl(hanning(s*2+2),max=0.8))
frampup = r[0:s]
frampdn = r[s+1:s*2+1]
r = bytarr(256)
g = bytarr(256)
b = bytarr(256)
st = 20
; violet
i = 0
r[i*s+st:(i+1)*s+st] = frampup * 0.4
g[i*s+st:(i+1)*s+st] = 0
b[i*s+st:(i+1)*s+st] = frampup
; blue
i = 1
r[i*s+st:(i+1)*s+st] = frampdn * 0.4
g[i*s+st:(i+1)*s+st] = 0
b[i*s+st:(i+1)*s+st] = 255
; cyan
i = 2
r[i*s+st:(i+1)*s+st] = 0
g[i*s+st:(i+1)*s+st] = rampup
b[i*s+st:(i+1)*s+st] = 255
; green
i = 3
r[i*s+st:(i+1)*s+st] = 0
g[i*s+st:(i+1)*s+st] = 255
b[i*s+st:(i+1)*s+st] = frampdn
i = 4
r[i*s+st:(i+1)*s+st] = rampup
g[i*s+st:(i+1)*s+st] = 255
b[i*s+st:(i+1)*s+st] = 0
; yellow
i = 5
r[i*s+st:(i+1)*s+st] = 255
g[i*s+st:(i+1)*s+st] = rampdn
b[i*s+st:(i+1)*s+st] = 0
; red
i = 6
r[i*s+st:*] = 255
r[255] = 255
g[255] = 255
b[255] = 255
tvlct,r,g,b
return
end
