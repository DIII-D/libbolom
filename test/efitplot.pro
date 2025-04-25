pro efitplot,shot,time,sol,matrmin=matrmin,matrmax=matrmax,matzmin=matzmin,matzmax=matzmax,tl=tl,camera=camera,runid=runid,range=range,snowflake=snowflake,pixels=pixels,wavelength=wavelength,filename=filename,xrange=xrange,yrange=yrange,limiter=limiter,nographics=nographics,quiet=quiet
common xptcom,nextwindow
if n_elements(matrmin) eq 0 then matrmin = 100.0
if n_elements(matrmax) eq 0 then matrmax = 200.0
if n_elements(matzmin) eq 0 then matzmin = -140.0
if n_elements(matzmax) eq 0 then matzmax = -40.0
if isa(nographics) eq 0 then nographics = 0
if isa(xrange) eq 0 then xrange = [matrmin,matrmax]
if isa(yrange) eq 0 then yrange = [matzmin,matzmax]
if n_elements(snowflake) eq 0 then snowflake = 0
if n_elements(camera) eq 0 then camera = 'Tangtv'
if n_elements(tl) eq 0 then begin
    tl = 'DIII-D Emissivity '+strtrim(camera,2)+' Shot '+strtrim(string(shot),2)+'@'+strtrim(string(time),2)
    if n_elements(wavelength) gt 0 then tl = tl + "     Filter:"+wavelength
endif
if n_elements(range) eq 0 then crange = [0.0, max(sol)]
if n_elements(range) gt 0 then crange = range
csol = bytscl(sol,min=crange[0],max=crange[1])
lrange = [0,256]
if n_elements(pixels) eq 0 then pixels = 0
s = size(sol)
imxlen = s[1]
imylen = s[2]
img_r = findgen(imxlen) * (matrmax-matrmin)/(imxlen-1) + matrmin
img_z = findgen(imylen) * (matzmax - matzmin) / (imylen-1) + matzmin
img_r = img_r / 100.0
img_z = img_z / 100.0
delr = (img_r[1]-img_r[0]) * 0.5
delz = (img_z[1]-img_z[0]) * 0.5
getefit,shot,time,aaa,ggg,deltime=20,err=err,mode='mdsplus',runid=runid
if err ne 0 then begin
    print,'Problem reading efit for ',shot,'@',time,'status = ',err
    print,'Reconnecting and trying again'
    mdsdisconnect
    mdsconnect,"atlas.gat.com"
    getefit,shot,time,aaa,ggg,deltime=20,err=err,mode='mdsplus',runid=runid
    if err ne 0 then return
endif
_imxlen = float(imxlen)
_imylen = float(imylen)
if isa(filename) ne 0 then begin
    _scale = 8.0 / max([_imxlen,_imylen])
    _imxlen = _imxlen * _scale
    _imylen = _imylen * _scale
    set_plot,'ps',/interpolate,/copy
    device,/inches,/color,xsize=_imxlen,ysize=_imylen,/Helvetica,/TT_FONT,bits_per_pixel=8,decomposed=0,filename=filename,/landscape
endif else begin
    _scale = 1024.0 / max([_imxlen,_imylen])
    _imxlen = _imxlen * _scale
    _imylen = _imylen * _scale
    nextwindow = nextwindow+1
    winifopen,nextwindow,wins=wins,xpos=xpos,ypos=ypos,xsize=_imxlen,ysize=_imylen
endelse
lim = ggg.lim
if isa(limiter) gt 0 then lim=limiter
x = size(lim)
lim_points = x[2]
psin = (ggg.psirz - ggg.ssimag) / (ggg.ssibry - ggg.ssimag)
rmajor = gdp_calc_rmajor(psin,aaa,ggg)
; handle color visual (Pseudo 8bit = default, TrueColor 24bit =1)
colorvisual=0 ; default for 8bit
ncolors=!d.n_colors-2
device, decomposed=0
colorvisual=1
ncolors=255
;loadct,39,/silent ; colour table rainbow+white
loadctwm
savebg = !p.background
savecol = !p.color
!p.background=!d.n_colors-1 ; background = white
!p.color=1
xinc = 1.7 / (ggg.mw - 1)
yinc = 3.2 / (ggg.mh - 1)
psin_r = 0.84 + xinc * findgen(ggg.mw)
psin_z = -1.6 + yinc * findgen(ggg.mh)
; adjust c_colours to your terminal
lblue=(ncolors-1)*0.35 ; light blue for core surfaces
blue=(ncolors-1)*0.3 ; light blue for core surfaces
red=ncolors-1 ; red for separatrix and x-cm SOL
lred=(ncolors-1)*0.9 ; light(er) red for SOL half-centimeter
plot,lim[0,0:lim_points-1],lim[1,0:lim_points-1],$
yrange=yrange/100.0,xrange=xrange/100.0, $
thick=5,title=tl,xtitle='R Meters', $
ytitle='Z Meters',xgridstyle=1,ygridstyle=1, $
xticklen=1,yticklen=1,/isotropic,font=0, $
position=[0.05,0.05,0.85,0.95]
levels = (findgen(100)+1)/100*(lrange[1] - lrange[0]) + lrange[0]
if pixels eq 0 then begin
    contour,csol,img_r,img_z, $
    /fill, /overplot, xstyle=4, ystyle=4, $
    levels=levels
endif else begin
    for i = 0,imxlen-1 do begin
        rpos = img_r[i]
        for j = 0,imylen-1 do begin
            zpos = img_z[j]
            rvec=[rpos-delr,rpos+delr,rpos+delr,rpos-delr]
            zvec=[zpos-delz,zpos-delz,zpos+delz,zpos+delz]
            if csol[i,j] gt 5 then polyfill,rvec,zvec,color=csol[i,j]
        endfor
    endfor
endelse
if snowflake eq 1 then begin
    contour,rmajor(0:ggg.mw-1,0:ggg.mh-1),psin_r,psin_z,/noerase,/overplot,$
    levels=[-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,-0.004,-0.003,-0.002,-0.001,0.0,0.001,0.002,0.003,0.004,0.005,0.01,0.02,0.03,0.04,0.05], $
    c_thick=[1,2,1,2,1,1,1,1,1,1,3,1,1,1,1,1,2,1,2,1],$
    c_colors = [lblue,blue,lblue,blue,lblue,lblue,lblue,lblue,lblue,lblue,red,lred,lred,lred,lred,lred,lred,red,lred,red,lred]
endif else begin
    contour,rmajor[0:ggg.mw-1,0:ggg.mh-1],psin_r,psin_z,/noerase,/overplot,$
    levels=[-0.04,-0.02,0.0,0.02,0.04], $
    c_thick=[1,1,2,1,1],$
    c_colors = [blue,blue,red,red,red]
endelse
oplot,lim[0,0:lim_points-1],lim[1,0:lim_points-1],$
thick=1,color=200
print,min(levels),max(levels)
tcolorbar,/vertical,position=[0.95,0.00,0.98,0.90],range=crange,format='(g8.3)',charsize=1.0,NCOLORS=ncolors
empty
if nographics eq 0 then begin
    set_plot,'X'
    device,decomposed=1
    !p.background=savebg
    !p.color=savecol
    loadct,0,/silent
endif
end
