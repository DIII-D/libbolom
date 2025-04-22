pro linuxtest
common xptcom,nextwindow
nextwindow = 0
nchans = long(0);
shot = long(186961)
shot = long(202730)
stime = 2000.0
xlen = long(0);
ylen = long(0);
DIRBOLO="/fusion/projects/diagnostics/bolometers/libbolom/"
DIRBOLO="./"
lib = DIRBOLO+'./linux64/libbolom6565.linux64.so'
x = call_external(lib,'bolom_sizes',nchans,xlen,ylen)
x = call_external(lib,'bolom_debug',250L)


print,nchans,xlen,ylen

maxchans = nchans
efit = fltarr(65,65)
efit(*) = 1.0
proj = fltarr(nchans)




g = readg(shot,stime,mode='mdsplus')
aeq = reada(shot,stime,mode='mdsplus')
rimage = efit
zimage = efit
xarr=g.rgrid1+ findgen(g.mw)        *g.xdim/(g.mw-1)
yarr=g.zmid  +(findgen(g.mh)-(g.mh-1)/2)*g.zdim/(g.mh-1)
fluximage=(g.psirz(0:g.mw-1,0:g.mh-1)-g.ssimag)/(g.ssibry-g.ssimag)

  for i=0,g.mw-1 do begin
    rimage(i,*)=xarr(i)
  endfor
  for i=0,g.mh-1 do begin
    zimage(*,i)=yarr(i)
  endfor

  rhoimage  =rho_rz(g,rimage,zimage)
  rhoimage  =REFORM(rhoimage,g.mw,g.mh,/overwrite)
  fitimage  =efit
  errorimage=fltarr(g.mw,g.mh)

  testimage = abs(max(rhoimage) - rhoimage)
x = call_external(lib,'bolom_proj',shot,testimage,proj)
window,0, xsize=400,ysize=800
plot,proj
tvscl,congrid(testimage,400,800)
origproj = proj
efit = fltarr(65,65)
window,1, xsize=400,ysize=800
x = call_external(lib,'bolom_bproj',shot,efit,proj)
tvscl,congrid(efit,400,800)
origimage = testimage


insigma = fltarr(n_elements(proj)) + 0.05 * max(proj)
nchans = long(21)
chans=long([1,2,3,4,15,16,17,18,19,20,21,22,39,40,41,42,43,44,45,46,47])
nkpsi = long(6)
kpsi = [0.0,0.5,0.9,1.0,1.06,1.12]
tension= 10.0
kmax =long(20)   ;max number of knot points
rmaxis = g.rmaxis*100.0
zmaxis = g.zmaxis*100.0
rsep1  = aeq.d.rxpt1
zsep1  = aeq.d.zxpt1
rsep2  = aeq.d.rxpt2
zsep2  = aeq.d.zxpt2

chi = 0.0
in = proj
fitproj = fltarr(maxchans)
fit = fltarr(maxchans)

status =call_external(lib,'bolom_core_fit',   $
                shot, chans, nchans, in, insigma, kpsi,          $
                nkpsi, tension, rhoimage,                       $
                rmaxis, zmaxis, rsep1, zsep1, rsep2, zsep2,     $
                fit, fitimage, errorimage, chi)

corefitimage = fitimage
window,2, xsize=400,ysize=800
tvscl,congrid(fitimage,400,800)
corefit = fit
totalfit = fitimage

nchans = long(48)
chans=long([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,$
                        20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,$
                        37,38,39,40,41,42,43,44,45,46,47])
sm_horiz_d=float(20.0)*1.0e-2
sm_vert_d=float(20.0)*1.0e-2
sm_horiz_b=float(60.0)*1.0e-2
sm_vert_b=float(60.0)*1.0e-2
z_inl=float(0.0)
z_outl=float(-70.0)
flux_inl=float(1.06)
flux_outl=float(1.06)
flux_corl=float(0.98)
flux_privl=float(0.995)
chi = 0.0



divproj = in - fit
divlimage = fluximage * 0.0

status = call_external(lib,'bolom_lower_fit_cells',$
        shot,chans,nchans,divproj,insigma,z_inl,z_outl,flux_privl,$
        flux_corl,flux_inl,flux_outl,sm_horiz_d,sm_vert_d,sm_horiz_b,sm_vert_b,$
        fluximage,rmaxis,zmaxis,rsep1,zsep1,rsep2,zsep2,fit,divlimage,chi)

lowerfit = fit
lowerfitimage = divlimage
window,3, xsize=400,ysize=800
tvscl,congrid(divlimage,400,800)
totalfit = totalfit + divlimage

divproj = origproj - corefit - lowerfit
divproj = origproj

n_chans=long(48)
chans=long([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,$
                        20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,$
                        37,38,39,40,41,42,43,44,45,46,47])
sm_horiz_d=float(20.0)
sm_vert_d=float(20.0)
sm_horiz_b=float(60.0)
sm_vert_b=float(60.0)
z_inu=float(0.0)
z_outu=float(70.0)
flux_inu=float(1.06)
flux_outu=float(1.06)
flux_coru=float(0.98)
flux_privu=float(0.995)
divuimage = fluximage * 0.0

x = call_external(lib,'bolom_upper_fit_cells',$
        shot,chans,nchans,divproj,insigma,z_inu,z_outu,flux_privu,$
        flux_coru,flux_inu,flux_outu,sm_horiz_d,sm_vert_d,sm_horiz_b,sm_vert_b,$
        fluximage,rmaxis,zmaxis,rsep1,zsep1,rsep2,zsep2,fit,divuimage,chi)



upperfit = fit
upperfitimage = divuimage
window,4, xsize=400,ysize=800
tvscl,congrid(divuimage,400,800)
totalfit = totalfit + divuimage
window,5, xsize=400,ysize=800
tvscl,congrid(totalfit,400,800)
range=[0,2.0]
; the next two commands require "module load tangtv"
;efitplot,shot,stime,totalfit,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,range=range
;efitplot,shot,stime,testimage,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,range=range


save

  return

  end

