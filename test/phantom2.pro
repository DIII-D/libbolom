function phantom2_parms
parms = [[.1,.65,.5,-0.05,-0.05,0],$
        [.1,.65,.5,-0.05,-0.05,0],$
        [.1,.65,.5,-0.05,-0.05,0],$
        [.1,.65,.5,-0.05,-0.05,0],$
        [10,.22,.04,-0.5,-0.75,15]]
for i=0,3 do begin
      parms[1,i] = parms[1,i]/(i+1)
      parms[2,i] = parms[2,i]/(i+1)
endfor
return,parms
end

pro phantom2
common xptcom,nextwindow
nextwindow = 0
nchans = long(0);
shot = long(202934)
stime = 1000.0
xlen = long(0);
ylen = long(0);
DIRBOLO="/fusion/projects/diagnostics/bolometers/libbolom/"
DIRBOLO="../"
lib = DIRBOLO+'./linux64/libbolom6565.linux64.so'
x = call_external(lib,'bolom_sizes',nchans,xlen,ylen)
x = call_external(lib,'bolom_debug',250L)

chmax=long(48)  ;max number of channels for fitting
kmax =long(20)   ;max number of knot points
core ={core,$                                   ;core fit parameters
    n_chans:0L,chans:lonarr(chmax),$            ;channels used for core
    n_krho:0L,krho:fltarr(kmax),$               ;rho knot points
    tension:0.0}                                ;core spline tension

div={div,$                                      ;double null fit parameters
    n_chans:0L,chans:lonarr(chmax),$            ;channels for doublenull fit
    sm_horiz_d:0.0,$                            ;horizontal smoothing divertor
    sm_vert_d:0.0,$                             ;vertical smoothing divertor
    sm_horiz_b:0.0,$                            ;horizontal smoothing boundary
    sm_vert_b:0.0,$                             ;vertical smoothing boundary
    z_inl:0.0,$                                ;max z for lowin div rad
    z_outl:0.0,$                                ;max z for lowout div rad
    z_inu:0.0,$                                 ;min z for upin div rad
    z_outu:0.0,$                                ;min z for upout div rad
    flux_inl:0.0,$                              ;max flux lowinner div
    flux_inu:0.0,$                              ;max flux upinner div
    flux_outl:0.0,$                             ;max flux lowouter div
    flux_outu:0.0,$                             ;max flux upouter div
    flux_corl:0.0,$                             ;min flux lower core
    flux_coru:0.0,$                             ;min flux upper core
    flux_privl:0.0,$                            ;min flux lower private
    flux_privu:0.0}                             ;min flux upper private

param={param,core:core,div:div}

lsn_def=param
usn_def=param
dnd_def=param

lsn_def.core.n_chans=16
lsn_def.core.chans(0:15)=[1,2,3,14,15,16,17,18,19,20,38,39,40,41,42,47]
lsn_def.core.chans(16:*)=0
lsn_def.core.n_krho=5
lsn_def.core.krho(0:4)=[0.0,0.9,1.0,1.06,1.12]
lsn_def.core.krho(5:*)=0
lsn_def.core.tension=float(1.0)
lsn_def.div.n_chans=28
lsn_def.div.chans(0:27)=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,24,25,26,$
                         27,28,29,30,31,32,33,34,35]
lsn_def.div.chans(28:*)=0
lsn_def.div.sm_horiz_d=float(20.0)
lsn_def.div.sm_vert_d=float(20.0)
lsn_def.div.sm_horiz_b=float(60.0)
lsn_def.div.sm_vert_b=float(60.0)
lsn_def.div.z_inl=float(0.0)
lsn_def.div.z_outl=float(-70.0)
lsn_def.div.flux_inl=float(1.10)
lsn_def.div.flux_outl=float(1.10)
lsn_def.div.flux_corl=float(0.98)
lsn_def.div.flux_privl=float(0.99)

usn_def.core.n_chans=26
usn_def.core.chans(0:25)=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,24,25,26,27,28,29,$
                         30,31,32,33,34,35]
usn_def.core.chans(26:*)=0.0
usn_def.core.n_krho=5
usn_def.core.krho(0:4)=[0.0,0.9,1.0,1.06,1.12]
usn_def.core.krho(5:*)=0
usn_def.core.tension=float(1.0)
usn_def.div.n_chans=23
usn_def.div.chans(0:22)=[14,15,16,17,18,19,20,21,22,23,34,35,36,$
                        38,39,40,41,42,43,44,45,46,47]
usn_def.div.chans(23:*)=0
usn_def.div.sm_horiz_d=float(20.0)
usn_def.div.sm_vert_d=float(20.0)
usn_def.div.sm_horiz_b=float(60.0)
usn_def.div.sm_vert_b=float(60.0)
usn_def.div.z_inu=float(0.0)
usn_def.div.z_outu=float(70.0)
usn_def.div.flux_inu=float(1.10)
usn_def.div.flux_outu=float(1.10)
usn_def.div.flux_coru=float(0.98)
usn_def.div.flux_privu=float(0.99)

dnd_def.core.n_chans=17
dnd_def.core.chans(0:16)=[1,2,3,4,5,13,14,15,16,17,32,33,34,35,36,38,47]
dnd_def.core.chans(17:*)=0
dnd_def.core.n_krho=5
dnd_def.core.krho(0:4)=[0.0,0.9,1.0,1.06,1.12]
dnd_def.core.krho(5:*)=0
dnd_def.core.tension=float(1.0)
dnd_def.div.n_chans=47
dnd_def.div.chans(0:46)=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,$
                        20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,$
                        38,39,40,41,42,43,44,45,46,47]
dnd_def.div.chans(47)=0
dnd_def.div.sm_horiz_d=float(20.0)
dnd_def.div.sm_vert_d=float(20.0)
dnd_def.div.sm_horiz_b=float(60.0)
dnd_def.div.sm_vert_b=float(60.0)
dnd_def.div.z_inl=float(0.0)
dnd_def.div.z_outl=float(-70.0)
dnd_def.div.z_inu=float(0.0)
dnd_def.div.z_outu=float(70.0)
dnd_def.div.flux_inl=float(1.10)
dnd_def.div.flux_inu=float(1.10)
dnd_def.div.flux_outl=float(1.10)
dnd_def.div.flux_outu=float(1.10)
dnd_def.div.flux_corl=float(0.98)
dnd_def.div.flux_coru=float(0.98)
dnd_def.div.flux_privl=float(0.99)
dnd_def.div.flux_privu=float(0.99)

core = lsn_def.core
lsn = lsn_def.div
usn = []


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

print,g.mw,g.mh
  rhoimage  =rho_rz(g,rimage,zimage)
  rhoimage  =REFORM(rhoimage,g.mw,g.mh,/overwrite)
  fitimage  =efit
  errorimage=fltarr(g.mw,g.mh)

  bigimage = phantom(650,650,ellipses=phantom2_parms())
  testimage = congrid(bigimage,65,65)

efitplot,shot,stime,testimage,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,range=range,tl='Phantom 2'
status = call_external(lib,'bolom_proj',shot,testimage,proj)

save,file='phantom2.dat'

  return

  end

