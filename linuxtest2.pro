pro linuxtest2
common xptcom,nextwindow
nextwindow = 0
nchans = long(0);
shot = long(186961)
shot = long(202730)
stime = 2000.0
xlen = long(0);
ylen = long(0);
DIRBOLO="/fusion/projects/diagnostics/bolometers/libbolom/"
lib = DIRBOLO+'./linux64/libbolom6565.linux64.so'
x = call_external(lib,'bolom_sizes',nchans,xlen,ylen)
x = call_external(lib,'bolom_debug',250L)
print,nchans,xlen,ylen

efit = fltarr(65,65)
efit(*) = 1.0
proj = fltarr(nchans)

x = call_external(lib,'bolom_proj',shot,efit,proj)
plot,proj
efit = fltarr(65,65)
x = call_external(lib,'bolom_bproj',shot,efit,proj)
tvscl,congrid(efit,400,800)


chmax=long(nchans)  ;max number of channels for fitting
kmax =long(20)   ;max number of knot points
core ={core,$                                   ;core fit parameters
    n_chans:chmax,chans:lindgen(chmax),$            ;channels used for core
    n_krho:10L,krho:fltarr(kmax),$               ;rho knot points
    tension:1.0}

core.krho[0:core.n_krho-1] = findgen(core.n_krho) / (core.n_krho-1) 
insigma = fltarr(n_elements(proj)) + 0.05 * max(proj)

g = readg(shot,stime,mode='mdsplus')
aeq = reada(shot,stime,mode='mdsplus')
rimage = efit
zimage = efit
xarr=g.rgrid1+ findgen(g.mw)        *g.xdim/(g.mw-1)
yarr=g.zmid  +(findgen(g.mh)-(g.mh-1)/2)*g.zdim/(g.mh-1)

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

nchans = core.n_chans
chans  = core.chans(0:nchans-1)
nkpsi  = core.n_krho
kpsi   = core.krho(0:nkpsi-1)
tension= core.tension
rmaxis = g.rmaxis*100.0
zmaxis = g.zmaxis*100.0
rsep1  = aeq.d.rxpt1
zsep1  = aeq.d.zxpt1
rsep2  = aeq.d.rxpt2
zsep2  = aeq.d.zxpt2

chi = 0.0
in = proj
fitproj = fltarr(nchans)
fit = fltarr(nchans)

status =call_external(lib,'bolom_core_fit',   $
                shot, chans, nchans, in, insigma, kpsi,          $
                nkpsi, tension, rhoimage,                       $
                rmaxis, zmaxis, rsep1, zsep1, rsep2, zsep2,     $
                fit, fitimage, errorimage, chi)
print,'here'
print,min(fitimage),max(fitimage)
tvscl,congrid(fitimage,400,800)
save

range = [0.0,2.0]
; the next two commands require "module load tangtv"
;efitplot,shot,stime,fitimage,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,range=range
;efitplot,shot,stime,efit,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,range=range
  return
  end

