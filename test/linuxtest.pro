pro linuxtest,savefile,mse,production=production,referencefile=referencefile,saveref=saveref
restore,savefile
common xptcom,nextwindow
nextwindow = 0
if isa(production) eq 0 then begin
    DIRBOLO="../"
endif else begin
    DIRBOLO=getenv("DIRBOLO")
endelse
lib = DIRBOLO+'./linux64/libbolom6565.linux64.so'

testplots = 0 ; turn on to debug test
print,'Using library :',lib
maxchans = nchans
efit = fltarr(65,65)
efit(*) = 1.0




origproj = proj
efit = fltarr(65,65)
x = call_external(lib,'bolom_bproj',shot,efit,proj)
if testplots eq 1 then begin
   window,nextwindow, xsize=800,ysize=400,title="Channels"
   nextwindow = nextwindow+1
   plot,proj
   window,nextwindow, xsize=400,ysize=800,title="Phantom"
   nextwindow = nextwindow+1
   tvscl,congrid(testimage,400,800)
   window,nextwindow, xsize=400,ysize=800,title="Phantom Backproject"
   nextwindow = nextwindow+1
   tvscl,congrid(efit,400,800)
endif
origimage = testimage


insigma = fltarr(n_elements(proj)) + 0.05 * max(proj)
nchans=core.n_chans
chans=core.chans(0:nchans-1)
nkpsi=core.n_krho
kpsi=core.krho(0:nkpsi-1)
tension=core.tension
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
if testplots eq 1 then begin
   window,nextwindow, xsize=400,ysize=800,title="Core Fit"
   nextwindow = nextwindow+1
   tvscl,congrid(fitimage,400,800)
endif
corefit = fit
totalfit = fitimage


sm_horiz_d=float(20.0)*1.0e-2
sm_vert_d=float(20.0)*1.0e-2
sm_horiz_b=float(60.0)*1.0e-2
sm_vert_b=float(60.0)*1.0e-2
z_inl=float(0.0)
z_outl=float(-70.0)
flux_inl=float(1.10)
flux_outl=float(1.10)
flux_corl=float(0.98)
flux_privl=float(0.99)
chi = 0.0



divproj = in - fit
divlimage = fluximage * 0.0

if isa(lsn) ne 0 then begin
nchans    = lsn.n_chans
chans= lsn.chans(0:nchans-1)
status = call_external(lib,'bolom_lower_fit_cells',$
        shot,chans,nchans,divproj,insigma,z_inl,z_outl,flux_privl,$
        flux_corl,flux_inl,flux_outl,sm_horiz_d,sm_vert_d,sm_horiz_b,sm_vert_b,$
        fluximage,rmaxis,zmaxis,rsep1,zsep1,rsep2,zsep2,fit,divlimage,chi)
endif

lowerfit = fit
lowerfitimage = divlimage
if testplots eq 1 then begin
   window,nextwindow, xsize=400,ysize=800,title="Lower Divertor Fit"
   nextwindow = nextwindow+1
   tvscl,congrid(divlimage,400,800)
endif
totalfit = totalfit + divlimage

divproj = origproj - corefit - lowerfit
divproj = origproj

sm_horiz_d=float(20.0)*1.0e-2
sm_vert_d=float(20.0)*1.0e-2
sm_horiz_b=float(60.0)*1.0e-2
sm_vert_b=float(60.0)*1.0e-2
z_inu=float(0.0)
z_outu=float(70.0)
flux_inu=float(1.10)
flux_outu=float(1.10)
flux_coru=float(0.98)
flux_privu=float(0.99)
divuimage = fluximage * 0.0

if isa(usn) ne 0 then begin
nchans    = usn.n_chans
chans= usn.chans(0:nchans-1)
x = call_external(lib,'bolom_upper_fit_cells',$
        shot,chans,nchans,divproj,insigma,z_inu,z_outu,flux_privu,$
        flux_coru,flux_inu,flux_outu,sm_horiz_d,sm_vert_d,sm_horiz_b,sm_vert_b,$
        fluximage,rmaxis,zmaxis,rsep1,zsep1,rsep2,zsep2,fit,divuimage,chi)
endif





upperfit = fit
upperfitimage = divuimage
totalfit = totalfit + divuimage
if testplots eq 1 then begin
   window,nextwindow, xsize=400,ysize=800,title="Upper Divertor Fit"
   nextwindow = nextwindow+1
   tvscl,congrid(divuimage,400,800)
   window,nextwindow, xsize=400,ysize=800,title="Total Fit"
   nextwindow = nextwindow+1
   tvscl,congrid(totalfit,400,800)
   ; the next two commands require "module load tangtv"
   range = [min(testimage),max(testimage)]
   efitplot,shot,stime,totalfit,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,tl='Fit',range=range
   efitplot,shot,stime,testimage,matrmin=84.0,matrmax=254.0,matzmin=-160.0,matzmax=160.0,tl='Phantom',range=range
endif


mse = 0
reference = totalfit
if isa(saveref) ne 0 then begin
   if isa(referencefile) ne 0 then begin
        save,file=referencefile,reference
   endif else begin
        save,file='reference.dat',reference
   endelse
endif else begin
   if isa(referencefile) ne 0 then begin
      restore,referencefile
      ; Calculate the squared differences
      squared_diff = (totalfit - reference)^2
      
      ; Compute the mean of these squared differences
      mse = total(squared_diff) / n_elements(totalfit)
   endif 
endelse

  return

  end

