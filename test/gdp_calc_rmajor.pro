FUNCTION gdp_calc_rmajor,PSIN,a,g
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; C A L C _ R M A J O R
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Routine : CALC_RMAJOR
; Descrip : Calculate Rmajor given normalized psi values of a disgnostic.
;           Rmajor is calculated by splie-fitting the given normalized psi to
;           the efit grid normalized psi @ the z mag axis at R > r mag axis,
;           then subtracting R(normizedpsi=1.0), which is R sep.
;           Prior to this using this function, you must have first read in
;           the efitdata using READ_EFIT,shot,time.  This will fill in
;           the common block defined in EFITCOMMON.PRO properly.
; Author  : Michael D. Brown, LLNL, 10/09/92
;           Modified 1/29/93, MDB, No longer used AUTOCOM.PRO.
;           modified for new read efit routines 17Feb98 (GDP)
;
; Parameters:
;
; PSIN : Specify the normalized psi values for the calculation.
; a    : Structure containing A0 parameters
; g    : Structure containing G0 parameters

if n_elements(a.shot) eq 0 or n_elements(g.time) eq 0 then begin
  message,/cont,'You must first read in EFIT data (READ_EFIT,shot,time)'
  return,0
endif

sz = size(psin)

i_rout=where(g.r(0:g.mw-1) ge g.rmaxis)     ; r indexes outside the mag axis.
nrout = n_elements(i_rout)

; Select R(outer) and normalized psi(outer).
refit_outer=g.r(i_rout)         ; r values outside the mag axis.
;OLD without structures
;psin_ax_outer=calcpsi(refit_outer,g.zmaxis,a,g,/grid,/norm)
psin_ax_outer=gdp_calcpsi(refit_outer,g.zmaxis,a,g,/grid,/norm)
; psi @ z=zmaxis outside.

; qa that the psi values are of increasing order and do the calculation.
if total(sort(psin_ax_outer) eq indgen(nrout)) eq nrout then begin
  ; calculate r at the sep.
  rsep = float(interpol(refit_outer,psin_ax_outer,[1.0]),0)
  ; calculate r's for given psi values.
  rmajor=interpol(refit_outer,psin_ax_outer,psin)
  ; final result -- subtract of rsep offset. 
  rmajorsep=rmajor-rsep
  ; reform back to 2-dim array if necessary.
  if sz(0) eq 2 then rmajorsep=reform(temporary(rmajorsep),sz(1),sz(2))
endif else begin
  rmajorsep=fltarr(nrout)
  print,'(calc_rmajor) error: psin_ax_outer not increasing in value.'
endelse

RETURN,rmajorsep
END


