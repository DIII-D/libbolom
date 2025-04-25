function gdp_calcpsi,r,z,a,g,grid=grid,norm=norm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; routine: calcpsi
; author:  Michael D. Brown, LLNL
; date:    01/25/93
; desc:    Given r and z (or array of r and array of z) calculate and return
;          a corresponding scalar or array of psi  OR normalize psi (NORM=1).
;          Note that the eqdsk files must have been previously read in by
;          the routine READ_EFIT, which fills in the common block defined
;          in the file EFITCOMMON.PRO.
; Parameters:
;          r      - scalar or array of radius values, in meters.
;          z      - scalar or array of vertical z values, in meters.
;          a      - structure containing A0 parameters
;          g      - structure containing G0 parameters
; Key words:
;          GRID   - Specify /GRID or GRID=1 to specify that the calculated
;                   psi normalized be returned as if r and z represented
;                   a grid of r's and z's rather than just ordered pairs.
;          NORM   - Specify /NORM or NORM=1 to return NORMALIZED psi values.
; Returned:
;          Scalar, vector, or matrix of psi values, depending
;          on input types and keywords specified.
;         
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



; qa that efit data has been read
if n_elements(a.shot) eq 0 or n_elements(g.shot) eq 0 then begin
  message,/cont,'You must first read in the efit data with READ_EFIT,shot,time'
  return,[0.0]
endif

; qa input parameters.
if n_elements(grid) eq 0 then grid=0
if n_elements(norm) eq 0 then norm=0
if (n_elements(r) ne n_elements(z) and grid eq 0) or  $
    n_elements(r) eq 0 or n_elements(z) eq 0 then return,[0.0]

; calculate psi at given r's and z's.
temppsi = bicubicspline(g.psirz(0:g.mw-1,0:g.mh-1),g.r(0:g.mw-1), $
                        g.z(0:g.mh-1),r,z,grid=grid,ierr=ierr)
if ierr ne 0 then temppsi=[0.0]

; convert to scaler or leave as vector, as required by inputs
rinfo=size(r) & zinfo=size(z)
if rinfo(0) eq 0 and zinfo(0) eq 0 then temppsi=temppsi(0)

; normalize psi values if required.
if norm ne 0 then temppsi=(temppsi-g.ssimag)/(g.ssibry-g.ssimag)
return,temppsi
END
