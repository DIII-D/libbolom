; Translated from:
;## Copyright (C) 2010 Alex Opie <lx_op@orcon.net.nz>
;##
;## This program is free software; you can redistribute it and/or modify it
;## under the terms of the GNU General Public License as published by
;## the Free Software Foundation; either version 3 of the License, or (at
;## your option) any later version.
;##
;## This program is distributed in the hope that it will be useful, but
;## WITHOUT ANY WARRANTY; without even the implied warranty of
;## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
;## General Public License for more details.
;##
;## You should have received a copy of the GNU General Public License
;## along with this program; see the file COPYING. If not, see
;## <http:
;
; by Bill Meyer <meyer8@llnl.gov> 3/25/2014
;
;
;
function shepp_logan_parms
return,[[ 2, .69, .92, 0, 0, 0], $
[-.98, .6624, .8740, 0, -.0184, 0], $
[-.02, .1100, .3100, .22, 0, -18], $
[-.02, .1600, .4100, -.22, 0, 18], $
[ .01, .2100, .2500, 0, .35, 0], $
[ .01, .0460, .0460, 0, .1, 0], $
[ .02, .0460, .0460, 0, -.1, 0], $
[ .01, .0460, .0230, -.08, -.605, 0], $
[ .01, .0230, .0230, 0, -.606, 0], $
[ .01, .0230, .0460, .06, -.605, 0]]
end
function mod_shepp_logan_parms
return,[[ 1, .69, .92, 0, 0, 0], $
[-.80, .6624, .8740, 0, -.0184, 0], $
[-.20, .1100, .3100, .22, 0, -18], $
[-.20, .1600, .4100, -.22, 0, 18], $
[ .10, .2100, .2500, 0, .35, 0], $
[ .10, .0460, .0460, 0, .1, 0], $
[ .10, .0460, .0460, 0, -.1, 0], $
[ .10, .0460, .0230, -.08, -.605, 0], $
[ .10, .0230, .0230, 0, -.606, 0], $
[ .10, .0230, .0460, .06, -.605, 0]]
end
function pos_shepp_logan_parms
return,[[ 1.8, .69, .92, 0, 0, 0], $
[ .00, .6624, .8740, 0, -.0184, 0], $
[ .60, .1100, .3100, .22, 0, -18], $
[ .60, .1600, .4100, -.22, 0, 18], $
[ .90, .2100, .2500, 0, .35, 0], $
[ .90, .0460, .0460, 0, .1, 0], $
[ .90, .0460, .0460, 0, -.1, 0], $
[ .90, .0460, .0230, -.08, -.605, 0], $
[ .90, .0230, .0230, 0, -.606, 0], $
[ .90, .0230, .0460, .06, -.605, 0]]
end
function flow_shepp_logan_parms
return,[[ 1.8, .69, .92, 0, 0, 0], $
[ .00, .6624, .8740, 0, -.0184, 0], $
[ -2.60, .1100, .3100, .22, 0, -18], $
[ .60, .1600, .4100, -.22, 0, 18], $
[ .90, .2100, .2500, 0, .35, 0], $
[ .90, .0460, .0460, 0, .1, 0], $
[ -1.90, .0460, .0460, 0, -.1, 0], $
[ .90, .0460, .0230, -.08, -.605, 0], $
[ .90, .0230, .0230, 0, -.606, 0], $
[ .90, .0230, .0460, .06, -.605, 0]]
end
function phantom,n,m,type=type,help=help,ellipses=ellipses
phan = fltarr(n,m)
if n_elements(type) eq 0 and isa(ellipses) eq 0 then type = 'Modified Pos Shepp-Logan'
if n_elements(type) eq 0 then type = 'none'
if n_elements(help) ne 0 then goto,help
SWITCH type of
    'Shepp-Logan': Begin
        ellipses = shepp_logan_parms()
        BREAK
    end
    'Modified Shepp-Logan': Begin
        ellipses = mod_shepp_logan_parms()
        BREAK
    end
    'Modified Pos Shepp-Logan': Begin
        ellipses = pos_shepp_logan_parms()
        BREAK
    end
    'Flow Shepp-Logan': Begin
        ellipses = flow_shepp_logan_parms()
        BREAK
    end
    else: Begin
        BREAK
    end
ENDSWITCH
s = size(ellipses)
t = (findgen(n)/(n-1.0) * 2.0) - 1.0
t = reform(t,n,1)
expand,t,n,m,xgrid
t = (findgen(n)/(n-1.0) * 2.0) - 1.0
t = transpose(t)
expand,t,n,m,ygrid
if s[2] eq 0 then begin
    print,'No ellipses found for type ',type
    return,phan
end
for k = 1,s[2] do begin
    ellip = k-1
    I = ellipses[0,ellip]
    a2 = ellipses[1,ellip]^2
    b2 = ellipses[2,ellip]^2
    x0 = ellipses[3,ellip]
    y0 = ellipses[4,ellip]
    phi = ellipses[5,ellip] * !pi / 180.0
    x = xgrid - x0
    y = ygrid - y0
    cos_p = cos(phi)
    sin_p = sin(phi)
    locs = (((x * cos_p + y * sin_p)^2) / a2 + ((y * cos_p - x * sin_p)^2) / b2)
    j = where(locs le 1,num)
    if num gt 0 then begin
        phan(j) = phan(j) + I
    endif
endfor
return,phan
help:
print," Create 2-D phantom for tomography testing. "
print," <variable> = phantom(columns,rows,type=type,help=help)"
print
print," Type is 'Shepp-Logan'"
print,"         'Modified Shepp-Logan' - higher contrast"
print,"         'Modified Pos Shepp-Logan' - higher contrast with only positive values"
print
return,phan
end
