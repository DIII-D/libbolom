#include <stdio.h>
#include "bolom.h"
extern int debug;

bolom_core_fit_python_wrapper(shot,chans,nchans,inproj,sigma,kpsi,nkpsi,tension,efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,errorimage,fret)
    float *inproj,*sigma,*kpsi,*fitproj,*fret,*tension;
    float *efit,*fitimage,*errorimage;
    int *chans,*nchans,*nkpsi,*shot;
    float *rax,*zax,*rxpt1,*zxpt1,*rxpt2,*zxpt2;
{
  // Wraper by SRH so that the same routines can be called from OMFITprofiles using python
  // because of difficulties using the existing wrapper for IDL
  bolomfit_core(*shot,chans,*nchans,inproj,sigma,kpsi,*nkpsi,tension,efit, 
		rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,errorimage,fret); 
  return(0);
}
