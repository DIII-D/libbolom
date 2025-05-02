#include <stdio.h>
#include "bolom.h"
extern int debug;

bolom_debug(argc,argv)
int argc;
void *argv[];
{
    int *idebug;
    if(argc != 1){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_debug\",num) \n");
        fprintf(stderr,"       where:  num - if iterations modulo num equals zero some debug\n");
        fprintf(stderr,"                     info is printed, 0 to turn debug info off\n");
        fprintf(stderr,"\n");
        return(1);
    }
    idebug = argv[0];
    debug = *idebug;

    return(0);

}

bolom_sizes(argc,argv)
int argc;
void *argv[];
{
    int *x,*y,*ichans;
    if(argc != 3){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_sizes\",chans,\n");
        fprintf(stderr,"                                   image_x,image_y) \n");
        fprintf(stderr,"       where:  chans - number of input soft x-ray channels\n");
        fprintf(stderr,"               image_x - x dimension of solution image/cos terms\n");
        fprintf(stderr,"               image_y - y dimension of solution image/cos terms\n");
        fprintf(stderr,"               all are long integer\n");
        fprintf(stderr,"\n");
        return(1);
    }
	fprintf(stderr,"bolom sizes\n");
    ichans = argv[0];
    x = argv[1];
    y = argv[2];

    *ichans = CHANS;
    *x = XLEN;
    *y = YLEN;

    return(0);

}
#if defined(CORE)
bolom_core_fit(argc,argv)
int argc;
void *argv[];
{
    float *inproj,*sigma,*kpsi,*fitproj,*fret,*tension;
    float *efit,*fitimage,*errorimage;
    int *chans,*nchans,*nkpsi,*shot;
		float *rax,*zax,*rxpt1,*zxpt1,*rxpt2,*zxpt2;

    if(argc != 19){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_core_fit\",shot,chans,\n");
        fprintf(stderr,"                                   nchans,inproj,sigma,kpsi,nkpsi,tension,\n");
        fprintf(stderr,"                                   efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,\n");
        fprintf(stderr,"                                   fitproj,fitimage,errorimage,chi\n");
        fprintf(stderr,"       where:  shot - shot number of raw data determines which geometry matrix is used (long integer)\n");
        fprintf(stderr,"               chans - array of channel numbers to fit (long integer)\n");
        fprintf(stderr,"               nchans - number of channels in array chans (long integer)\n");
        fprintf(stderr,"               inproj - raw bolometer line integrated data in W (real)\n");
        fprintf(stderr,"               sigma - corresponding sigma values for raw data (real)\n");
        fprintf(stderr,"               kpsi - array of normalized psi values at which to iterate\n");
        fprintf(stderr,"                      the W/cm^3 for fit (real 0.0 - 1.0)\n");
        fprintf(stderr,"               nkpsi - number of kpsi values (long integer)\n");
        fprintf(stderr,"               tension - B-spline smoothing parameter (real)\n");
        fprintf(stderr,"               efit - 2d array of normalized psi values on efit 33x65\n");
        fprintf(stderr,"                      grid (real 0.0 - 1.0)\n");
        fprintf(stderr,"               rax,zax - coordinates of axis in cm\n");
        fprintf(stderr,"               rxpt1,zxpt1 - coordinates of lower xpt in cm\n");
        fprintf(stderr,"               rxpt2,zxpt2 - coordinates of upper xpt in cm\n");
        fprintf(stderr,"               fitproj - line integral of fit image (real)\n");
        fprintf(stderr,"               fitimage - 2d array of W/cm^3 bolometer power on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               error - 2d array of error values on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               fret - chi^2 returned\n");
        fprintf(stderr,"\n");
        return(1);
    }
    shot = argv[0];
    chans = argv[1];
    nchans = argv[2];
    inproj = argv[3];
    sigma = argv[4];
    kpsi = argv[5];
    nkpsi = argv[6];
    tension = argv[7];
    efit = argv[8];
    rax = argv[9];
    zax = argv[10];
    rxpt1 = argv[11];
    zxpt1 = argv[12];
    rxpt2 = argv[13];
    zxpt2 = argv[14];
    fitproj = argv[15];
    fitimage = argv[16];
		errorimage = argv[17];
    fret = argv[18];
    bolomfit_core(*shot,chans,*nchans,inproj,sigma,kpsi,*nkpsi,tension,efit,
        rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,errorimage,fret);

    return(0);
}
#endif
#if defined(LOWER)
bolom_lower_fit_cells(argc,argv)
int argc;
void *argv[];
{
    float *inproj,*sigma,*fitproj,*fret;
    float *efit,*fitimage;
    int *chans,*nchans,*shot;
		float *rax,*zax,*rxpt1,*zxpt1,*rxpt2,*zxpt2;
		float *zinner,*zouter,*minprivate,*mincore;
		float *maxinner,*maxouter,*drwtd,*dzwtd,*drwtb,*dzwtb;

    if(argc != 25){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_lower_fit_cells\",shot,chans,\n");
        fprintf(stderr,"                                   nchans,inproj,sigma,zinner,zouter,minprivate,\n");
        fprintf(stderr,"                                   mincore,maxinner,maxouter,drwtd,dzwtd,drwtb,dzwtb,\n");
        fprintf(stderr,"                                   efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,\n");
        fprintf(stderr,"                                   fitproj,fitimage,errorimage,chi\n");
        fprintf(stderr,"       where:  shot - shot number of raw data determines which geometry matrix is used (long integer)\n");
        fprintf(stderr,"               chans - array of channel numbers to fit (long integer)\n");
        fprintf(stderr,"               nchans - number of channels in array chans (long integer)\n");
        fprintf(stderr,"               inproj - raw bolometer line integrated data in W (real)\n");
        fprintf(stderr,"               sigma - corresponding sigma values for raw data (real)\n");
        fprintf(stderr,"              zinner - Z value of highest cell inside of Xpt (real, units cm)\n");
        fprintf(stderr,"             zouter - Z value of highest cell outside of Xpt (real, units cm)\n");
        fprintf(stderr,"          minprivate - minimum normalized Psi value of cell to be used located in private flux region\n");
        fprintf(stderr,"             mincore - minimum normalized Psi value of cell to be used located in core region\n");
        fprintf(stderr,"            maxinner - maximum normalized Psi value of cell to be used located inside of separatrix\n");
        fprintf(stderr,"            maxouter - maximum normalized Psi value of cell to be used located outside of separatrix\n");
        fprintf(stderr,"            drwtd - weighting factor for smoothness constraint in R direction\n");
        fprintf(stderr,"            dzwtd - weighting factor for smoothness constraint in Z direction\n");
        fprintf(stderr,"            drwtb - weighting factor for smoothness constraint in R direction\n");
        fprintf(stderr,"            dzwtb - weighting factor for smoothness constraint in Z direction\n");
        fprintf(stderr,"               efit - 2d array of normalized psi values on efit 33x65\n");
        fprintf(stderr,"                      grid (real 0.0 - 1.0)\n");
        fprintf(stderr,"               rax,zax - coordinates of axis in cm\n");
        fprintf(stderr,"               rxpt1,zxpt1 - coordinates of lower xpt in cm\n");
        fprintf(stderr,"               rxpt2,zxpt2 - coordinates of upper xpt in cm\n");
        fprintf(stderr,"               fitproj - line integral of fit image (real)\n");
        fprintf(stderr,"               fitimage - 2d array of W/cm^3 bolometer power on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               fret - chi^2 returned\n");
        fprintf(stderr,"\n");
        return(1);
    }
    shot = argv[0];
    chans = argv[1];
    nchans = argv[2];
    inproj = argv[3];
    sigma = argv[4];
    zinner = argv[5];
    zouter = argv[6];
    minprivate = argv[7];
    mincore = argv[8];
    maxinner = argv[9];
    maxouter = argv[10];
    drwtd = argv[11];
    dzwtd = argv[12];
    drwtb = argv[13];
    dzwtb = argv[14];
    efit = argv[15];
    rax = argv[16];
    zax = argv[17];
    rxpt1 = argv[18];
    zxpt1 = argv[19];
    rxpt2 = argv[20];
    zxpt2 = argv[21];
    fitproj = argv[22];
    fitimage = argv[23];
    fret = argv[24];
    bolomfit_lower_cells(*shot,chans,*nchans,inproj,sigma,zinner,zouter,
				minprivate,mincore,maxinner,maxouter,drwtd,dzwtd,drwtb,dzwtb,
        efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,fret);

    return(0);

}
#endif
#if defined(UPPER)
bolom_upper_fit_cells(argc,argv)
int argc;
void *argv[];
{
    float *inproj,*sigma,*fitproj,*fret;
    float *efit,*fitimage;
    int *chans,*nchans,*shot;
		float *rax,*zax,*rxpt1,*zxpt1,*rxpt2,*zxpt2;
		float *zinner,*zouter,*minprivate,*mincore;
		float *maxinner,*maxouter,*drwtd,*dzwtd,*drwtb,*dzwtb;

    if(argc != 25){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_upper_fit_cells\",shot,chans,\n");
        fprintf(stderr,"                                   nchans,inproj,sigma,zinner,zouter,minprivate,\n");
        fprintf(stderr,"                                   mincore,maxinner,maxouter,drweight,dzweight,\n");
        fprintf(stderr,"                                   efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,\n");
        fprintf(stderr,"                                   fitproj,fitimage,errorimage,chi\n");
        fprintf(stderr,"       where:  shot - shot number of raw data determines which geometry matrix is used (long integer)\n");
        fprintf(stderr,"               chans - array of channel numbers to fit (long integer)\n");
        fprintf(stderr,"               nchans - number of channels in array chans (long integer)\n");
        fprintf(stderr,"               inproj - raw bolometer line integrated data in W (real)\n");
        fprintf(stderr,"               sigma - corresponding sigma values for raw data (real)\n");
        fprintf(stderr,"              zinner - Z value of highest cell inside of Xpt (real, units cm)\n");
        fprintf(stderr,"             zouter - Z value of highest cell outside of Xpt (real, units cm)\n");
        fprintf(stderr,"          minprivate - minimum normalized Psi value of cell to be used located in private flux region\n");
        fprintf(stderr,"             mincore - minimum normalized Psi value of cell to be used located in core region\n");
        fprintf(stderr,"            maxinner - maximum normalized Psi value of cell to be used located inside of separatrix\n");
        fprintf(stderr,"            maxouter - maximum normalized Psi value of cell to be used located outside of separatrix\n");
        fprintf(stderr,"            drwtd - weighting factor for smoothness constraint in R direction\n");
        fprintf(stderr,"            dzwtd - weighting factor for smoothness constraint in Z direction\n");
        fprintf(stderr,"            drwtb - weighting factor for smoothness constraint in R direction\n");
        fprintf(stderr,"            dzwtb - weighting factor for smoothness constraint in Z direction\n");
        fprintf(stderr,"               efit - 2d array of normalized psi values on efit 33x65\n");
        fprintf(stderr,"                      grid (real 0.0 - 1.0)\n");
        fprintf(stderr,"               rax,zax - coordinates of axis in cm\n");
        fprintf(stderr,"               rxpt1,zxpt1 - coordinates of lower xpt in cm\n");
        fprintf(stderr,"               rxpt2,zxpt2 - coordinates of upper xpt in cm\n");
        fprintf(stderr,"               fitproj - line integral of fit image (real)\n");
        fprintf(stderr,"               fitimage - 2d array of W/cm^3 bolometer power on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               fret - chi^2 returned\n");
        fprintf(stderr,"\n");
        return(1);
    }
    shot = argv[0];
    chans = argv[1];
    nchans = argv[2];
    inproj = argv[3];
    sigma = argv[4];
    zinner = argv[5];
    zouter = argv[6];
    minprivate = argv[7];
    mincore = argv[8];
    maxinner = argv[9];
    maxouter = argv[10];
    drwtd = argv[11];
    dzwtd = argv[12];
    drwtb = argv[13];
    dzwtb = argv[14];
    efit = argv[15];
    rax = argv[16];
    zax = argv[17];
    rxpt1 = argv[18];
    zxpt1 = argv[19];
    rxpt2 = argv[20];
    zxpt2 = argv[21];
    fitproj = argv[22];
    fitimage = argv[23];
    fret = argv[24];
    bolomfit_upper_cells(*shot,chans,*nchans,inproj,sigma,zinner,zouter,
				minprivate,mincore,maxinner,maxouter,drwtd,dzwtd,drwtb,dzwtb,
        efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,fret);

    return(0);

}
#endif
bolom_proj(argc,argv)
int argc;
void *argv[];
{
		int *shot;
    float *fitproj;
    float *fitimage;

    if(argc != 3){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_proj\",shot,image,proj)\n");
        fprintf(stderr,"       where:  shot - shot number determines which geometry matrix is used (long integer)\n");
        fprintf(stderr,"               image - 2d array of W/cm^3 bolometer power on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               fitproj - line integral of fit image (real)\n");
        fprintf(stderr,"\n");
        return(1);
    }
		shot = argv[0];
    fitimage = argv[1];
    fitproj = argv[2];
    bolomproj(*shot,fitimage,fitproj);

    return(0);

}
bolom_bproj(argc,argv)
int argc;
void *argv[];
{
		int *shot;
    float *fitproj;
    float *fitimage;

    if(argc != 3){
        fprintf(stderr,"\n");
        fprintf(stderr,"IDL Usage   : stat = call_external(\"libbolom.so\",\"bolom_bproj\",shot,image,proj)\n");
        fprintf(stderr,"       where:  shot - shot number determines which geometry matrix is used (long integer)\n");
        fprintf(stderr,"               image - 2d array of W/cm^3 bolometer power on efit\n");
        fprintf(stderr,"                          grid (real)\n");
        fprintf(stderr,"               fitproj - line integral of fit image (real)\n");
        fprintf(stderr,"\n");
        return(1);
    }
		shot = argv[0];
    fitimage = argv[1];
    fitproj = argv[2];
    bolombackproj(*shot,fitimage,fitproj);

    return(0);

}

