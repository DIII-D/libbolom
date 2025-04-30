#define EMIN 0.0
#define EMAX 10.0

#include "bolom.h"
extern int debug;
#include <math.h>
#include <stdlib.h>
#include <nrutil.h>
#include <signal.h>
#include <stdio.h>
#include "boundext.h"
extern float *gmatrix;
extern int *pow_its;
extern int *gpos;
extern int gmatrix_len;
extern int gpos_len;
#include <setjmp.h>
static jmp_buf sjbuf ;
static int funcevals;

#define PROJ (gpos[k])
#define CELL (gpos[j])

static void (*prev_handler)();
static float *undump;
static float *indata;
static float *insigma;
static int *fitchans;
static int nfitchans;
static float *rfitpts,*zfitpts;
static int nrfitpts,nzfitpts;
static int done;
float psi_norm[XLEN*YLEN];




bolomfit_upper(shot,chans,nchans,inproj,sigma,r,nr,z,nz,xtens,ytens,efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,ftol,fret,iter)
int shot;
float *inproj,*sigma,*r,*z,*fitproj,ftol,*fret;
float xtens,ytens;
image fitimage;
image efit;
int *chans,nchans,nr,nz,*iter;
float rax,zax,rxpt1,zxpt1,rxpt2,zxpt2;

{

    float spuppereval();
    float *p,**xi;
    int dumpum_spupper();
    FILE *f;
    int i,j,k;
    int n;

    struct sigaction vec,ovec;
    funcevals = 1;

    fitchans = chans;
    nfitchans = nchans;
    indata = inproj;
    insigma = sigma;
    rfitpts = r;
    zfitpts = z;
    nrfitpts = nr;
    nzfitpts = nz;
    n = nr * nz;
    done = 0;
    efitgeom(efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,psi_norm);
		if(!bolom_set_gmatrix(shot))return(0);
    setjmp(sjbuf);
    if(done)return;
    done = 1;

#if defined(TRAPS)
    vec.sa_handler = dumpum_spupper;
    sigfillset(&vec.sa_mask);
    vec.sa_flags = 0;
    sigaction(SIGINT,&vec,&ovec);
    prev_handler = ovec.sa_handler;
#endif

    /*write_projdat("upperorig.dat",inproj);*/

    p = vector(1,n);
    undump = p;
    for(i = 1; i <= n;++i)p[i] = 0.;
    p[1] = 1e-4;
    xi = matrix(1,n,1,n);

    for(i = 1; i <= n ; ++i)
        memset(&xi[i][1],0,n * sizeof(float));

    for(i = 1; i <= n ; ++i){
        xi[i][i] = 1e-4;
    }

    if((f = fopen("uppercoefs.dat","r")) != NULL){
        for(i = 1;i <= n; ++i){
            if(fscanf(f,"%g\n",&undump[i]) != 1)break;
        }
        fclose(f);
    }



    memset(fitimage,0,sizeof(float)*XLEN*YLEN);
    spgridimage_init(xtens,ytens);

    if(ftol > 0.0)
        powell(p,xi,n,ftol,iter,fret,spuppereval);

#if defined(TRAPS)
    sigaction(SIGINT,&ovec,0);
#endif

    f = fopen("uppercoefs.dat","w");
    for(i = 1;i <= n; ++i){
        fprintf(f,"%g\n",undump[i]);
    }
    fflush(f);
    fclose(f);
    spgridimage(fitimage,p);
    spgridimage_done();
    bolomproj(shot,fitimage,fitproj);
    free_vector(p,1,n);
    free_matrix(xi,1,n,1,n);



}

float spuppereval(x)
float x[];
{

    register int i,j,k,l;
    float proj[CHANS];
    float val;
    float fsanity_spupper();
    float wall;
    char *tmp;

    ++funcevals;
		++(*pow_its);
    wall = fsanity_spupper(x);
    if(debug)
        spupperproj(proj,x,((funcevals%debug) == 0));
    else
        spupperproj(proj,x,0);
    val = 0.0;
    for(i = 0; i < nfitchans; ++i){
        val += pow((indata[fitchans[i]] - 
            proj[fitchans[i]]),2.0) / pow(insigma[fitchans[i]],2.0);
    }
    if(nfitchans > (nrfitpts * nzfitpts))
        val = val * wall  / (float)(nfitchans - (nrfitpts * nzfitpts));
    else
        val = val * wall;

    if(debug && (funcevals % debug) == 0) printf("after %d iterations chi = %g\n",funcevals,val);
    return(val);

}
spupperproj(proj,x,f)
float proj[],x[];
int f;
{

    register int i,j,k,l;
    image dimage;
    int ent;
    int cells;
    int xlen,ylen;
    int len;
    int dn;
    float val;
    float gfunc();

    cells = XLEN * YLEN;
    len = gmatrix_len;

    memset(proj,0,CHANS * sizeof(float));
    memset(dimage,0,cells *  sizeof(float));
    spgridimage(dimage,x);
    for(i = 0,j=0,k=1;i < len;++i,j+=2,k+=2){
        proj[PROJ] += gmatrix[i] * dimage[CELL];
    }

    if(f)write_projdat("iter.dat",proj);
    if(f)write_imagedat("iter.sdt",dimage);

}

sanity_spupper(x)
float x[];
{
    int i,j;
    int dn;
    int n;
    n = nrfitpts * nzfitpts;
    for(j=1; j<=  n; ++j){
        if(x[j] < EMIN) x[j] = EMIN;
        if(x[j] > EMAX) x[j] = EMAX;
    }
}
static      float lastx[300],thisx[300];
float fsanity_spupper(x)
float x[];
{
    int i,j;
    int dn;
    int len;
    float wall;
    len = nrfitpts * nzfitpts;
    memcpy(thisx,x,len * sizeof(float));
    sanity_spupper(thisx);
    wall = 1.0;
    for(i=1; i <=len; ++i){
        if(thisx[i] != x[i]) {
            wall += 10.0;
            if(fabs(thisx[i] - x[i]) > fabs(lastx[i] - x[i])) wall += 1.00;
            if(fabs(thisx[i] - x[i]) < fabs(lastx[i] - x[i])) wall -= 1.00;
            lastx[i] = x[i];
        }
    }
    if(wall < 1.0) wall = 1.0;
    /*printf("wall = %g\n",wall);*/
    return(wall);
}
dumpum_spupper()
{
    FILE *f;
    int i;
    int n;
    struct sigaction vec,ovec;
    image gimage;
    n = nrfitpts * nzfitpts;
    printf("dumping unknowns in uppercoefs.dat\n");
    f = fopen("uppercoefs.dat","w");
    for(i = 1;i <= n; ++i){
        fprintf(f,"%g\n",undump[i]);
    }
    fflush(f);
    fclose(f);
    if(prev_handler){
        vec.sa_handler = prev_handler;
	sigfillset(&vec.sa_mask);
        vec.sa_flags = 0;
        sigaction(SIGINT,&vec,0);
        prev_handler();
        longjmp(sjbuf,0) ;
    }  else
        exit(0);
}

