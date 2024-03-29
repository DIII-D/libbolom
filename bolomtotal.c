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
extern int *gpos;
extern int gmatrix_len;
extern int gpos_len;
#include <setjmp.h>
static jmp_buf sjbuf ;
static int funcevals;
extern int *pow_its;
static int ilxp,iuxp;

#define PROJ (gpos[k])
#define CELL (gpos[j])

static void (*prev_handler)();
static float *undump;
static float *indata;
static float *insigma;
static int *fitchans;
static int nfitchans;
static float *rfitpts,*zfitpts,*corefitpts;
static float *urfitpts,*uzfitpts;
static int nrfitpts,nzfitpts,ncorefitpts;
static int nurfitpts,nuzfitpts;
static int numunk;
static int done;
static float xinc,yinc;
float psi_norm[XLEN*YLEN],*inefit;

void tot_spgridimage();




float sptotaleval(x)
float x[];
{

    register int i,j,k,l;
    float proj[CHANS];
    float val;
    float fsanity_sptotal();
    float wall;
    char *tmp;

    ++funcevals;
    ++(*pow_its);
    wall = fsanity_sptotal(x);
    if(debug)
        sptotalproj(proj,x,((funcevals%debug) == 0));
    else
        sptotalproj(proj,x,0);
    val = 0.0;
    for(i = 0; i < nfitchans; ++i){
        val += pow((indata[fitchans[i]] - 
            proj[fitchans[i]]),2.0) / pow(insigma[fitchans[i]],2.0);
    }
    if(nfitchans > numunk)
        val = val * wall  / (float)(nfitchans - numunk);
    else
        val = val * wall;

    if(debug && (funcevals % debug) == 0) printf("after %d iterations chi = %g\n",funcevals,val);
    return(val);

}

void tot_spgridimage_done()
{
    if(ncorefitpts > 0)spcoreimage_done();
    if(nrfitpts > 0 && nzfitpts > 0)splowerimage_done();
    if(nurfitpts > 0 && nuzfitpts > 0)spupperimage_done();
}
sptotalproj(proj,x,f)
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
    tot_spgridimage(dimage,x);
    for(i = 0,j=0,k=1;i < len;++i,j+=2,k+=2){
        proj[PROJ] += gmatrix[i] * dimage[CELL];
    }

    if(f)write_projdat("iter.dat",proj);
    if(f)write_imagedat("iter.sdt",dimage);

}
sanity_sptotal(x)
float x[];
{
    int i,j;
    int dn;
    for(j=1; j<=  numunk; ++j){
        if(x[j] < EMIN) x[j] = EMIN;
        if(x[j] > EMAX) x[j] = EMAX;
    }
}
static      float lastx[300],thisx[300];
float fsanity_sptotal(x)
float x[];
{
    int i,j;
    int dn;
    int len;
    float wall;
    len = numunk;
    memcpy(thisx,x,len * sizeof(float));
    sanity_sptotal(thisx);
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
dumpum_sptotal()
{
    FILE *f;
    int i;
    int n;
    struct sigaction vec,ovec;
    image gimage;
    n = numunk;
    printf("dumping unknowns in totalcoefs.dat\n");
    f = fopen("totalcoefs.dat","w");
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

static int nxpoint,nypoint;
static float *x,*y,**z,*tz,*dy,**xp,*yp,*ltemp;
static float sx,sy,slp1,slpn;
static float xe,xb,ye,yb,xb2,xe2;

void splowerimage_init(xtension,ytension)
float xtension,ytension;
{
    int i;
    nxpoint = nrfitpts;
    nypoint = nzfitpts;
    x = vector(0,nxpoint);
    for(i =0 ; i < nxpoint; ++i)x[i] = rfitpts[i];
    /*for(i =0 ; i < nxpoint; ++i)printf("rfitpts %d %g\n",i,x[i]);*/
    y = vector(0,nypoint);
    for(i =0 ; i < nypoint; ++i)y[i] = zfitpts[i];
    /*for(i =0 ; i < nypoint; ++i)printf("zfitpts %d %g\n",i,y[i]);*/
    sx = xtension;
    sy = ytension;
    slp1 = 0.0;
    slpn = 0.0;
    tz = vector(0,nypoint);
    yp = vector(0,nypoint);
    ltemp = vector(0,nypoint*nxpoint);
    z = matrix(0,nypoint,0,nxpoint);
    xp = matrix(0,nypoint,0,nxpoint);
    for(i = 0, xe = rfitpts[0];i <  nxpoint; ++i)if(xe < rfitpts[i])xe = rfitpts[i];
    for(i = 0, xb = rfitpts[0];i <  nxpoint; ++i)if(xb > rfitpts[i])xb = rfitpts[i];
    for(i = 0, xe2 = xe;i <  nxpoint; ++i)if(xe2 > rfitpts[i] && rfitpts[i] >= 0.0)xe2 = rfitpts[i];
    for(i = 0, xb2 = xb;i <  nxpoint; ++i)if(xb2 < rfitpts[i] && rfitpts[i] <= 0.0)xb2 = rfitpts[i];
    for(i = 0, ye = zfitpts[0];i <  nypoint; ++i)if(ye < zfitpts[i])ye = zfitpts[i];
    for(i = 0, yb = zfitpts[0];i <  nypoint; ++i)if(yb > zfitpts[i])yb = zfitpts[i];


}
void splowerimage(im,iter)
float *im,*iter;
{
    register int i,j,k;
    register int l,m,n;
    float fsp;
    float val;
    float tx,ty,ltx;
    float ymax;
    float curv2();

    int it;
    int num;
    FILE *f;



    num = nypoint * nxpoint;
    for(j = 0; j < nypoint; ++j)
        for(i =0 ; i < nxpoint; ++i){
            z[j][i] = iter[i + 1 + j * nxpoint];
        }


    slp1 = 0.0;
    slpn = 0.0;
    for(k = 0; k < nypoint; ++k){
        curv1(&nxpoint,x,z[k],&slp1,&slpn,xp[k],ltemp,&sx);
    }
    for(i = 0; i < XLEN; ++i){
        for(j = 0,it = 1; j < YLEN; ++j){
            tx = psi_norm[i + j * XLEN];
            if(tx < xb  || tx > xe)continue;
            if(tx > xb2 && tx < xe2) continue;
            ty = YMIN + j * yinc;
            if(ty < yb  || ty > ye)continue;
            for(k = 0,it=1; k < nypoint; ++k){
                tz[k] = curv2(&tx,&nxpoint,x,z[k],xp[k],&sx,&it);
            }
            slp1 = 0.0;
            slpn = 0.0;
            curv1(&nypoint,y,tz,&slp1,&slpn,yp,ltemp,&sy);
            it = 1;
            im[j * XLEN + i] = curv2(&ty,&nypoint,y,tz,yp,&sx,&it);
        }
    }

    num = XLEN * YLEN;
    for(i = 0; i < num; ++i)if(im[i] < 0.0 || bound_flag[i]) im[i] = 0.0;

}
void splowerimage_done()
{
    free_matrix(xp,0,nypoint,0,nxpoint);
    free_matrix(z,0,nypoint,0,nxpoint);
    free_vector(ltemp,0,nypoint*nxpoint);
    free_vector(yp,0,nypoint);
    free_vector(tz,0,nypoint);
    free_vector(y,0,nypoint);
    free_vector(x,0,nxpoint);

}

static float *cx,*cxp,*ctemp;
static float cs;
void spcoreimage_init(tension)
float tension;
{
    int i;
    float inc;

    cx = vector(0,ncorefitpts);
    cxp = vector(0,ncorefitpts);
    cs = tension;
    ctemp = vector(0,ncorefitpts);
    for(i =0 ; i < ncorefitpts; ++i)cx[i] = corefitpts[i];
    /*for(i =0 ; i < ncorefitpts; ++i)printf("core %d %g\n",i,cx[i]);*/

}
void spcoreimage(im,iter)
float *im,*iter;
{
    register int i,j,k;
    register int l,m,n;
    float fsp;
    float val;
    float tx,ty,ltx;
    float ymax;
    float curv2();

    int it;
    int num;
    FILE *f;

    slp1 = 0.0;
    slpn = 0.0;
    curv1(&ncorefitpts,cx,&iter[1],&slp1,&slpn,cxp,ctemp,&cs);
    it = 1;
    for(i = 0; i < XLEN; ++i){
        for(j = 0; j < YLEN; ++j){
            if(j < ilxp || j > iuxp)continue;
            tx = inefit[i + j * XLEN];
            if(tx < corefitpts[0]  || tx > corefitpts[ncorefitpts-1])continue;
            im[j * XLEN + i] = curv2(&tx,&ncorefitpts,cx,&iter[1],cxp,&cs,&it);
            it = 0;
        }
    }

    num = XLEN * YLEN;
    for(i = 0; i < num; ++i)if(im[i] < 0.0 || bound_flag[i]) im[i] = 0.0;

}
void spcoreimage_done()
{
    free_vector(cxp,0,ncorefitpts);
    free_vector(cx,0,ncorefitpts);
    free_vector(ctemp,0,ncorefitpts);

}
static int nuxpoint,nuypoint;
static float *ux,*uy,**uz,*utz,*udy,**xpu,*uyp,*utemp;
static float usx,usy,slp1,slpn;
static float uxe,uxb,uye,uyb,uxb2,uxe2;

void spupperimage_init(xtension,ytension)
float xtension,ytension;
{
    int i;
    nuxpoint = nurfitpts;
    nuypoint = nuzfitpts;
    ux = vector(0,nuxpoint);
    for(i =0 ; i < nuxpoint; ++i)ux[i] = urfitpts[i];
    uy = vector(0,nuypoint);
    for(i =0 ; i < nuypoint; ++i)uy[i] = uzfitpts[i];
    usx = xtension;
    usy = ytension;
    slp1 = 0.0;
    slpn = 0.0;
    utz = vector(0,nuypoint);
    uyp = vector(0,nuypoint);
    utemp = vector(0,nuypoint*nuxpoint);
    uz = matrix(0,nuypoint,0,nuxpoint);
    xpu = matrix(0,nuypoint,0,nuxpoint);
    for(i = 0, uxe = urfitpts[0];i <  nuxpoint; ++i)if(uxe < urfitpts[i])uxe = urfitpts[i];
    for(i = 0, uxb = urfitpts[0];i <  nuxpoint; ++i)if(uxb > urfitpts[i])uxb = urfitpts[i];
    for(i = 0, uxe2 = uxe;i <  nuxpoint; ++i)if(uxe2 > urfitpts[i] && urfitpts[i] >= 0.0)uxe2 = urfitpts[i];
    for(i = 0, uxb2 = uxb;i <  nuxpoint; ++i)if(uxb2 < urfitpts[i] && urfitpts[i] <= 0.0)uxb2 = urfitpts[i];
    for(i = 0, uye = uzfitpts[0];i <  nuypoint; ++i)if(uye < uzfitpts[i])uye = uzfitpts[i];
    for(i = 0, uyb = uzfitpts[0];i <  nuypoint; ++i)if(uyb > uzfitpts[i])uyb = uzfitpts[i];


}
void spupperimage(im,iter)
float *im,*iter;
{
    register int i,j,k;
    register int l,m,n;
    float fsp;
    float val;
    float tx,ty,ltx;
    float ymax;
    float curv2();

    int it;
    int num;
    FILE *f;



    num = nuypoint * nuxpoint;
    for(j = 0; j < nuypoint; ++j)
        for(i =0 ; i < nuxpoint; ++i){
            uz[j][i] = iter[i + 1 + j * nuxpoint];
        }


    slp1 = 0.0;
    slpn = 0.0;
    for(k = 0; k < nuypoint; ++k){
        curv1(&nuxpoint,ux,uz[k],&slp1,&slpn,xpu[k],utemp,&usx);
    }
    for(i = 0; i < XLEN; ++i){
        for(j = 0,it = 1; j < YLEN; ++j){
            tx = psi_norm[i + j * XLEN];
            if(tx < uxb  || tx > uxe)continue;
            if(tx > uxb2 && tx < uxe2) continue;
            ty = YMIN + j * yinc;
            if(ty < uyb  || ty > uye)continue;
            for(k = 0,it=1; k < nuypoint; ++k){
                utz[k] = curv2(&tx,&nuxpoint,ux,uz[k],xpu[k],&usx,&it);
            }
            slp1 = 0.0;
            slpn = 0.0;
            curv1(&nuypoint,uy,utz,&slp1,&slpn,uyp,utemp,&usy);
            it = 1;
            im[j * XLEN + i] = curv2(&ty,&nuypoint,uy,utz,uyp,&usx,&it);
        }
    }

    num = XLEN * YLEN;
    for(i = 0; i < num; ++i)if(im[i] < 0.0 || bound_flag[i]) im[i] = 0.0;

}
void spupperimage_done()
{
    free_matrix(xpu,0,nuypoint,0,nuxpoint);
    free_matrix(uz,0,nuypoint,0,nuxpoint);
    free_vector(utz,0,nuypoint);
    free_vector(uyp,0,nuypoint);
    free_vector(utemp,0,nuypoint*nuxpoint);
    free_vector(uy,0,nuypoint);
    free_vector(ux,0,nuxpoint);

}

void tot_spgridimage_init(ctens,lxtens,lytens,uxtens,uytens)
float ctens,lxtens,lytens,uxtens,uytens;
{

    if(ncorefitpts > 0)spcoreimage_init(ctens);
    if(nrfitpts > 0 && nzfitpts > 0)splowerimage_init(lxtens,lytens);
    if(nurfitpts > 0 && nuzfitpts > 0)spupperimage_init(uxtens,uytens);

}

void tot_spgridimage(im,iter)
float *im,*iter;
{
    image lowerimage,coreimage,upperimage;
    int i;
    memset(lowerimage,0,sizeof(float)*XLEN*YLEN);
    memset(upperimage,0,sizeof(float)*XLEN*YLEN);
    memset(coreimage,0,sizeof(float)*XLEN*YLEN);
    if(ncorefitpts > 0){
        spcoreimage(coreimage,iter);
        for(i = 0; i < XLEN*YLEN; ++i)im[i] += coreimage[i];
    }
    if(nrfitpts > 0 && nzfitpts > 0){
        splowerimage(lowerimage,&iter[ncorefitpts]);
        for(i = 0; i < XLEN*YLEN; ++i)im[i] += lowerimage[i];
    }
    if(nurfitpts > 0 && nuzfitpts > 0){
        spupperimage(upperimage,&iter[ncorefitpts + nrfitpts *nzfitpts]);
        for(i = 0; i < XLEN*YLEN; ++i)im[i] += upperimage[i];
    }
}

bolomfit_global(shot,chans,nchans,inproj,sigma,c,nc,r,nr,z,nz,ur,nur,uz,nuz,ctens,xtens,ytens,uxtens,uytens,efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,ftol,fret,iter)
int shot;
float *inproj,*sigma,*c,*r,*z,*ur,*uz,*fitproj,ftol,*fret;
float ctens,xtens,ytens,uxtens,uytens;
image fitimage;
image efit;
int *chans,nchans,nc,nr,nz,nur,nuz;
int *iter;
float rax,zax,rxpt1,zxpt1,rxpt2,zxpt2;

{

    float sptotaleval();
    float *p,**xi;
    int dumpum_sptotal();
    FILE *f;
    int i,j,k;
    int n;
    int allbound;

    struct sigaction vec,ovec;
    funcevals = 1;

    fitchans = chans;
    nfitchans = nchans;
    indata = inproj;
    insigma = sigma;
    corefitpts = c;
    rfitpts = r;
    zfitpts = z;
    urfitpts = ur;
    uzfitpts = uz;
    ncorefitpts = nc;
    inefit = efit;
    nrfitpts = nr;
    nzfitpts = nz;
    nurfitpts = nur;
    nuzfitpts = nuz;
    numunk = nr * nz +nur * nuz + nc;
    n = numunk + 1;
    done = 0;
    efitgeom(efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,psi_norm);
    write_imagedat("psi_norm.sdt",psi_norm);
    xinc = (XMAX - XMIN) / (float)(XLEN-1);
    yinc = (YMAX - YMIN) / (float)(YLEN-1);
    if(rxpt1 < XMIN || rxpt1 > XMAX ||
        zxpt1 < YMIN || zxpt1 > YMAX){
        ilxp = 0;
    } else ilxp = (zxpt1 - YMIN) / yinc;
    if(rxpt2 < XMIN || rxpt2 > XMAX ||
        zxpt2 < YMIN || zxpt2 > YMAX){
        iuxp = YLEN-1;
    } else iuxp = (zxpt2 - YMIN) / yinc;
    if(!bolom_set_gmatrix(shot))return(0);

    setjmp(sjbuf);
    if(done)return;
    done = 1;
#if defined(TRAPS)
    vec.sa_handler = dumpum_sptotal;
    sigfillset(&vec.sa_mask);
    vec.sa_flags = 0;
    sigaction(SIGINT,&vec,&ovec);
    prev_handler = ovec.sa_handler;
#endif

    write_projdat("totalorig.dat",inproj);

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

    if((f = fopen("totalcoefs.dat","r")) != NULL){
        for(i = 1;i <= n; ++i){
            if(fscanf(f,"%g\n",&undump[i]) != 1)break;
        }
        fclose(f);
    }



    memset(fitimage,0,sizeof(float)*XLEN*YLEN);
    tot_spgridimage_init(ctens,xtens,ytens,uxtens,uytens);

    if(ftol > 0.0)
        powell(p,xi,n,ftol,iter,fret,sptotaleval);

    fprintf(stderr,"iter = %d\n",*iter);
#if defined(TRAPS)
    sigaction(SIGINT,&ovec,0);
#endif

    f = fopen("totalcoefs.dat","w");
    for(i = 1;i <= n; ++i){
        fprintf(f,"%g\n",undump[i]);
    }
    fflush(f);
    fclose(f);
    tot_spgridimage(fitimage,p);
    tot_spgridimage_done();
    bolomproj(shot,fitimage,fitproj);
    free_vector(p,1,n);
    free_matrix(xi,1,n,1,n);



}
