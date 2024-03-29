
/* 20110323 whm changed lsei to define, defaults to lsei_  */
/* 20031031 tbt changed lsei to lsei_  */
#if defined(nounderscore)
#define LSEI lsei
#else
#define LSEI lsei_
#endif

#define PREC float
#include <stdio.h>
/*#include <VIEW/user.h>*/
#include <nrutil.h>
#include <math.h>
#define EMIN 0.0
#define EMAX 20.0
static int debug  = 0;
#include "bolom.h"

#include <math.h>
#include <stdlib.h>
#include <nrutil.h>
#include <signal.h>
#include <stdio.h>
#include "boundext.h"
#include <setjmp.h>
static jmp_buf sjbuf ;
static int ilxp;
static float xinc,yinc;
extern float *gmatrix;
extern int *gpos;
extern int gmatrix_len;
extern int gpos_len;
#define PROJ (gpos[k])
#define CELL (gpos[j])
extern int grad[XLEN*YLEN];
extern image psi_norm;

static void (*prev_handler)();
static float *undump;
static float *insigma;
static float *fitpts;
static int nfitpts;
static image inefit;
static int done;
static int iuix,ouix,ibndy;


bolomfit_upper_cells(shot,chans,nchans,inproj,sigma,zinner,zouter,minprivate,mincore,maxinner,maxouter,drwtd,dzwtd,drwtb,dzwtb,efit,rax,zax,
rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,fret)
unsigned int shot;
float *inproj,*sigma,*fitproj,*fret;
float zinner,zouter,minprivate,mincore,maxinner,maxouter;
float drwtd,dzwtd,drwtb,dzwtb;
image efit,fitimage;
int *chans,nchans;
float rax,zax,rxpt1,zxpt1,rxpt2,zxpt2;

{

    float *p,**xi;
    int iter;
    int splcore();
    FILE *f;
    int i,j,k;
    int ws_len,ip_len;
    PREC c1,c2,c3,c4;
    PREC *amat,*cmat,*imat,*timat;
    PREC sum;
    PREC *work;
    int na,nma,eqn,toteqn;
    int l_m,l_n,l_p,l_lda,l_ldb;
    int lwork,info;
    float basisel(),basisv();
    int lda,numunk,numlow;
    int iparm;
    float *g,*h,*w;
    int mdw,me,ma,mg;
    int *ip;
    float *ws;
    float prgopt[100],rnorme,rnorml;
    int mode;
    int allbound;
    int here;


    PREC *u;
    PREC *tz,*uz;
    PREC *sol;

    struct sigaction vec,ovec;
    efitgeom(efit,rax,zax,rxpt1,zxpt1,rxpt2,zxpt2,psi_norm);
    if(!bolom_set_gmatrix(shot))return(0);
    write_imagedat("psi_norm.sdt",psi_norm);


    done = 0;
    setjmp(sjbuf);
    if(done)return;
    done = 1;

#if defined(TRAPS)
    vec.sa_handler = splcore;
    sigfillset(&vec.sa_mask);
    vec.sa_flags = 0;
    sigaction(SIGINT,&vec,&ovec);
    prev_handler = ovec.sa_handler;
#endif

    /*write_projdat("coreorig.dat",inproj);*/

    xinc = (XMAX - XMIN) / (float)(XLEN-1);
    yinc = (YMAX - YMIN) / (float)(YLEN-1);
    drwtd *= yinc / xinc;
    drwtb *= yinc / xinc;
    printf("rw=%g zw=%g\n",drwtd,dzwtd);

    if(rxpt2 < XMIN || rxpt2 > XMAX ||
        zxpt2 < YMIN || zxpt2 > YMAX){
        ilxp = 0;
    } else ilxp = (zxpt2 - YMIN) / yinc;
    ilxp *= XLEN;

    if(zinner < YMIN || zinner > YMAX){
        iuix = 0;
    } else iuix = (zinner - YMIN) / yinc;
    iuix *= XLEN;

    if(zouter < YMIN || zouter > YMAX){
        ouix = 0;
    } else ouix = (zouter - YMIN) / yinc;
    ouix *= XLEN;

    memset(grad,0L,sizeof(grad));
    here = 0;


    printf("ilxp = %d  iuix = %d  ouix = %d\n",ilxp,iuix,ouix);

#if defined(TRAPS)
    sigaction(SIGINT,&ovec,0);
#endif




    na = XLEN * YLEN;
    /*printf("na:%d nchans:%d\n",na,nchans);*/
    numlow = 0;

/*  Do a flat field back projection so we exclude cells with no measured response */
    for(i = 0; i < CHANS; ++i) fitproj[i] = 1;
    bolombackproj(shot,fitimage,fitproj);

    ibndy=(iuix<ouix)?iuix:ouix;
    for(i = ibndy+1;i <na;++i){
        if(bound_flag[i] || fitimage[i] == 0.0)continue;
        if(psi_norm[i] <= 0 && efit[i] >= minprivate &&
            efit[i] <= maxinner && i > ilxp){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
        if(psi_norm[i] > 0 && efit[i] >= minprivate &&
            efit[i] <= maxouter && i > ilxp){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
        if(psi_norm[i] <= 0 && efit[i] >= mincore &&
            efit[i] <= maxinner && i > iuix && i <= ilxp){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
        if(psi_norm[i] > 0 && efit[i] >= mincore &&
            efit[i] <= maxouter && i > ouix && i <= ilxp){
            grad[i] = numlow + 1.0;
            ++numlow;
        }
    }

    lda = numlow;
    toteqn = 2*lda + nchans;
    printf("numlow = %d lda = %d toteqn = %d\n",numlow,lda,toteqn);
    cmat = calloc(sizeof(PREC),na);

    imat = calloc(sizeof(PREC) * 2,lda * toteqn);
    timat = calloc(sizeof(PREC) * 2,lda * toteqn);

    amat = calloc(sizeof(PREC) * 2,na * (CHANS + 1));
    for(i = 0,j=0,k=1;i < gmatrix_len;++i,j+=2,k+=2){
        if(CELL < na && PROJ < CHANS ){
            amat[PROJ * na + CELL] = gmatrix[i];
        }
    }
    for(i=0;i<nchans;++i){
        for(j = 0,k = 0;k < na;++k){
            if(grad[k]){
                imat[i * lda + j] = amat[chans[i]*na+k]/sigma[chans[i]];
                ++j;
            }
        }
    }
    /*radial smoothing*/
    for(j=0,k=ibndy;k < ilxp;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans)*lda+j] += -2 * drwtb;
            if((k%XLEN)>0 && grad[k-1] > 0.0){
                imat[(j+nchans)*lda+(j-1)] += 1 * drwtb;
            }
            if((k%XLEN)<(XLEN-1) && grad[k+1] > 0.0){
                imat[(j+nchans)*lda+(j+1)] += 1 * drwtb;
            }
            ++j;
        }
    }
    for(k=ilxp;k < na;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans)*lda+j] += -2 * drwtd;
            if((k%XLEN)>0 && grad[k-1] > 0.0){
                imat[(j+nchans)*lda+(j-1)] += 1 * drwtd;
            }
            if((k%XLEN)<(XLEN-1) && grad[k+1] > 0.0){
                imat[(j+nchans)*lda+(j+1)] += 1 * drwtd;
            }
            ++j;
        }
    }

    /*vertical smoothing*/
    for(j=0,k=ibndy;k < ilxp;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans+numlow)*lda+j] += -2 * dzwtb;
            if(k>XLEN && grad[k-XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k-XLEN]-1] += 1 * dzwtb;
            }
            if(k<(na-XLEN) && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k+XLEN]-1] += 1 * dzwtb;
            }
            ++j;
        }
    }
    for(k=ilxp;k < na;++k){
        if(grad[k] > 0.0){
            imat[(j+nchans+numlow)*lda+j] += -2 * dzwtd;
            if(k>XLEN && grad[k-XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k-XLEN]-1] += 1 * dzwtd;
            }
            if(k<(na-XLEN) && grad[k+XLEN] > 0.0) {
                imat[(j+nchans+numlow)*lda+grad[k+XLEN]-1] += 1 * dzwtd;
            }
            ++j;
        }
    }


    /*    for(i = nchans;i < toteqn; ++i){
        for(k=0;k < lda;++k){
            printf("%g ",imat[i*lda+k]);
        }
        printf("\n");
    }*/


    eqn = 0;
    /*u = calloc(sizeof(PREC)*2,eqn * lda);*/

    /*uz = calloc(sizeof(PREC)*2,eqn);*/
    nma = 0;


    /*printf("constraint matrix is %d x %d\n",numunk,nma);*/

    sol = calloc(sizeof(PREC) * 2,lda);
    tz = calloc(sizeof(PREC) * 2,toteqn);
    for(i = 0; i < nchans; ++i)tz[i] = inproj[chans[i]]/sigma[chans[i]];

    for(i = 0; i < toteqn ; ++i){
        for(j = 0; j < lda; ++j){
            timat[j * toteqn + i] = imat[i *lda + j];
        }
    }

    lwork = 10*(lda + nma + toteqn);
    work = calloc(sizeof(PREC) * 2,lwork);
    l_m = toteqn;
    l_n = lda;
    l_p = nma;
    l_lda = toteqn;
    l_ldb = eqn;

    mg = lda;
    me = nma;
    ma = toteqn;
    g = calloc(sizeof(float),mg * mg);
    h = calloc(sizeof(float),mg);
    for(j = 0; j < mg; ++j){
        g[j * mg + j] = 1.0;
    }

    mdw = me + ma + mg;
    w = calloc(sizeof(float),mdw * (lda+1));

    if(me > 0){
        for(i = 0; i < me; ++i){
            for(j = 0; j < lda; ++j){
                w[i + j * mdw] = u[i + j * eqn];
            }
        }
        for(i = 0; i < me; ++i){
            w[i + lda * mdw] = uz[i];
        }
    }

    if(ma > 0){
        for(i = 0; i < toteqn; ++i){
            for(j = 0; j < lda; ++j){
                w[i + me + j * mdw] = imat[j + i * lda];
            }
        }
        for(i = 0;i < toteqn; ++i){
            w[i + me + lda * mdw] = tz[i];
        }
    }

    if(mg > 0){
        for(i = 0;i < mg; ++i){
            for(j = 0; j < lda; ++j){
                w[i + me + ma + j * mdw] = g[i * lda + j];
            }
        }
        for(i = 0;i < mg; ++i){
            w[i + me + ma + lda * mdw] = h[i];
        }
    }




    k = ((ma+mg) > lda)?ma+mg:lda;
    ip_len = mg+2*lda+2;
    ws_len = 2*(me+lda)+(mg+2)*(lda+7)+k;
    ip = calloc(sizeof(int),ip_len);
    ws = calloc(sizeof(float),ws_len);
    ip[0] = ws_len;
    ip[1] = ip_len;
    prgopt[0] = 4;
    prgopt[1] = 4;
    prgopt[2] = 1.0e-4;
    prgopt[3] = 7;
    prgopt[4] = 5;
    prgopt[5] = 1.0e-4;
    prgopt[6] = 9;
    prgopt[7] = 1;
    prgopt[8] = 1;
    LSEI(w,&mdw,&me,&ma,&mg,&lda,prgopt,sol,&rnorme,&rnorml,&mode,ws,ip);
    printf("lsei mode = %d\n",mode);
    memset(fitimage,0,XLEN*YLEN*sizeof(float));
    for(j = 0,k = 0;k < na;++k){
        if(grad[k]){
            fitimage[k] = sol[j];
            ++j;
        }
    }


    /*for(i = 0 ; i < lda; ++i)printf("sol %d %g\n",i,sol[i]);*/

    /*write_imagedat("core.sdt",fitimage);*/

    bolomclip(shot,fitimage,chans,nchans);
    bolomproj(shot,fitimage,fitproj);
    *fret = 0.0;
    for(i = 0; i < nchans; ++i){
        *fret += pow((inproj[chans[i]] - 
            fitproj[chans[i]]),2.0) / pow(sigma[chans[i]],2.0);
    }
    if(nchans > lda)
        *fret = *fret / (float)(nchans - lda);
    /*fprintf(stderr,"core chi = %g\n",*fret);*/
    free(cmat);
    free(imat);
    free(timat);
    free(amat);
    free(sol);
    free(tz);
    free(work);
    free(g);
    free(h);
    free(w);
    free(ws);

}
/*matherr(x)
register struct exception *x;
{
    if(strcmp(x->name,"sqrt") == 0){
        if(x->arg1 < 0.0)x->retval = sqrt(-x->arg1);
        else x->retval = 0;
        return(1);
    } else if(strcmp(x->name,"log") == 0){
        if(x->arg1 < 0.0)x->retval = -MAXFLOAT;
        else if(x->arg1 == 0.0)x->retval = -MAXFLOAT;
        else x->retval = -MAXFLOAT;
        return(1);
    } else if(strcmp(x->name,"pow") == 0){
        if(x->arg1 < 0.0 && (floor(x->arg2) - x->arg2) != 0.0 )
            x->retval = pow(x->arg1,floor(x->arg2));
        else if(x->arg1 == 0.0)x->retval = 0;
        else x->retval = -MAXFLOAT;
        return(1);
    }
} */
