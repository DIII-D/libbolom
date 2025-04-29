/*                                                                                */
/*    whm 20110323 changed names to define, defaults to end with underscore  */
/*    tbt 20031031 Added "_" to end of basisv  for gamini.gat.com compatibility.  */
/*    tbt 20031031 Added "_" to end of basisel for gamini.gat.com compatibility.  */
/*    tbt 20031031 Added "_" to end of lsei    for gamini.gat.com compatibility.  */
/*    tbt 20031031 Added "_" to end of bacnst  for gamini.gat.com compatibility.  */

#if defined(nounderscore)
#define BASISV basisv
#define BASISEL basisel
#define LSEI lsei
#define BACNST bacnst
#else
#define BASISV basisv_
#define BASISEL basisel_
#define LSEI lsei_
#define BACNST bacnst_
#endif

#define PREC float
#include <stdio.h>
/*#include <VIEW/user.h>*/
/*#include <nrutil.h>*/
#include <math.h>
#define EMIN 0.0
#define EMAX 20.0
int debug  = 0;
#include "bolom.h"
#define POSITIVE_EMISSION

#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include "boundext.h"
#include <setjmp.h>
static jmp_buf sjbuf ;
static int ilxp,iuxp;
static float xinc,yinc;
extern float *gmatrix;
extern int *gpos;
extern int gmatrix_len;
extern int gpos_len;
#define PROJ (gpos[k])
#define CELL (gpos[j])

static void (*prev_handler)();
static float *undump;
static float *insigma;
static float *fitpts;
static int nfitpts;
static image inefit;
static int done;



bolomfit_core(shot,chans,nchans,inproj,sigma,kpsi,nkpsi,tens,efit,rax,zax,
rxpt1,zxpt1,rxpt2,zxpt2,fitproj,fitimage,errorimage,fret)
unsigned int shot;
float *inproj,*sigma,*kpsi,*fitproj,*fret,tens;
image efit,fitimage,errorimage;
int *chans,nchans,nkpsi;
float rax,zax,rxpt1,zxpt1,rxpt2,zxpt2;

{

    float *p,**xi;
    int iter;
    int splcore();
    FILE *f;
    int i,j,k;
    PREC T;
    PREC c1,c2,c3,c4;
    PREC *amat,*cmat,*imat,*timat;
    PREC sum;
    PREC *work;
    int na,nma,eqn;
    int l_m,l_n,l_p,l_lda,l_ldb;
    int lwork,info;
    float tproj[71];
    float BASISEL(),BASISV();
    int numunk;
    int iparm;
    float *g,*h,*w;
    int mdw,me,ma,mg;
    int ip[4096];
    float *ws;
    float prgopt[100],rnorme,rnorml;
    int mode;
    int allbound;


    PREC *u;
    PREC *tz,*uz;
    PREC *sol;

#if defined(TRAPS)
    struct sigaction vec,ovec;
		if(!bolom_set_gmatrix(shot))return(0);


    done = 0;
    setjmp(sjbuf);
    if(done)return;
    done = 1;

    vec.sa_handler = splcore;
    sigfillset(&vec.sa_mask);
    vec.sa_flags = 0;
    sigaction(SIGINT,&vec,&ovec);
    prev_handler = ovec.sa_handler;
#endif

    /*write_projdat("coreorig.dat",inproj);*/

    for(i=0;i<71;++i)tproj[i] = 0.0;
    for(i=0;i<nchans;++i){
        tproj[chans[i]] = inproj[chans[i]];
    }
    /*for(i=0;i<71;++i)inproj[i] = tproj[i];*/
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
    ilxp *= XLEN;
    iuxp *= XLEN;


#if defined(TRAPS)
    sigaction(SIGINT,&ovec,0);
#endif







    //printf("number of knots= %d\n",nkpsi);

    na = XLEN * YLEN;
    //printf("na:%d nchans:%d\n",na,nchans);
    numunk = 2 * nkpsi;
    cmat = calloc(sizeof(PREC),numunk * na);

    for(k = 0; k < na; ++k){
        for(j = 0; j < numunk; ++j){
            iparm = j + 1;
            cmat[k * numunk + j] = BASISEL(&iparm,(float *)&(efit[k]),
                &T,kpsi,&nkpsi);
        }
    }

    imat = calloc(sizeof(PREC) * 2,numunk * nchans);
    timat = calloc(sizeof(PREC) * 2,numunk * nchans);

    amat = calloc(sizeof(PREC) * 2,na * (CHANS + 1));
    for(i = 0,j=0,k=1;i < gmatrix_len;++i,j+=2,k+=2){
        if(CELL < na && PROJ < CHANS ){
            if(bound_flag[CELL])continue;
            if(CELL < ilxp || CELL > iuxp)continue;
            amat[PROJ * na + CELL] = gmatrix[i];
        }
    }

    for(i=0;i<nchans;++i){
        for(j=0;j<numunk;++j){
            for(k=0,sum=0.0;k<na;++k) sum += amat[chans[i]*na+k] * cmat[k*numunk+j];
            imat[i * numunk + j] = sum;  /*  passed for a fortran rout, row major  */
        }
    }


    eqn =  3*nkpsi + 6;
    u = calloc(sizeof(PREC)*2,eqn * numunk);

    uz = calloc(sizeof(PREC)*2,eqn);
    nma = 0;
    BACNST(&numunk,&T,kpsi,&nkpsi,u,&eqn,&nma,uz);


    //printf("least squares matrix is %d x %d\n",numunk,nchans);
    //printf("constraint matrix is %d x %d\n",numunk,nma);

    sol = calloc(sizeof(PREC) * 2,numunk);
    tz = calloc(sizeof(PREC) * 2,nchans);
    for(i = 0; i < nchans; ++i)tz[i] = inproj[chans[i]];

    for(i = 0; i < nchans ; ++i){
        for(j = 0; j < numunk; ++j){
            timat[j * nchans + i] = imat[i *numunk + j];
        }
    }

    lwork = 10*(numunk + nma + nchans);
    work = calloc(sizeof(PREC) * 2,lwork);
    l_m = nchans;
    l_n = numunk;
    l_p = nma;
    l_lda = nchans;
    l_ldb = eqn;
#ifndef  POSITIVE_EMISSION
    sgglse(&l_m,&l_n,&l_p,timat,&l_lda,u,&l_ldb,tz,uz,sol,work,&lwork,&info);
    if(info != 0){
        /*printf("info = %d\n",info);*/
        return;
    }
    //printf("lwork = %d work(1) = %g\n",lwork,work[0]);
#else

    mg = numunk / 2;
    me = nma;
    ma = nchans;
    g = calloc(sizeof(float),mg * numunk);
    h = calloc(sizeof(float),mg);
    for(j = 0; j < mg; ++j){
        g[j * numunk + 2*j] = 1.0;
    }

    mdw = me + ma + mg;
    w = calloc(sizeof(float),mdw * (numunk+1));

    if(me > 0){
        for(i = 0; i < me; ++i){
            for(j = 0; j < numunk; ++j){
                w[i + j * mdw] = u[i + j * eqn];
            }
        }
        for(i = 0; i < me; ++i){
            w[i + numunk * mdw] = uz[i];
        }
    }

    if(ma > 0){
        for(i = 0; i < ma; ++i){
            for(j = 0; j < numunk; ++j){
                w[i + me + j * mdw] = imat[j + i * numunk];
            }
        }
        for(i = 0;i < ma; ++i){
            w[i + me + numunk * mdw] = tz[i];
        }
    }

    if(mg > 0){
        for(i = 0;i < mg; ++i){
            for(j = 0; j < numunk; ++j){
                w[i + me + ma + j * mdw] = g[i * numunk + j];
            }
        }
        for(i = 0;i < mg; ++i){
            w[i + me + ma + numunk * mdw] = h[i];
        }
    }




    ws = calloc(sizeof(float),128000);
    ip[0] = 128000;
    ip[1] = 4096;
    memset(ip,0,sizeof(ip));
    prgopt[0] = 4;
    prgopt[1] = 4;
    prgopt[2] = 1.0e-4;
    prgopt[3] = 7;
    prgopt[4] = 5;
    prgopt[5] = 1.0e-4;
    prgopt[6] = 9;
    prgopt[7] = 1;
    prgopt[8] = 1;
    LSEI(w,&mdw,&me,&ma,&mg,&numunk,prgopt,sol,&rnorme,&rnorml,&mode,ws,ip);








#endif



    memset(fitimage,0,XLEN*YLEN*sizeof(float));
    for(k = 0; k < na; ++k){
        if(bound_flag[k])continue;
        if(k < ilxp || k > iuxp)continue;
        if(efit[k] >= kpsi[0] && efit[k] <= kpsi[nkpsi-1])
            fitimage[k] = BASISV(&numunk,&efit[k],&T,sol,
                kpsi,&nkpsi);
        else fitimage[k] = 0.0;
        /*if(fitimage[k] < 0.0)fitimage[k] = 0.0;*/
    }
    memset(errorimage,0,XLEN*YLEN*sizeof(float));
#ifdef  POSITIVE_EMISSION
    for(i = 0; i < numunk; ++i){
				if(w[i+i*mdw] >= 0.0)sol[i] = sqrt(w[i + i*mdw]);
				else sol[i] = 0.0;
    }
		T=10.0;
    for(k = 0; k < na; ++k){
        if(bound_flag[k])continue;
        if(k < ilxp || k > iuxp)continue;
        if(efit[k] >= kpsi[0] && efit[k] <= kpsi[nkpsi-1])
            errorimage[k] = BASISV(&numunk,&efit[k],&T,sol,
                kpsi,&nkpsi);
        else errorimage[k] = 0.0;
        /*if(errorimage[k] < 0.0)errorimage[k] = 0.0;*/
    }
#endif
    /*write_imagedat("core.sdt",fitimage);*/

    bolomclip(shot,fitimage,chans,nchans);
    bolomproj(shot,fitimage,fitproj);
    *fret = 0.0;
    for(i = 0; i < nchans; ++i){
        *fret += pow((inproj[chans[i]] - 
            fitproj[chans[i]]),2.0) / pow(sigma[chans[i]],2.0);
    }
    if(nchans > numunk)
        *fret = *fret / (float)(nchans - numunk);
    //fprintf(stderr,"core chi = %g\n",*fret);
    //printf("core chi = %g\n",*fret);
    free(cmat);
    free(imat);
    free(timat);
    free(amat);
    free(u);
    free(uz);
    free(sol);
    free(tz);
    free(work);
    free(g);
    free(h);
    free(w);
    free(ws);

}
splcore()
{
    FILE *f;
    int i;
    image gimage;
    struct sigaction vec,ovec;
    printf("dumping unknowns in corecoefs.dat\n");
    f = fopen("corecoefs.dat","w");
    for(i = 1;i <= nfitpts; ++i){
        fprintf(f,"%g\n",undump[i]);
    }
    fflush(f);
    close(f);
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
    } else if(strcmp(x->name,"exp") == 0){
        if(x->arg1 < 0.0)
            x->retval = 0.0;
        else x->retval = MAXFLOAT;
        return(1);
    } else {
        x->retval = MAXFLOAT;
    }
    return(0);
}*/
c()
{

}
