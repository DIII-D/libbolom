
#include "bolom.h"

static int nxpoint,nypoint;
static float *x,*y,**z,*tz,*dy,**xp,*yp,*temp;
static float sx,sy,slp1,slpn;
static float xe,xb,ye,yb,xb2,xe2;
extern int nrfitpts,nzfitpts;

static spgridimage_init(xtension,ytension)
float xtension,ytension;
{
    int i;
    nxpoint = nrfitpts;
    nypoint = nzfitpts;
    x = vector(0,nxpoint);
    for(i =0 ; i < nxpoint; ++i)x[i] = rfitpts[i];
    y = vector(0,nypoint);
    for(i =0 ; i < nypoint; ++i)y[i] = zfitpts[i];
    sx = xtension;
    sy = ytension;
    slp1 = 0.0;
    slpn = 0.0;
    tz = vector(0,nypoint);
    yp = vector(0,nypoint);
    temp = vector(0,nypoint*nxpoint);
    z = matrix(0,nypoint,0,nxpoint);
    xp = matrix(0,nypoint,0,nxpoint);
    for(i = 0, xe = rfitpts[0];i <  nxpoint; ++i)if(xe < rfitpts[i])xe = rfitpts[i];
    for(i = 0, xb = rfitpts[0];i <  nxpoint; ++i)if(xb > rfitpts[i])xb = rfitpts[i];
    for(i = 0, xe2 = xe;i <  nxpoint; ++i)if(xe2 > rfitpts[i] && rfitpts[i] >= 0.0)xe2 = rfitpts[i];
    for(i = 0, xb2 = xb;i <  nxpoint; ++i)if(xb2 < rfitpts[i] && rfitpts[i] <= 0.0)xb2 = rfitpts[i];
    for(i = 0, ye = zfitpts[0];i <  nypoint; ++i)if(ye < zfitpts[i])ye = zfitpts[i];
    for(i = 0, yb = zfitpts[0];i <  nypoint; ++i)if(yb > zfitpts[i])yb = zfitpts[i];


}
static spgridimage(im,iter)
float *im,*iter;
{
    register int i,j,k;
    register int l,m,n;
    float fsp;
    float val;
    float xinc,yinc;
    float tx,ty,ltx;
    float ymax;
    float curv2();

    int it;
    int num;
    FILE *f;


    xinc = 170.0 / 32.0;
    yinc = 320.0 / 64.0;

    num = nypoint * nxpoint;
    for(j = 0; j < nypoint; ++j)
        for(i =0 ; i < nxpoint; ++i){
            z[j][i] = iter[i + 1 + j * nxpoint];
        }


    for(k = 0; k < nypoint; ++k){
        curv1(&nxpoint,x,z[k],&slp1,&slpn,xp[k],temp,&sx);
    }
    for(i = 0; i < XLEN; ++i){
        for(j = 0,it = 1; j < YLEN; ++j){
            tx = psi_norm[i + j * XLEN];
            if(tx < xb  || tx > xe)continue;
            if(tx > xb2 && tx < xe2) continue;
            ty = -160.0 + j * yinc;
            if(ty < yb  || ty > ye)continue;
            for(k = 0,it=1; k < nypoint; ++k){
                tz[k] = curv2(&tx,&nxpoint,x,z[k],xp[k],&sx,&it);
            }
            curv1(&nypoint,y,tz,&slp1,&slpn,yp,temp,&sy);
            it = 1;
            im[j * XLEN + i] = curv2(&ty,&nypoint,y,tz,yp,&sx,&it);
        }
    }

    num = XLEN * YLEN;
    for(i = 0; i < num; ++i)if(im[i] < 0.0 || bound_flag[i]) im[i] = 0.0;

}
static spgridimage_done()
{
    free_matrix(xp,0,nypoint,0,nxpoint);
    free_matrix(z,0,nypoint,0,nxpoint);
    free_vector(tz,0,nypoint);
    free_vector(yp,0,nypoint);
    free_vector(temp,0,nypoint*nxpoint);
    free_vector(y,0,nypoint);
    free_vector(x,0,nxpoint);

}
