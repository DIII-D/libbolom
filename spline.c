/*      spline - interpolate points in a tabulated function

        Usage...
                spline [file][options]
                options are:
                  -a  [step [start]] automatic abscissas 
                  -b                 break interpolation after each label 
                  -c                 general curve 
                  -n  num            interpolate over num intervals 
                  -q                 increase intervals fourfold 
                  -t  num            tension in interpolating curve 
                  -x  [min [max]]    interpolate from min to max 
                  -xl                take logs of x values before interpolating 
                  -yl                take logs of y values before interpolating 

        Notes...
                This program fits a smooth curve through a given set of points,
                using the splines under tension introduced by Schweikert [J. 
                Math.  and Physics v 45 pp 312-317 (1966)].  They are similar
                to cubic splines, but are less likely to introduce spurious
                inflection points.

        History...
                3 Jun 86        -s switch allows specified slopes at ends.
                2 Jun 86        Ver 1.3: Fixed several embarrassing bugs in -b option.
                                        -q switch quadruples number of intervals.
                6 Apr 86        -b switch breaks interpolation at labels.
                21 Sep 85       Added -xl and -yl switches for taking logs
                23 Nov 85       Passing lines starting with ';' unchanged, otherwise
                                        ignoring them.
                22 Jan 90       added missing * in FILE *file;
                                (suggested by Jean-Marie de Montarby)

        Author...
                Copyright (c) 1985, 1986  James R. Van Zandt (jrv@mbunix.mitre.org)

                All rights reserved.
                This program may be copied for personal, non-profit use only.

                Based on algorithms by A.  K.  Cline [Communications of the ACM
                v 17 n 4 pp 218-223 (Apr 74)].

*/

#include <stdio.h>
#include <math.h>

#define VERSION "1.3"
char *RIGHT=            "Copyright (c) 1985, 1986  James R. Van Zandt";

#define ENTRIES 200
#define MAXLABELS 50

double *x, *y, *temp, *xp, *yp, *path;
double slp1=0.,slpn=0.,sigma=1.;
double abscissa=0.;
double abscissa_step=1.;
double slope1=0.,slopen=0.;     /* slopes supplied by user */
double xmin=0.;
double xmax=0.;
double curv2(), mylog();

int xlog=0;             /* nonzero if taking log of x values */
int ylog=0;             /* nonzero if taking log of y values */
int debugging=0;
int breaking=0;         /* nonzero if breaking at labels */
int labels=0;           /* number of labels in input data */
int quadruple=0;        /* nonzero if user asking for 4 times the resolution */
int automatic_abscissas=0;
int slopes=0;           /* nonzero if slopes supplied by user */
int abscissa_arguments=0;
int curve=0;            /* nonzero if general curve permitted */
int x_arguments=0;
int n;                          /* number of entries in x, y, yp, and temp */
int index_array[MAXLABELS];     /* indexes into x and y */
int *p_data=index_array;
int total=100;

FILE *file;

char *label_array[MAXLABELS];
char **p_text=label_array;

main(argc,argv) int argc; 
char **argv;
{       
    int nn,origin;
    read_data(argc,argv);
    if(breaking && labels)
    {
        origin=0;
        while(labels--)
        {
            p_data[0] -= origin;
            nn=p_data[0]+1;
            if(quadruple) total=nn*4-3;
            if(nn) curv0(nn,total);
            origin += nn;
            path += nn;
            x += nn;
            y += nn;
            p_data++;
            p_text++;
        }
    }
    else
    {
        if(quadruple) total=n*4-3;
        curv0(n,total);
    }
}

curv0(n,requested)
int n,requested;
{       
    int i, j, each, stop, seg=0, regressing=0;
    double s,ds,xx,yy;
    for(i=1; i<n; i++) if(path[i]<=path[i-1]) regressing++;
    if (regressing)
    {
        fprintf(stderr,"ERROR: independent variable not strictly increasing\n");
        return;
    }
    if (curve && (xmin<0. || xmax>1.))
    {
        fprintf(stderr,"ERROR: independent variable not in range 0 to 1\n");
        return;
    }
    if(curve) {
        curv1(path,x,xp,n); 
        curv1(path,y,yp,n);
    }
    else {
        slp1=slope1; 
        slpn=slopen; 
        curv1(x,y,yp,n);
    }
    s=path[0];
    seg=0;
    if(requested<n) requested=n;
    if(x_arguments==2)                      /* specific range was requested */
    {
        ds=(xmax-xmin)/requested;
        if(curve)
        {
            ds*=(path[n-1]-path[0]);
            s=xmin*(path[n-1]-path[0])+path[0];
        }
        else s=xmin;
        for(i=requested+1; i; i--)
        {
            value(s,&xx,&yy,n);
            printf("%g %g ",xx,yy);
            if(j==p_data[seg]) puts(p_text[seg++]);
            putchar('\n');
            s+=ds;
        }
    }
    else /* spline for entire curve was requested */
    {
        for (i=j=0; j<n-1; j++)
        {
            stop=requested*(j+1)/(n-1);
            ds=(path[j+1]-path[j])/(stop-i);
            for ( ; i<stop; i++)
            {
                value(s,&xx,&yy,n);
                printf("%g %g ",xx,yy);
                if(j==p_data[seg]) puts(p_text[seg++]);
                putchar('\n');
                s+=ds;
            }
        }
        xx=x[n-1]; 
        if(xlog) xx=exp(xx);
        yy=y[n-1]; 
        if(ylog) yy=exp(yy);
        printf("%g %g ",xx,yy);
        if(j==p_data[seg]) puts(p_text[seg++]);
        putchar('\n');
    }
}

value(s,q,r,n) double s,*q,*r; 
int n;
{       
    if(curve)
    {
        *q=curv2(path,x,xp,s,n);
        *r=curv2(path,y,yp,s,n);
    }
    else
    {
        *q=s;
        *r=curv2(x,y,yp,s,n);
    }
    if(xlog) *q=exp(*q);
    if(ylog) *r=exp(*r);
}

read_data(argc,argv) int argc; 
char **argv;
{       
    int i,j,length;
    double xx,yy,d,*pd,sum;
    char *s,*t;
#define BUFSIZE 200
    static char buf[BUFSIZE];

    x=path=malloc(ENTRIES*sizeof(double));
    y=malloc(ENTRIES*sizeof(double));
    if(x==0 || y==0) {
        fprintf(stderr,"can\'t allocate buffer"); 
        exit();
    }
    if(argc>1 && streq(argv[1],"?")) help();
    if(argc<=1 || *argv[1]=='-') file=stdin;
    else
    {
        if(argc>1)
        {
            file=fopen(argv[1],"r");
            if(file==0) {
                printf("file %s not found\n",argv[1]); 
                exit();
            }
            argc--; 
            argv++;
        }
        else help();
    }
    argc--; 
    argv++;
    while(argc>0)
    {
        i=get_parameter(argc,argv);
        argv=argv+i; 
        argc=argc-i;
    }
    if(sigma<=0.)
    {
        fprintf(stderr,"ERROR: tension must be positive\n");
        exit(1);
    }
    if(slopes && curve)
    {
        fprintf(stderr,"ERROR: slopes can\'t be specified for general curve\n");
        exit(1);
    }
    sigma *= -1.;
    if(xlog && !curve)
    {
        if(x_arguments>1) xmax=mylog(xmax);
        if(x_arguments>=1) xmin=mylog(xmin);
    }
    if(automatic_abscissas && abscissa_arguments<2 && x_arguments>=1)
        abscissa=xmin;
    p_data[0]=-1;
    i=0;
    while(i<ENTRIES)
    {
        if(fgets(buf,BUFSIZE,file)==0)break;
        t=buf;
        while(*t && isspace(*t)) t++;
        if(*t == 0) continue;           /* skip blank lines */
        buf[strlen(buf)-1]=0;           /* zero out the line feed */
        if(buf[0]==';') {
            printf("%s\n",buf); 
            continue;
        }  /* skip comment */
        if(automatic_abscissas)
        {
            x[i]=abscissa;
            abscissa+=abscissa_step;
            sscanf(buf,"%F",&y[i]);
            if(ylog) y[i]=mylog(y[i]);
        }
        else
        {
            sscanf(buf,"%F %F",&x[i],&y[i]);
            if(xlog) x[i]=mylog(x[i]);
            if(ylog) y[i]=mylog(y[i]);
        }
        s=buf;                      /* start looking for label */
        while(*s==' ')s++;                      /* skip first number */
        while(*s && (*s!=' '))s++;
        if(!automatic_abscissas)        /* skip second number */
        {
            while(*s==' ')s++;
            while(*s && (*s!=' '))s++;
        }
        while(*s==' ')s++;
        if((length=strlen(s))&&(labels<MAXLABELS))
        {
            if(*s=='\"')           /* label in quotes */
            {
                t=s+1;
                while(*t && (*t!='\"')) t++;
                t++;
            }
            else /* label without quotes */
            {
                t=s;
                while(*t && (*t!=' '))t++;
            }
            *t=0;
            length=t-s;
            p_data[labels]=i;
            p_text[labels]=malloc(length+1);
            if(p_text[labels]) strcpy(p_text[labels++],s);
        }
        i++;
    }
    n=i;
    if(breaking && (!labels || p_data[labels-1]!=n-1))
    {
        p_data[labels]=i-1;
        if(p_text[labels]=malloc(1)) *p_text[labels++]=0;
    }
    yp=malloc(n*sizeof(double));
    temp=malloc(n*sizeof(double));
    if(temp==0 || yp==0) {
        fprintf(stderr,"can\'t allocate buffer"); 
        exit();
    }
    if(curve)
    {
        xp=malloc(n*sizeof(double));
        path=malloc(n*sizeof(double));
        if(xp==0|| path==0) {
            fprintf(stderr,"can\'t allocate buffer"); 
            exit();
        }
        path[0]=sum=0.;
        for (i=1; i<n; i++)
        {
            xx=x[i]-x[i-1];
            yy=y[i]-y[i-1];
            path[i]=(sum+=sqrt(xx*xx + yy*yy));
        }
        /*              for(i=0; i<n; i++)
                        printf("path[%d]=%g  x[%d]=%g \n",i,path[i],i,x[i]); */
    }
}


/* get_parameter - process one command line option
                (return # parameters used) */
get_parameter(argc,argv) int argc; 
char **argv;
{       
    int i;
    if(streq(*argv,"-a"))
    {
        i=get_double(argc,argv,2,&abscissa_step,&abscissa,&abscissa);
        abscissa_arguments=i-1;
        automatic_abscissas=1;
        return i;
    }
    else if(streq(*argv,"-b")) {
        breaking=1; 
        return 1;
    }
    else if(streq(*argv,"-c")) {
        curve=1; 
        return 1;
    }
    else if(streq(*argv,"-d")) {
        debugging=1; 
        return 1;
    }
    else if(streq(*argv,"-n"))
    {
        if((argc>1) && numeric(argv[1])) total=atoi(argv[1]);
        return 2;
    }
    else if(streq(*argv,"-q")) {
        quadruple=1; 
        return 1;
    }
    else if(streq(*argv,"-t"))
    {
        return(get_double(argc,argv,1,&sigma,&abscissa,&abscissa));
    }
    else if(streq(*argv,"-s"))
    {
        slopes++;
        return(get_double(argc,argv,2,&slope1,&slopen,&slopen));
    }
    else if(streq(*argv,"-x"))
    {
        i=get_double(argc,argv,2,&xmin,&xmax,&xmax);
        x_arguments=i-1;
        return i;
    }
    else if(streq(*argv,"-xl")) {
        xlog++; 
        return 1;
    }
    else if(streq(*argv,"-yl")) {
        ylog++; 
        return 1;
    }
    else gripe(argv);
}

get_double(argc,argv,permitted,a,b,c)
int argc,permitted; 
char **argv; 
double *a,*b,*c;
{       
    int i=1;
    if((permitted--)>0 && (argc>i) && numeric(argv[i])) *a=atof(argv[i++]);
    if((permitted--)>0 && (argc>i) && numeric(argv[i])) *b=atof(argv[i++]);
    if((permitted--)>0 && (argc>i) && numeric(argv[i])) *c=atof(argv[i++]);
    return i;
}

int streq(a,b) char *a,*b;
{       
    while(*a)
    {
        if(*a!=*b)return 0;
        a++; 
        b++;
    }
    return 1;
}

gripe_arg(s) char *s;
{       
    fprintf(stderr,"argument missing for switch %s",s);
    help();
}

gripe(argv) char **argv;
{       
    fprintf(stderr,*argv); 
    fprintf(stderr," isn\'t a legal argument \n\n");
    help();
}

help()
{       
    fprintf(stderr,"spline   version %s",VERSION);
    fprintf(stderr,"\nusage: spline [file][options]\n");
    fprintf(stderr,"options are:\n");
    fprintf(stderr,"  -a  [step [start]] automatic abscissas \n");
    fprintf(stderr,"  -b                 break interpolation after each label \n");
    fprintf(stderr,"  -c                 general curve \n");
    fprintf(stderr,"  -n  num            interpolate over num intervals \n");
    fprintf(stderr,"  -q                 increase intervals fourfold \n");
    fprintf(stderr,"  -s  [num [num]]    specify slopes at beginning and end \n");
    fprintf(stderr,"  -t  num            tension in interpolating curve \n");
    fprintf(stderr,"  -x  [min [max]]    interpolate from min to max \n");
    fprintf(stderr,"  -xl                take logs of x values before interpolating \n");
    fprintf(stderr,"  -yl                take logs of y values before interpolating \n");
    exit();
}

numeric(s) char *s;
{       
    char c;
    while(c=*s++)
    {
        if((c<='9' && c>='0') || c=='+' || c=='-' || c=='.') continue;
        return 0;
    }
    return 1;
}

curv1(x,y,yp,n) double *x,*y,*yp; 
int n;
{       
    int i;
    double c1,c2,c3,deln,delnm1,delnn,dels,delx1,delx2,delx12;
    double diag1,diag2,diagin,dx1,dx2,exps;
    double sigmap,sinhs,sinhin,slpp1,slppn,spdiag;

    delx1=x[1] - x[0];
    dx1=(y[1] - y[0])/delx1;
    if(slopes) {
        slpp1=slp1; 
        slppn=slpn;
    }
    else
    {
        if(n!=2)
        {
            delx2= x[2] - x[1];
            delx12= x[2] - x[0];
            c1= -(delx12 + delx1)/delx12/delx1;
            c2= delx12/delx1/delx2;
            c3= -delx1/delx12/delx2;
            slpp1= c1*y[0] + c2*y[1] + c3*y[2];
            deln= x[n-1] - x[n-2];
            delnm1= x[n-2] - x[n-3];
            delnn= x[n-1] - x[n-3];
            c1= (delnn + deln)/delnn/deln;
            c2= -delnn/deln/delnm1;
            c3= deln/delnn/delnm1;
            slppn= c3*y[n-3] + c2*y[n-2] + c1*y[n-1];
        }
        else yp[0]=yp[1]=0.;
    }
    /* denormalize tension factor */
    sigmap=fabs(sigma)*(n-1)/(x[n-1]-x[0]);
    /* set up right hand side and tridiagonal system for
                                                   yp and perform forward elimination                           */
    dels=sigmap*delx1;
    exps=exp(dels);
    sinhs=.5*(exps-1./exps);
    sinhin=1./(delx1*sinhs);
    diag1=sinhin*(dels*.5*(exps+1./exps)-sinhs);
    diagin=1./diag1;
    yp[0]=diagin*(dx1-slpp1);
    spdiag=sinhin*(sinhs-dels);
    temp[0]=diagin*spdiag;
    if(n!=2)
    {
        for(i=1; i<=n-2; i++)
        {
            delx2=x[i+1]-x[i];
            dx2=(y[i+1]-y[i])/delx2;
            dels=sigmap*delx2;
            exps=exp(dels);
            sinhs=.5*(exps-1./exps);
            sinhin=1./(delx2*sinhs);
            diag2=sinhin*(dels*(.5*(exps+1./exps))-sinhs);
            diagin=1./(diag1+diag2-spdiag*temp[i-1]);
            yp[i]=diagin*(dx2-dx1-spdiag*yp[i-1]);
            spdiag=sinhin*(sinhs-dels);
            temp[i]=diagin*spdiag;
            dx1=dx2;
            diag1=diag2;
        }
    }
    diagin=1./(diag1-spdiag*temp[n-2]);
    yp[n-1]=diagin*(slppn-dx2-spdiag*yp[n-2]);
    /* perform back substitution */
    for (i=n-2; i>=0; i--) yp[i] -= temp[i]*yp[i+1];
}


double curv2(x,y,yp,t,n) double *x,*y,*yp,t; 
int n;
{       
    int i,j;
    static int i1=1;
    double del1,del2,dels,exps,exps1,s,sigmap,sinhd1,sinhd2,sinhs;

    s=x[n-1]-x[0];
    sigmap=fabs(sigma)*(n-1)/s;
#ifdef WORK
    for (j=2; j; j--)               /* want: x[i-1] <= t < x[i], 0 < i <= n */
    {
        for (i=i1; i<n; i++)
        {
            if(x[i]>t) break;
        }
        if(i==n) i=n-1;
        if(x[i-1]<=t || t<=x[0]) break;
        i1=1;
    }
#endif
    i=i1;
    if(i>n) i=1;
    while(i<n && t>=x[i]) i++;
    while(i>1 && x[i-1]>t) i--;
    i1=i;
    del1=t-x[i-1];
    del2=x[i]-t;
    dels=x[i] - x[i-1];
    exps1=exp(sigmap*del1); 
    sinhd1=.5*(exps1-1./exps1);
    exps= exp(sigmap*del2); 
    sinhd2=.5*(exps-1./exps);
    exps=exps1*exps;        
    sinhs=.5*(exps-1./exps);
    return ((yp[i]*sinhd1 + yp[i-1]*sinhd2)/sinhs +
        ((y[i] - yp[i])*del1 + (y[i-1] - yp[i-1])*del2)/dels);
}

double mylog(x) double x;
{       
    if(x>0.) return log(x);
    fprintf(stderr,"can%'t take log of nonpositive number");
    exit(1);
    return 0.;
}
