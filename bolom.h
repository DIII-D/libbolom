#ifdef EFIT3365
#define YLEN 65
#define XLEN 33
#endif
#ifdef EFIT6565
#define YLEN 65
#define XLEN 65
#endif
#define XMIN 84.0
#define XMAX 254.0
#define YMIN -160.0
#define YMAX 160.0
#define CHANS 71
typedef float image[YLEN * XLEN];
#ifdef linux
#define sigvector sigvec
#endif
#define CORE
#define LOWER
#define UPPER
void efitgeom(float *,float,float,float,float,float,float,float *)

