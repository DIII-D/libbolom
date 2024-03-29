 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef long unsigned int size_t ;
 
 
 
 
 
typedef void * __gnuc_va_list ;
typedef char * va_list ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef unsigned char __u_char ;
typedef unsigned short __u_short ;
typedef unsigned int __u_int ;
typedef unsigned long __u_long ;
typedef unsigned long long int __u_quad_t ;
typedef long long int __quad_t ;
typedef signed char __int8_t ;
typedef unsigned char __uint8_t ;
typedef signed short int __int16_t ;
typedef unsigned short int __uint16_t ;
typedef signed int __int32_t ;
typedef unsigned int __uint32_t ;
typedef signed long long int __int64_t ;
typedef unsigned long long int __uint64_t ;
typedef __quad_t * __qaddr_t ;
typedef __u_quad_t __dev_t ;	
typedef __u_int __uid_t ;	
typedef __u_int __gid_t ;	
typedef __u_long __ino_t ;	
typedef __u_int __mode_t ;	
typedef __u_int __nlink_t ; 
typedef long int __off_t ;	
typedef __quad_t __loff_t ;	
typedef int __pid_t ;	
typedef int __ssize_t ;	
typedef long int __rlim_t ;	
typedef __quad_t __rlim64_t ;	
typedef __u_int __id_t ;	
typedef struct
  {
    int __val [ 2 ] ;
  } __fsid_t ;	
 
typedef int __daddr_t ;	
typedef char * __caddr_t ;
typedef long int __time_t ;
typedef long int __swblk_t ;	

typedef long int __clock_t ;
 
typedef unsigned long int __fd_mask ;
 
 
 
typedef struct
  {
    
 
    __fd_mask __fds_bits [ 1024 / ( 8 * sizeof ( __fd_mask ) ) ] ;
  } __fd_set ;
typedef int __key_t ;
 
typedef unsigned short int __ipc_pid_t ;
 
 
typedef __u_long __blkcnt_t ;
typedef __u_quad_t __blkcnt64_t ;
 
typedef long int __fsblkcnt_t ;
typedef __quad_t __fsblkcnt64_t ;
 
typedef __u_long __fsfilcnt_t ;
typedef __u_quad_t __fsfilcnt64_t ;
 
typedef __u_long __ino64_t ;
 
typedef __loff_t __off64_t ;
 
typedef int __t_scalar_t ;
typedef unsigned int __t_uscalar_t ;
 
typedef int __intptr_t ;
 
 
typedef struct _IO_FILE FILE ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef int wchar_t ;
 
 
 
typedef unsigned int wint_t ;
typedef int _G_int16_t ;
typedef int _G_int32_t ;
typedef unsigned int _G_uint16_t ;
typedef unsigned int _G_uint32_t ;
 
 
 
 
 
 
 
 
struct _IO_jump_t ; struct _IO_FILE ;
 
typedef void _IO_lock_t ;
 
struct _IO_marker {
  struct _IO_marker * _next ;
  struct _IO_FILE * _sbuf ;
  
 
   
  int _pos ;
} ;
struct _IO_FILE {
  int _flags ;	
   
   
  char * _IO_read_ptr ;	
  char * _IO_read_end ;	
  char * _IO_read_base ;	
  char * _IO_write_base ;	
  char * _IO_write_ptr ;	
  char * _IO_write_end ;	
  char * _IO_buf_base ;	
  char * _IO_buf_end ;	
   
  char * _IO_save_base ; 
  char * _IO_backup_base ; 
  char * _IO_save_end ; 
  struct _IO_marker * _markers ;
  struct _IO_FILE * _chain ;
  int _fileno ;
  int _blksize ;
  __off_t _old_offset ; 
   
  unsigned short _cur_column ;
  signed char _vtable_offset ;
  char _shortbuf [ 1 ] ;
   
  _IO_lock_t * _lock ;
  __off64_t _offset ;
   
  int _unused2 [ 16 ] ;
} ;
typedef struct _IO_FILE _IO_FILE ;
struct _IO_FILE_plus ;
extern struct _IO_FILE_plus _IO_2_1_stdin_ ;
extern struct _IO_FILE_plus _IO_2_1_stdout_ ;
extern struct _IO_FILE_plus _IO_2_1_stderr_ ;
 
 
typedef __ssize_t __io_read_fn ( void * __cookie , char * __buf , size_t __nbytes ) ;
 
typedef __ssize_t __io_write_fn ( void * __cookie , const char * __buf , size_t __n ) ;
 
typedef int __io_seek_fn ( void * __cookie , __off_t __pos , int __w ) ;
 
typedef int __io_close_fn ( void * __cookie ) ;
extern int __underflow ( _IO_FILE * ) ;
extern int __uflow ( _IO_FILE * ) ;
extern int __overflow ( _IO_FILE * , int ) ;
extern int _IO_getc ( _IO_FILE * __fp ) ;
extern int _IO_putc ( int __c , _IO_FILE * __fp ) ;
extern int _IO_feof ( _IO_FILE * __fp ) ;
extern int _IO_ferror ( _IO_FILE * __fp ) ;
extern int _IO_peekc_locked ( _IO_FILE * __fp ) ;
 
extern void _IO_flockfile ( _IO_FILE * ) ;
extern void _IO_funlockfile ( _IO_FILE * ) ;
extern int _IO_ftrylockfile ( _IO_FILE * ) ;
extern int _IO_vfscanf ( _IO_FILE * , const char * , __gnuc_va_list , int * ) ;
extern int _IO_vfprintf ( _IO_FILE * , const char * , __gnuc_va_list ) ;
extern __ssize_t _IO_padn ( _IO_FILE * , int , __ssize_t ) ;
extern size_t _IO_sgetn ( _IO_FILE * , void * , size_t ) ;
extern __off64_t _IO_seekoff ( _IO_FILE * , __off64_t , int , int ) ;
extern __off64_t _IO_seekpos ( _IO_FILE * , __off64_t , int ) ;
extern void _IO_free_backup_area ( _IO_FILE * ) ;
 
typedef __off_t fpos_t ;
 
 
 
 
 
 
 
 
extern FILE * stdin ;	
extern FILE * stdout ;	
extern FILE * stderr ;	
 
 
extern int remove ( const char * __filename ) ;
 
extern int rename ( const char * __old , const char * __new ) ;
 
extern FILE * tmpfile ( void ) ;
 
extern char * tmpnam ( char * __s ) ;
 
extern char * tmpnam_r ( char * __s ) ;
 
extern char * tempnam ( const char * __dir , const char * __pfx ) ;
 
extern int fclose ( FILE * __stream ) ;
 
extern int fflush ( FILE * __stream ) ;
 
extern int fflush_unlocked ( FILE * __stream ) ;
 
extern FILE * fopen ( const char * __filename , const char * __modes ) ;
 
extern FILE * freopen ( const char * __filename , const char * __modes , FILE * __stream ) ;
 
extern FILE * fdopen ( int __fd , const char * __modes ) ;
 
extern void setbuf ( FILE * __stream , char * __buf ) ;
 
extern int setvbuf ( FILE * __stream , char * __buf , int __modes , size_t __n ) ;
 
extern void setbuffer ( FILE * __stream , char * __buf , size_t __size ) ;
 
extern void setlinebuf ( FILE * __stream ) ;
 
extern int fprintf ( FILE * __stream , const char * __format , ... ) ;
 
extern int printf ( const char * __format , ... ) ;
 
extern int sprintf ( char * __s , const char * __format , ... ) ;
 
extern int vfprintf ( FILE * __s , const char * __format , __gnuc_va_list __arg ) ;
 
extern int vprintf ( const char * __format , __gnuc_va_list __arg ) ;
 
extern int vsprintf ( char * __s , const char * __format , __gnuc_va_list __arg ) ;
 
extern int snprintf ( char * __s , size_t __maxlen , const char * __format , ... )
     ;
extern int __vsnprintf ( char * __s , size_t __maxlen , const char * __format , __gnuc_va_list __arg )
     ;
extern int vsnprintf ( char * __s , size_t __maxlen , const char * __format , __gnuc_va_list __arg )
     ;
 
extern int fscanf ( FILE * __stream , const char * __format , ... ) ;
 
extern int scanf ( const char * __format , ... ) ;
 
extern int sscanf ( const char * __s , const char * __format , ... ) ;
 
extern int fgetc ( FILE * __stream ) ;
extern int getc ( FILE * __stream ) ;
 
extern int getchar ( void ) ;
 
 
extern int getc_unlocked ( FILE * __stream ) ;
extern int getchar_unlocked ( void ) ;
 
extern int fgetc_unlocked ( FILE * __stream ) ;
 
extern int fputc ( int __c , FILE * __stream ) ;
extern int putc ( int __c , FILE * __stream ) ;
 
extern int putchar ( int __c ) ;
 
 
extern int fputc_unlocked ( int __c , FILE * __stream ) ;
 
extern int putc_unlocked ( int __c , FILE * __stream ) ;
extern int putchar_unlocked ( int __c ) ;
 
extern int getw ( FILE * __stream ) ;
 
extern int putw ( int __w , FILE * __stream ) ;
 
extern char * fgets ( char * __s , int __n , FILE * __stream ) ;
 
extern char * gets ( char * __s ) ;
 
extern int fputs ( const char * __s , FILE * __stream ) ;
 
extern int puts ( const char * __s ) ;
 
extern int ungetc ( int __c , FILE * __stream ) ;
 
extern size_t fread ( void * __ptr , size_t __size , size_t __n , FILE * __stream ) ;
 
extern size_t fwrite ( const void * __ptr , size_t __size , size_t __n , FILE * __s ) ;
 
extern size_t fread_unlocked ( void * __ptr , size_t __size , size_t __n , FILE * __stream ) ;
extern size_t fwrite_unlocked ( const void * __ptr , size_t __size , size_t __n , FILE * __stream ) ;
 
extern int fseek ( FILE * __stream , long int __off , int __whence ) ;
 
extern long int ftell ( FILE * __stream ) ;
 
extern void rewind ( FILE * __stream ) ;
 
 
typedef __off_t off_t ;
 
extern int fgetpos ( FILE * __stream , fpos_t * __pos ) ;
 
extern int fsetpos ( FILE * __stream , const fpos_t * __pos ) ;
 
extern void clearerr ( FILE * __stream ) ;
 
extern int feof ( FILE * __stream ) ;
 
extern int ferror ( FILE * __stream ) ;
 
extern void clearerr_unlocked ( FILE * __stream ) ;
extern int feof_unlocked ( FILE * __stream ) ;
extern int ferror_unlocked ( FILE * __stream ) ;
 
extern void perror ( const char * __s ) ;
 
extern int sys_nerr ;
extern const char * const sys_errlist [ ] ;
 
extern int fileno ( FILE * __stream ) ;
 
extern int fileno_unlocked ( FILE * __stream ) ;
 
extern FILE * popen ( const char * __command , const char * __modes ) ;
 
extern int pclose ( FILE * __stream ) ;
 
extern char * ctermid ( char * __s ) ;
 
 
extern void flockfile ( FILE * __stream ) ;
 
extern int ftrylockfile ( FILE * __stream ) ;
 
extern void funlockfile ( FILE * __stream ) ;
 
 
 
 
 
 
 
 
 
 
 
 
static union { unsigned char __c [ 8 ] ; double __d ; } __huge_val = { { 0 , 0 , 0 , 0 , 0 , 0 , 0xf0 , 0x7f } } ;
 
 
 
 
 
 
 
 
 
extern double acos ( double __x ) ; extern double __acos ( double __x ) ;
 
extern double asin ( double __x ) ; extern double __asin ( double __x ) ;
 
extern double atan ( double __x ) ; extern double __atan ( double __x ) ;
 
extern double atan2 ( double __y , double __x ) ; extern double __atan2 ( double __y , double __x ) ;
 
extern double cos ( double __x ) ; extern double __cos ( double __x ) ;
 
extern double sin ( double __x ) ; extern double __sin ( double __x ) ;
 
extern double tan ( double __x ) ; extern double __tan ( double __x ) ;
 
 
extern double cosh ( double __x ) ; extern double __cosh ( double __x ) ;
 
extern double sinh ( double __x ) ; extern double __sinh ( double __x ) ;
 
extern double tanh ( double __x ) ; extern double __tanh ( double __x ) ;
 
extern double acosh ( double __x ) ; extern double __acosh ( double __x ) ;
 
extern double asinh ( double __x ) ; extern double __asinh ( double __x ) ;
 
extern double atanh ( double __x ) ; extern double __atanh ( double __x ) ;
 
 
extern double exp ( double __x ) ; extern double __exp ( double __x ) ;
 
extern double frexp ( double __x , int * __exponent ) ; extern double __frexp ( double __x , int * __exponent ) ;
 
extern double ldexp ( double __x , int __exponent ) ; extern double __ldexp ( double __x , int __exponent ) ;
 
extern double log ( double __x ) ; extern double __log ( double __x ) ;
 
extern double log10 ( double __x ) ; extern double __log10 ( double __x ) ;
 
extern double modf ( double __x , double * __iptr ) ; extern double __modf ( double __x , double * __iptr ) ;
 
extern double expm1 ( double __x ) ; extern double __expm1 ( double __x ) ;
 
extern double log1p ( double __x ) ; extern double __log1p ( double __x ) ;
 
extern double logb ( double __x ) ; extern double __logb ( double __x ) ;
 
 
extern double pow ( double __x , double __y ) ; extern double __pow ( double __x , double __y ) ;
 
extern double sqrt ( double __x ) ; extern double __sqrt ( double __x ) ;
 
extern double hypot ( double __x , double __y ) ; extern double __hypot ( double __x , double __y ) ;
 
extern double cbrt ( double __x ) ; extern double __cbrt ( double __x ) ;
 
 
extern double ceil ( double __x ) ; extern double __ceil ( double __x ) ;
 
extern double fabs ( double __x ) ; extern double __fabs ( double __x ) ;
 
extern double floor ( double __x ) ; extern double __floor ( double __x ) ;
 
extern double fmod ( double __x , double __y ) ; extern double __fmod ( double __x , double __y ) ;
 
extern int __isinf ( double __value ) ;
 
extern int __finite ( double __value ) ;
 
extern int isinf ( double __value ) ;
 
extern int finite ( double __value ) ;
 
extern double infnan ( int __error ) ; extern double __infnan ( int __error ) ;
 
extern double drem ( double __x , double __y ) ; extern double __drem ( double __x , double __y ) ;
 
extern double significand ( double __x ) ; extern double __significand ( double __x ) ;
 
extern double copysign ( double __x , double __y ) ; extern double __copysign ( double __x , double __y ) ;
 
extern int __isnan ( double __value ) ;
 
extern int isnan ( double __value ) ;
 
extern double j0 ( double ) ; extern double __j0 ( double ) ;
extern double j1 ( double ) ; extern double __j1 ( double ) ;
extern double jn ( int , double ) ; extern double __jn ( int , double ) ;
extern double y0 ( double ) ; extern double __y0 ( double ) ;
extern double y1 ( double ) ; extern double __y1 ( double ) ;
extern double yn ( int , double ) ; extern double __yn ( int , double ) ;
 
extern double erf ( double ) ; extern double __erf ( double ) ;
extern double erfc ( double ) ; extern double __erfc ( double ) ;
extern double lgamma ( double ) ; extern double __lgamma ( double ) ;
extern double tgamma ( double ) ; extern double __tgamma ( double ) ;
 
extern double gamma ( double ) ; extern double __gamma ( double ) ;
 
extern double lgamma_r ( double , int * __signgamp ) ; extern double __lgamma_r ( double , int * __signgamp ) ;
 
extern double rint ( double __x ) ; extern double __rint ( double __x ) ;
 
extern double nextafter ( double __x , double __y ) ; extern double __nextafter ( double __x , double __y ) ;
 
extern double remainder ( double __x , double __y ) ; extern double __remainder ( double __x , double __y ) ;
 
extern double scalb ( double __x , double __n ) ; extern double __scalb ( double __x , double __n ) ;
 
extern double scalbn ( double __x , int __n ) ; extern double __scalbn ( double __x , int __n ) ;
 
extern int ilogb ( double __x ) ; extern int __ilogb ( double __x ) ;
 
 
 
 
 
extern float acosf ( float __x ) ; extern float __acosf ( float __x ) ;
 
extern float asinf ( float __x ) ; extern float __asinf ( float __x ) ;
 
extern float atanf ( float __x ) ; extern float __atanf ( float __x ) ;
 
extern float atan2f ( float __y , float __x ) ; extern float __atan2f ( float __y , float __x ) ;
 
extern float cosf ( float __x ) ; extern float __cosf ( float __x ) ;
 
extern float sinf ( float __x ) ; extern float __sinf ( float __x ) ;
 
extern float tanf ( float __x ) ; extern float __tanf ( float __x ) ;
 
 
extern float coshf ( float __x ) ; extern float __coshf ( float __x ) ;
 
extern float sinhf ( float __x ) ; extern float __sinhf ( float __x ) ;
 
extern float tanhf ( float __x ) ; extern float __tanhf ( float __x ) ;
 
extern float acoshf ( float __x ) ; extern float __acoshf ( float __x ) ;
 
extern float asinhf ( float __x ) ; extern float __asinhf ( float __x ) ;
 
extern float atanhf ( float __x ) ; extern float __atanhf ( float __x ) ;
 
 
extern float expf ( float __x ) ; extern float __expf ( float __x ) ;
 
extern float frexpf ( float __x , int * __exponent ) ; extern float __frexpf ( float __x , int * __exponent ) ;
 
extern float ldexpf ( float __x , int __exponent ) ; extern float __ldexpf ( float __x , int __exponent ) ;
 
extern float logf ( float __x ) ; extern float __logf ( float __x ) ;
 
extern float log10f ( float __x ) ; extern float __log10f ( float __x ) ;
 
extern float modff ( float __x , float * __iptr ) ; extern float __modff ( float __x , float * __iptr ) ;
 
extern float expm1f ( float __x ) ; extern float __expm1f ( float __x ) ;
 
extern float log1pf ( float __x ) ; extern float __log1pf ( float __x ) ;
 
extern float logbf ( float __x ) ; extern float __logbf ( float __x ) ;
 
 
extern float powf ( float __x , float __y ) ; extern float __powf ( float __x , float __y ) ;
 
extern float sqrtf ( float __x ) ; extern float __sqrtf ( float __x ) ;
 
extern float hypotf ( float __x , float __y ) ; extern float __hypotf ( float __x , float __y ) ;
 
extern float cbrtf ( float __x ) ; extern float __cbrtf ( float __x ) ;
 
 
extern float ceilf ( float __x ) ; extern float __ceilf ( float __x ) ;
 
extern float fabsf ( float __x ) ; extern float __fabsf ( float __x ) ;
 
extern float floorf ( float __x ) ; extern float __floorf ( float __x ) ;
 
extern float fmodf ( float __x , float __y ) ; extern float __fmodf ( float __x , float __y ) ;
 
extern int __isinff ( float __value ) ;
 
extern int __finitef ( float __value ) ;
 
extern int isinff ( float __value ) ;
 
extern int finitef ( float __value ) ;
 
extern float infnanf ( int __error ) ; extern float __infnanf ( int __error ) ;
 
extern float dremf ( float __x , float __y ) ; extern float __dremf ( float __x , float __y ) ;
 
extern float significandf ( float __x ) ; extern float __significandf ( float __x ) ;
 
extern float copysignf ( float __x , float __y ) ; extern float __copysignf ( float __x , float __y ) ;
 
extern int __isnanf ( float __value ) ;
 
extern int isnanf ( float __value ) ;
 
extern float j0f ( float ) ; extern float __j0f ( float ) ;
extern float j1f ( float ) ; extern float __j1f ( float ) ;
extern float jnf ( int , float ) ; extern float __jnf ( int , float ) ;
extern float y0f ( float ) ; extern float __y0f ( float ) ;
extern float y1f ( float ) ; extern float __y1f ( float ) ;
extern float ynf ( int , float ) ; extern float __ynf ( int , float ) ;
 
extern float erff ( float ) ; extern float __erff ( float ) ;
extern float erfcf ( float ) ; extern float __erfcf ( float ) ;
extern float lgammaf ( float ) ; extern float __lgammaf ( float ) ;
extern float tgammaf ( float ) ; extern float __tgammaf ( float ) ;
 
extern float gammaf ( float ) ; extern float __gammaf ( float ) ;
 
extern float lgammaf_r ( float , int * __signgamp ) ; extern float __lgammaf_r ( float , int * __signgamp ) ;
 
extern float rintf ( float __x ) ; extern float __rintf ( float __x ) ;
 
extern float nextafterf ( float __x , float __y ) ; extern float __nextafterf ( float __x , float __y ) ;
 
extern float remainderf ( float __x , float __y ) ; extern float __remainderf ( float __x , float __y ) ;
 
extern float scalbf ( float __x , float __n ) ; extern float __scalbf ( float __x , float __n ) ;
 
extern float scalbnf ( float __x , int __n ) ; extern float __scalbnf ( float __x , int __n ) ;
 
extern int ilogbf ( float __x ) ; extern int __ilogbf ( float __x ) ;
 
 
 
 
 
extern long double acosl ( long double __x ) ; extern long double __acosl ( long double __x ) ;
 
extern long double asinl ( long double __x ) ; extern long double __asinl ( long double __x ) ;
 
extern long double atanl ( long double __x ) ; extern long double __atanl ( long double __x ) ;
 
extern long double atan2l ( long double __y , long double __x ) ; extern long double __atan2l ( long double __y , long double __x ) ;
 
extern long double cosl ( long double __x ) ; extern long double __cosl ( long double __x ) ;
 
extern long double sinl ( long double __x ) ; extern long double __sinl ( long double __x ) ;
 
extern long double tanl ( long double __x ) ; extern long double __tanl ( long double __x ) ;
 
 
extern long double coshl ( long double __x ) ; extern long double __coshl ( long double __x ) ;
 
extern long double sinhl ( long double __x ) ; extern long double __sinhl ( long double __x ) ;
 
extern long double tanhl ( long double __x ) ; extern long double __tanhl ( long double __x ) ;
 
extern long double acoshl ( long double __x ) ; extern long double __acoshl ( long double __x ) ;
 
extern long double asinhl ( long double __x ) ; extern long double __asinhl ( long double __x ) ;
 
extern long double atanhl ( long double __x ) ; extern long double __atanhl ( long double __x ) ;
 
 
extern long double expl ( long double __x ) ; extern long double __expl ( long double __x ) ;
 
extern long double frexpl ( long double __x , int * __exponent ) ; extern long double __frexpl ( long double __x , int * __exponent ) ;
 
extern long double ldexpl ( long double __x , int __exponent ) ; extern long double __ldexpl ( long double __x , int __exponent ) ;
 
extern long double logl ( long double __x ) ; extern long double __logl ( long double __x ) ;
 
extern long double log10l ( long double __x ) ; extern long double __log10l ( long double __x ) ;
 
extern long double modfl ( long double __x , long double * __iptr ) ; extern long double __modfl ( long double __x , long double * __iptr ) ;
 
extern long double expm1l ( long double __x ) ; extern long double __expm1l ( long double __x ) ;
 
extern long double log1pl ( long double __x ) ; extern long double __log1pl ( long double __x ) ;
 
extern long double logbl ( long double __x ) ; extern long double __logbl ( long double __x ) ;
 
 
extern long double powl ( long double __x , long double __y ) ; extern long double __powl ( long double __x , long double __y ) ;
 
extern long double sqrtl ( long double __x ) ; extern long double __sqrtl ( long double __x ) ;
 
extern long double hypotl ( long double __x , long double __y ) ; extern long double __hypotl ( long double __x , long double __y ) ;
 
extern long double cbrtl ( long double __x ) ; extern long double __cbrtl ( long double __x ) ;
 
 
extern long double ceill ( long double __x ) ; extern long double __ceill ( long double __x ) ;
 
extern long double fabsl ( long double __x ) ; extern long double __fabsl ( long double __x ) ;
 
extern long double floorl ( long double __x ) ; extern long double __floorl ( long double __x ) ;
 
extern long double fmodl ( long double __x , long double __y ) ; extern long double __fmodl ( long double __x , long double __y ) ;
 
extern int __isinfl ( long double __value ) ;
 
extern int __finitel ( long double __value ) ;
 
extern int isinfl ( long double __value ) ;
 
extern int finitel ( long double __value ) ;
 
extern long double infnanl ( int __error ) ; extern long double __infnanl ( int __error ) ;
 
extern long double dreml ( long double __x , long double __y ) ; extern long double __dreml ( long double __x , long double __y ) ;
 
extern long double significandl ( long double __x ) ; extern long double __significandl ( long double __x ) ;
 
extern long double copysignl ( long double __x , long double __y ) ; extern long double __copysignl ( long double __x , long double __y ) ;
 
extern int __isnanl ( long double __value ) ;
 
extern int isnanl ( long double __value ) ;
 
extern long double j0l ( long double ) ; extern long double __j0l ( long double ) ;
extern long double j1l ( long double ) ; extern long double __j1l ( long double ) ;
extern long double jnl ( int , long double ) ; extern long double __jnl ( int , long double ) ;
extern long double y0l ( long double ) ; extern long double __y0l ( long double ) ;
extern long double y1l ( long double ) ; extern long double __y1l ( long double ) ;
extern long double ynl ( int , long double ) ; extern long double __ynl ( int , long double ) ;
 
extern long double erfl ( long double ) ; extern long double __erfl ( long double ) ;
extern long double erfcl ( long double ) ; extern long double __erfcl ( long double ) ;
extern long double lgammal ( long double ) ; extern long double __lgammal ( long double ) ;
extern long double tgammal ( long double ) ; extern long double __tgammal ( long double ) ;
 
extern long double gammal ( long double ) ; extern long double __gammal ( long double ) ;
 
extern long double lgammal_r ( long double , int * __signgamp ) ; extern long double __lgammal_r ( long double , int * __signgamp ) ;
 
extern long double rintl ( long double __x ) ; extern long double __rintl ( long double __x ) ;
 
extern long double nextafterl ( long double __x , long double __y ) ; extern long double __nextafterl ( long double __x , long double __y ) ;
 
extern long double remainderl ( long double __x , long double __y ) ; extern long double __remainderl ( long double __x , long double __y ) ;
 
extern long double scalbl ( long double __x , long double __n ) ; extern long double __scalbl ( long double __x , long double __n ) ;
 
extern long double scalbnl ( long double __x , int __n ) ; extern long double __scalbnl ( long double __x , int __n ) ;
 
extern int ilogbl ( long double __x ) ; extern int __ilogbl ( long double __x ) ;
 
extern int signgam ;
 
 
typedef enum
{
  _IEEE_ = - 1 ,	
  _SVID_ ,	
  _XOPEN_ ,	
  _POSIX_ ,
  _ISOC_	
} _LIB_VERSION_TYPE ;
 
extern _LIB_VERSION_TYPE _LIB_VERSION ;
 
struct exception
  {
    int type ;
    char * name ;
    double arg1 ;
    double arg2 ;
    double retval ;
  } ;
extern int matherr ( struct exception * __exc ) ;
 
 
 
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
union __convert_long_double {
  unsigned __convert_long_double_i [ 4 ] ;
  long double __convert_long_double_d ;
} ;
    
    
    
    
    
    
 
 
 
 
 
 
int debug = 0 ;
typedef float image [ 65 * 33 ] ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef struct
  {
    int quot ;	
    int rem ;	
  } div_t ;
 
typedef struct
  {
    long int quot ;	
    long int rem ;	
  } ldiv_t ;
 
 
 
extern size_t __ctype_get_mb_cur_max ( void ) ;
 
extern double atof ( const char * __nptr ) ;
 
extern int atoi ( const char * __nptr ) ;
 
extern long int atol ( const char * __nptr ) ;
 
extern double strtod ( const char * __nptr , char * * __endptr ) ;
 
extern long int strtol ( const char * __nptr , char * * __endptr , int __base ) ;
 
extern unsigned long int strtoul ( const char * __nptr , char * * __endptr , int __base ) ;
 
extern double __strtod_internal ( const char * __nptr , char * * __endptr , int __group ) ;
extern float __strtof_internal ( const char * __nptr , char * * __endptr , int __group ) ;
extern long double __strtold_internal ( const char * __nptr , char * * __endptr , int __group ) ;
extern long int __strtol_internal ( const char * __nptr , char * * __endptr , int __base , int __group ) ;
extern unsigned long int __strtoul_internal ( const char * __nptr , char * * __endptr , int __base , int __group ) ;
 
extern char * l64a ( long int __n ) ;
 
extern long int a64l ( const char * __s ) ;
 
 
 
 
 
typedef __u_char u_char ;
typedef __u_short u_short ;
typedef __u_int u_int ;
typedef __u_long u_long ;
typedef __quad_t quad_t ;
typedef __u_quad_t u_quad_t ;
typedef __fsid_t fsid_t ;
typedef __loff_t loff_t ;
typedef __ino_t ino_t ;
typedef __dev_t dev_t ;
typedef __gid_t gid_t ;
typedef __mode_t mode_t ;
typedef __nlink_t nlink_t ;
typedef __uid_t uid_t ;
typedef __pid_t pid_t ;
typedef __id_t id_t ;
typedef __ssize_t ssize_t ;
typedef __daddr_t daddr_t ;
typedef __caddr_t caddr_t ;
typedef __key_t key_t ;
 
 
 
 
 
typedef __time_t time_t ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef unsigned long int ulong ;
typedef unsigned short int ushort ;
typedef unsigned int uint ;
 
 
typedef	char int8_t ;
typedef	short int int16_t ;
typedef	int int32_t ;
 typedef long long int int64_t ;
 
typedef	unsigned char u_int8_t ;
typedef	unsigned short int u_int16_t ;
typedef	unsigned int u_int32_t ;
 typedef unsigned long long int u_int64_t ;
typedef int register_t ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
typedef int __sig_atomic_t ;
 
typedef struct
  {
    unsigned long int __val [ ( 1024 / ( 8 * sizeof ( unsigned long int ) ) ) ] ;
  } __sigset_t ;
 
 
 
 
 
struct timespec
  {
    long int tv_sec ;	
    long int tv_nsec ;	
  } ;
 
struct timeval ;
typedef __fd_mask fd_mask ;
 
typedef __fd_set fd_set ;
 
 
 
 
extern int __select ( int __nfds , __fd_set * __readfds , __fd_set * __writefds , __fd_set * __exceptfds , struct timeval * __timeout ) ;
extern int select ( int __nfds , __fd_set * __readfds , __fd_set * __writefds , __fd_set * __exceptfds , struct timeval * __timeout ) ;
 
 
 
 
 
 
typedef __blkcnt_t blkcnt_t ;	
typedef __fsblkcnt_t fsblkcnt_t ; 
typedef __fsfilcnt_t fsfilcnt_t ; 
 
 
extern int32_t random ( void ) ;
 
extern void srandom ( unsigned int __seed ) ;
 
extern void * initstate ( unsigned int __seed , void * __statebuf , size_t __statelen ) ;
 
extern void * setstate ( void * __statebuf ) ;
 
struct random_data
  {
    int32_t * fptr ;	
    int32_t * rptr ;	
    int32_t * state ;	
    int rand_type ;	
    int rand_deg ;	
    int rand_sep ;	
    int32_t * end_ptr ;	
  } ;
extern int random_r ( struct random_data * __buf , int32_t * __result ) ;
extern int srandom_r ( unsigned int __seed , struct random_data * __buf ) ;
extern int initstate_r ( unsigned int __seed , void * __statebuf , size_t __statelen , struct random_data * __buf ) ;
extern int setstate_r ( void * __statebuf , struct random_data * __buf ) ;
 
extern int rand ( void ) ;
 
extern void srand ( unsigned int __seed ) ;
 
extern int rand_r ( unsigned int * __seed ) ;
 
 
extern double drand48 ( void ) ;
extern double erand48 ( unsigned short int __xsubi [ 3 ] ) ;
 
extern long int lrand48 ( void ) ;
extern long int nrand48 ( unsigned short int __xsubi [ 3 ] ) ;
 
extern long int mrand48 ( void ) ;
extern long int jrand48 ( unsigned short int __xsubi [ 3 ] ) ;
 
extern void srand48 ( long int __seedval ) ;
extern unsigned short int * seed48 ( unsigned short int __seed16v [ 3 ] ) ;
extern void lcong48 ( unsigned short int __param [ 7 ] ) ;
 
struct drand48_data
  {
    unsigned short int x [ 3 ] ;	
    unsigned short int a [ 3 ] ;	
    unsigned short int c ;	
    unsigned short int old_x [ 3 ] ; 
    int init ;	
  } ;
 
extern int drand48_r ( struct drand48_data * __buffer , double * __result ) ;
extern int erand48_r ( unsigned short int __xsubi [ 3 ] , struct drand48_data * __buffer , double * __result ) ;
 
extern int lrand48_r ( struct drand48_data * __buffer , long int * __result ) ;
extern int nrand48_r ( unsigned short int __xsubi [ 3 ] , struct drand48_data * __buffer , long int * __result ) ;
 
extern int mrand48_r ( struct drand48_data * __buffer , long int * __result ) ;
extern int jrand48_r ( unsigned short int __xsubi [ 3 ] , struct drand48_data * __buffer , long int * __result ) ;
 
extern int srand48_r ( long int __seedval , struct drand48_data * __buffer ) ;
extern int seed48_r ( unsigned short int __seed16v [ 3 ] , struct drand48_data * __buffer ) ;
extern int lcong48_r ( unsigned short int __param [ 7 ] , struct drand48_data * __buffer ) ;
 
extern void * malloc ( size_t __size ) ;
 
extern void * calloc ( size_t __nmemb , size_t __size ) ;
 
extern void * realloc ( void * __ptr , size_t __size ) ;
 
extern void free ( void * __ptr ) ;
 
extern void cfree ( void * __ptr ) ;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
extern void * __alloca ( size_t __size ) ;
extern void * alloca ( size_t __size ) ;
extern void * __builtin_alloca ( size_t __size ) ;
 
extern void * valloc ( size_t __size ) ;
 
extern void abort ( void ) ;
 
extern int atexit ( void ( * __func ) ( void ) ) ;
 
extern int __on_exit ( void ( * __func ) ( int __status , void * __arg ) , void * __arg ) ;
extern int on_exit ( void ( * __func ) ( int __status , void * __arg ) , void * __arg ) ;
 
extern void exit ( int __status ) ;
 
extern char * getenv ( const char * __name ) ;
 
extern char * __secure_getenv ( const char * __name ) ;
 
 
extern int putenv ( char * __string ) ;
 
extern int setenv ( const char * __name , const char * __value , int __replace ) ;
 
extern void unsetenv ( const char * __name ) ;
 
extern int clearenv ( void ) ;
 
extern char * mktemp ( char * __template ) ;
 
extern int mkstemp ( char * __template ) ;
 
extern int system ( const char * __command ) ;
 
extern char * realpath ( const char * __name , char * __resolved ) ;
 
typedef int ( * __compar_fn_t ) ( const void * , const void * ) ;
 
extern void * bsearch ( const void * __key , const void * __base , size_t __nmemb , size_t __size , __compar_fn_t __compar ) ;
 
extern void qsort ( void * __base , size_t __nmemb , size_t __size , __compar_fn_t __compar ) ;
 
extern int abs ( int __x ) ;
extern long int labs ( long int __x ) ;
 
 
extern div_t div ( int __numer , int __denom ) ;
extern ldiv_t ldiv ( long int __numer , long int __denom )
     ;
 
 
extern char * ecvt ( double __value , int __ndigit , int * __decpt , int * __sign ) ;
 
extern char * fcvt ( double __value , int __ndigit , int * __decpt , int * __sign ) ;
 
extern char * gcvt ( double __value , int __ndigit , char * __buf ) ;
 
extern char * qecvt ( long double __value , int __ndigit , int * __decpt , int * __sign ) ;
extern char * qfcvt ( long double __value , int __ndigit , int * __decpt , int * __sign ) ;
extern char * qgcvt ( long double __value , int __ndigit , char * __buf ) ;
 
extern int ecvt_r ( double __value , int __ndigit , int * __decpt , int * __sign , char * __buf , size_t __len ) ;
extern int fcvt_r ( double __value , int __ndigit , int * __decpt , int * __sign , char * __buf , size_t __len ) ;
extern int qecvt_r ( long double __value , int __ndigit , int * __decpt , int * __sign , char * __buf , size_t __len ) ;
extern int qfcvt_r ( long double __value , int __ndigit , int * __decpt , int * __sign , char * __buf , size_t __len ) ;
 
extern int mblen ( const char * __s , size_t __n ) ;
 
extern int mbtowc ( wchar_t * __pwc , const char * __s , size_t __n ) ;
 
extern int wctomb ( char * __s , wchar_t __wchar ) ;
 
extern size_t mbstowcs ( wchar_t * __pwcs , const char * __s , size_t __n ) ;
 
extern size_t wcstombs ( char * __s , const wchar_t * __pwcs , size_t __n ) ;
 
extern int rpmatch ( const char * __response ) ;
 
 
 
 
 
 
 
 
extern int __sigismember ( const __sigset_t * , int ) ;
extern int __sigaddset ( __sigset_t * , int ) ;
extern int __sigdelset ( __sigset_t * , int ) ;
 
typedef __sig_atomic_t sig_atomic_t ;
typedef __sigset_t sigset_t ;
 
 
 
 
 
 
 
typedef void ( * __sighandler_t ) ( int ) ;
 
extern __sighandler_t __sysv_signal ( int __sig , __sighandler_t __handler ) ;
 
extern __sighandler_t signal ( int __sig , __sighandler_t __handler ) ;
 
extern int kill ( __pid_t __pid , int __sig ) ;
 
extern int killpg ( __pid_t __pgrp , int __sig ) ;
 
extern int raise ( int __sig ) ;
 
extern __sighandler_t ssignal ( int __sig , __sighandler_t __handler ) ;
extern int gsignal ( int __sig ) ;
 
extern void psignal ( int __sig , const char * __s ) ;
 
extern int __sigpause ( int __sig_or_mask , int __is_sig ) ;
 
extern int sigpause ( int __mask ) ;
 
 
 
extern int sigblock ( int __mask ) ;
 
extern int sigsetmask ( int __mask ) ;
 
extern int siggetmask ( void ) ;
 
typedef __sighandler_t sig_t ;
 
 
 
 
 
 
typedef union sigval
  {
    int sival_int ;
    void * sival_ptr ;
  } sigval_t ;
typedef struct siginfo
  {
    int si_signo ;	
    int si_errno ;	
 
    int si_code ;	
    union
      {
	int _pad [ ( ( 128 / sizeof ( int ) ) - 3 ) ] ;
	  
	struct
	  {
	    __pid_t si_pid ;	
	    __uid_t si_uid ;	
	  } _kill ;
	 
	struct
	  {
	    unsigned int _timer1 ;
	    unsigned int _timer2 ;
	  } _timer ;
	 
	struct
	  {
	    __pid_t si_pid ;	
	    __uid_t si_uid ;	
	    sigval_t si_sigval ;	
	  } _rt ;
	 
	struct
	  {
	    __pid_t si_pid ;	
	    __uid_t si_uid ;	
	    int si_status ;	
	    __clock_t si_utime ;
	    __clock_t si_stime ;
	  } _sigchld ;
	 
	struct
	  {
	    void * si_addr ;	
	  } _sigfault ;
	 
	struct
	  {
	    int si_band ;	
	    int si_fd ;
	  } _sigpoll ;
      } _sifields ;
  } siginfo_t ;
 
 
enum
{
  SI_SIGIO = - 5 ,	
  SI_ASYNCIO ,	
  SI_MESGQ ,	
  SI_TIMER ,	
  SI_QUEUE ,	
  SI_USER	
} ;
 
enum
{
  ILL_ILLOPC = 1 ,	
  ILL_ILLOPN ,	
  ILL_ILLADR ,	
  ILL_ILLTRP ,	
  ILL_PRVOPC ,	
  ILL_PRVREG ,	
  ILL_COPROC ,	
  ILL_BADSTK	
} ;
 
enum
{
  FPE_INTDIV = 1 ,	
  FPE_INTOVF ,	
  FPE_FLTDIV ,	
  FPE_FLTOVF ,	
  FPE_FLTUND ,	
  FPE_FLTRES ,	
  FPE_FLTINV ,	
  FPE_FLTSUB	
} ;
 
enum
{
  SEGV_MAPERR = 1 ,	
  SEGV_ACCERR	
} ;
 
enum
{
  BUS_ADRALN = 1 ,	
  BUS_ADRERR ,	
  BUS_OBJERR	
} ;
 
enum
{
  TRAP_BRKPT = 1 ,	
  TRAP_TRACE	
} ;
 
enum
{
  CLD_EXITED = 1 ,	
  CLD_KILLED ,	
  CLD_DUMPED ,	
  CLD_TRAPPED ,	
  CLD_STOPPED ,	
  CLD_CONTINUED	
} ;
 
enum
{
  POLL_IN = 1 ,	
  POLL_OUT ,	
  POLL_MSG ,	
  POLL_ERR ,	
  POLL_PRI ,	
  POLL_HUP	
} ;
 
typedef struct sigevent
  {
    sigval_t sigev_value ;
    int sigev_signo ;
    int sigev_notify ;
    union
      {
	int _pad [ ( ( 64 / sizeof ( int ) ) - 3 ) ] ;
	struct
	  {
	    void ( * _function ) ( sigval_t ) ; 
	    void * _attribute ;	
	  } _sigev_thread ;
      } _sigev_un ;
  } sigevent_t ;
 
 
enum
{
  SIGEV_SIGNAL = 0 ,	
  SIGEV_NONE ,	
  SIGEV_THREAD	
} ;
 
extern int sigemptyset ( sigset_t * __set ) ;
 
extern int sigfillset ( sigset_t * __set ) ;
 
extern int sigaddset ( sigset_t * __set , int __signo ) ;
 
extern int sigdelset ( sigset_t * __set , int __signo ) ;
 
extern int sigismember ( const sigset_t * __set , int __signo ) ;
 
 
 
struct sigaction
  {
     
    union
      {
	 
	__sighandler_t sa_handler ;
	 
	void ( * sa_sigaction ) ( int , siginfo_t * , void * ) ;
      }
    __sigaction_handler ;
     
    __sigset_t sa_mask ;
     
    int sa_flags ;
     
    void ( * sa_restorer ) ( void ) ;
  } ;
 
 
 
 
extern int sigprocmask ( int __how , const sigset_t * __set , sigset_t * __oset ) ;
 
extern int sigsuspend ( const sigset_t * __set ) ;
 
extern int __sigaction ( int __sig , const struct sigaction * __act , struct sigaction * __oact ) ;
extern int sigaction ( int __sig , const struct sigaction * __act , struct sigaction * __oact ) ;
 
extern int sigpending ( sigset_t * __set ) ;
 
extern int sigwait ( const sigset_t * __set , int * __sig ) ;
 
extern int sigwaitinfo ( const sigset_t * __set , siginfo_t * __info ) ;
 
extern int sigtimedwait ( const sigset_t * __set , siginfo_t * __info , const struct timespec * __timeout ) ;
 
extern int sigqueue ( __pid_t __pid , int __sig , const union sigval __val ) ;
 
extern const char * const _sys_siglist [ 64 ] ;
extern const char * const sys_siglist [ 64 ] ;
 
struct sigvec
  {
    __sighandler_t sv_handler ;	
    int sv_mask ;	
    int sv_flags ;	
  } ;
 
 
extern int sigvec ( int __sig , const struct sigvec * __vec , struct sigvec * __ovec ) ;
 
 
 
 
struct _fpreg {
	unsigned short significand [ 4 ] ;
	unsigned short exponent ;
} ;
struct _fpstate {
	unsigned long cw ,
			sw ,
			tag ,
			ipoff ,
			cssel ,
			dataoff ,
			datasel ;
	struct _fpreg	_st [ 8 ] ;
	unsigned long	status ;
} ;
struct sigcontext {
	unsigned short gs , __gsh ;
	unsigned short fs , __fsh ;
	unsigned short es , __esh ;
	unsigned short ds , __dsh ;
	unsigned long edi ;
	unsigned long esi ;
	unsigned long ebp ;
	unsigned long esp ;
	unsigned long ebx ;
	unsigned long edx ;
	unsigned long ecx ;
	unsigned long eax ;
	unsigned long trapno ;
	unsigned long err ;
	unsigned long eip ;
	unsigned short cs , __csh ;
	unsigned long eflags ;
	unsigned long esp_at_signal ;
	unsigned short ss , __ssh ;
	struct _fpstate * fpstate ;
	unsigned long oldmask ;
	unsigned long cr2 ;
} ;
 
extern int sigreturn ( struct sigcontext * __scp ) ;
 
extern int siginterrupt ( int __sig , int __interrupt ) ;
 
 
struct sigstack
  {
    void * ss_sp ;	
    int ss_onstack ;	
  } ;
 
enum
{
  SS_ONSTACK = 1 ,
  SS_DISABLE
} ;
 
 
 
typedef struct sigaltstack
  {
    void * ss_sp ;
    int ss_flags ;
    size_t ss_size ;
  } stack_t ;
 
extern int sigstack ( struct sigstack * __ss , struct sigstack * __oss ) ;
 
extern int sigaltstack ( const struct sigaltstack * __ss , struct sigaltstack * __oss ) ;
 
 
extern int __libc_current_sigrtmin ( void ) ;
 
extern int __libc_current_sigrtmax ( void ) ;
 
 
extern int bound [ ] ;
extern int bound_len ;
extern int bound_flag [ ] ;
 
 
 
 
 
typedef int __jmp_buf [ 6 ] ;
 
 
 
 
typedef struct __jmp_buf_tag	
  {
    
 
    __jmp_buf __jmpbuf ;	
    int __mask_was_saved ;	
    __sigset_t __saved_mask ;	
  } jmp_buf [ 1 ] ;
 
extern int __sigsetjmp ( jmp_buf __env , int __savemask ) ;
 
 
 
extern void longjmp ( jmp_buf __env , int __val )
     ;
 
extern void _longjmp ( jmp_buf __env , int __val )
     ;
 
typedef jmp_buf sigjmp_buf ;
 
 
extern void siglongjmp ( sigjmp_buf __env , int __val )
     ;
static jmp_buf sjbuf ;
static int ilxp , iuxp ;
static float xinc , yinc ;
float * gmatrix ;
int * gpos ;
int gmatrix_len ;
int gpos_len ;
static void ( * prev_handler ) ( ) ;
static float * undump ;
static float * insigma ;
static float * fitpts ;
static int nfitpts ;
static image inefit ;
static int done ;
bolomfit_core ( shot , chans , nchans , inproj , sigma , kpsi , nkpsi , tens , efit , rax , zax ,
rxpt1 , zxpt1 , rxpt2 , zxpt2 , fitproj , fitimage , errorimage , fret )
unsigned int shot ;
float * inproj , * sigma , * kpsi , * fitproj , * fret , tens ;
image efit , fitimage , errorimage ;
int * chans , nchans , nkpsi ;
float rax , zax , rxpt1 , zxpt1 , rxpt2 , zxpt2 ;
{
    float * p , * * xi ;
    int iter ;
    int splcore ( ) ;
    FILE * f ;
    int i , j , k ;
    float T ;
    float c1 , c2 , c3 , c4 ;
    float * amat , * cmat , * imat , * timat ;
    float sum ;
    float * work ;
    int na , nma , eqn ;
    int l_m , l_n , l_p , l_lda , l_ldb ;
    int lwork , info ;
    float tproj [ 71 ] ;
    float basisel ( ) , basisv ( ) ;
    int numunk ;
    int iparm ;
    float * g , * h , * w ;
    int mdw , me , ma , mg ;
    int ip [ 4096 ] ;
    float * ws ;
    float prgopt [ 100 ] , rnorme , rnorml ;
    int mode ;
    int allbound ;
    float * u ;
    float * tz , * uz ;
    float * sol ;
    struct sigvec vec , ovec ;
		if ( ! bolom_set_gmatrix ( shot ) ) return ( 0 ) ;
    done = 0 ;
    __sigsetjmp ( ( sjbuf ) , 0 ) ;
    if ( done ) return ;
    done = 1 ;
    vec . sv_handler = splcore ;
    vec . sv_mask = 0xffff ;
    vec . sv_flags = 0 ;
    vec . sv_flags = 0 ;
    sigvector ( 2 , & vec , & ovec ) ;
    T = tens ;
    prev_handler = ovec . sv_handler ;
     
    for ( i = 0 ; i < 71 ; ++ i ) tproj [ i ] = 0.0 ;
    for ( i = 0 ; i < nchans ; ++ i ) {
        tproj [ chans [ i ] ] = inproj [ chans [ i ] ] ;
    }
     
    xinc = ( 254.0 - 84.0 ) / ( float ) ( 33 - 1 ) ;
    yinc = ( 160.0 - - 160.0 ) / ( float ) ( 65 - 1 ) ;
    if ( rxpt1 < 84.0 || rxpt1 > 254.0 ||
        zxpt1 < - 160.0 || zxpt1 > 160.0 ) {
        ilxp = 0 ;
    } else ilxp = ( zxpt1 - - 160.0 ) / yinc ;
    if ( rxpt2 < 84.0 || rxpt2 > 254.0 ||
        zxpt2 < - 160.0 || zxpt2 > 160.0 ) {
        iuxp = 65 - 1 ;
    } else iuxp = ( zxpt2 - - 160.0 ) / yinc ;
    ilxp *= 33 ;
    iuxp *= 33 ;
    sigvector ( 2 , & ovec , 0 ) ;
     
    na = 33 * 65 ;
     
    numunk = 2 * nkpsi ;
    cmat = calloc ( sizeof ( float ) , numunk * na ) ;
    for ( k = 0 ; k < na ; ++ k ) {
        for ( j = 0 ; j < numunk ; ++ j ) {
            iparm = j + 1 ;
            cmat [ k * numunk + j ] = basisel ( & iparm , ( float * ) & ( efit [ k ] ) ,
                & T , kpsi , & nkpsi ) ;
        }
    }
    imat = calloc ( sizeof ( float ) * 2 , numunk * nchans ) ;
    timat = calloc ( sizeof ( float ) * 2 , numunk * nchans ) ;
    amat = calloc ( sizeof ( float ) * 2 , na * ( 71 + 1 ) ) ;
    for ( i = 0 , j = 0 , k = 1 ; i < gmatrix_len ; ++ i , j += 2 , k += 2 ) {
        if ( ( gpos [ j ] ) < na && ( gpos [ k ] ) < 71 ) {
            if ( bound_flag [ ( gpos [ j ] ) ] ) continue ;
            if ( ( gpos [ j ] ) < ilxp || ( gpos [ j ] ) > iuxp ) continue ;
            amat [ ( gpos [ k ] ) * na + ( gpos [ j ] ) ] = gmatrix [ i ] ;
        }
    }
    for ( i = 0 ; i < nchans ; ++ i ) {
        for ( j = 0 ; j < numunk ; ++ j ) {
            for ( k = 0 , sum = 0.0 ; k < na ; ++ k ) sum += amat [ chans [ i ] * na + k ] * cmat [ k * numunk + j ] ;
            imat [ i * numunk + j ] = sum ; 
        }
    }
    eqn = 3 * nkpsi + 6 ;
    u = calloc ( sizeof ( float ) * 2 , eqn * numunk ) ;
    uz = calloc ( sizeof ( float ) * 2 , eqn ) ;
    nma = 0 ;
    bacnst ( & numunk , & T , kpsi , & nkpsi , u , & eqn , & nma , uz ) ;
     
     
    sol = calloc ( sizeof ( float ) * 2 , numunk ) ;
    tz = calloc ( sizeof ( float ) * 2 , nchans ) ;
    for ( i = 0 ; i < nchans ; ++ i ) tz [ i ] = inproj [ chans [ i ] ] ;
    for ( i = 0 ; i < nchans ; ++ i ) {
        for ( j = 0 ; j < numunk ; ++ j ) {
            timat [ j * nchans + i ] = imat [ i * numunk + j ] ;
        }
    }
    lwork = 10 * ( numunk + nma + nchans ) ;
    work = calloc ( sizeof ( float ) * 2 , lwork ) ;
    l_m = nchans ;
    l_n = numunk ;
    l_p = nma ;
    l_lda = nchans ;
    l_ldb = eqn ;
    mg = numunk / 2 ;
    me = nma ;
    ma = nchans ;
    g = calloc ( sizeof ( float ) , mg * numunk ) ;
    h = calloc ( sizeof ( float ) , mg ) ;
    for ( j = 0 ; j < mg ; ++ j ) {
        g [ j * numunk + 2 * j ] = 1.0 ;
    }
    mdw = me + ma + mg ;
    w = calloc ( sizeof ( float ) , mdw * ( numunk + 1 ) ) ;
    if ( me > 0 ) {
        for ( i = 0 ; i < me ; ++ i ) {
            for ( j = 0 ; j < numunk ; ++ j ) {
                w [ i + j * mdw ] = u [ i + j * eqn ] ;
            }
        }
        for ( i = 0 ; i < me ; ++ i ) {
            w [ i + numunk * mdw ] = uz [ i ] ;
        }
    }
    if ( ma > 0 ) {
        for ( i = 0 ; i < ma ; ++ i ) {
            for ( j = 0 ; j < numunk ; ++ j ) {
                w [ i + me + j * mdw ] = imat [ j + i * numunk ] ;
            }
        }
        for ( i = 0 ; i < ma ; ++ i ) {
            w [ i + me + numunk * mdw ] = tz [ i ] ;
        }
    }
    if ( mg > 0 ) {
        for ( i = 0 ; i < mg ; ++ i ) {
            for ( j = 0 ; j < numunk ; ++ j ) {
                w [ i + me + ma + j * mdw ] = g [ i * numunk + j ] ;
            }
        }
        for ( i = 0 ; i < mg ; ++ i ) {
            w [ i + me + ma + numunk * mdw ] = h [ i ] ;
        }
    }
    ws = calloc ( sizeof ( float ) , 128000 ) ;
    ip [ 0 ] = 128000 ;
    ip [ 1 ] = 4096 ;
    memset ( ip , 0 , sizeof ( ip ) ) ;
    prgopt [ 0 ] = 4 ;
    prgopt [ 1 ] = 4 ;
    prgopt [ 2 ] = 1.0e-4 ;
    prgopt [ 3 ] = 7 ;
    prgopt [ 4 ] = 5 ;
    prgopt [ 5 ] = 1.0e-4 ;
    prgopt [ 6 ] = 9 ;
    prgopt [ 7 ] = 1 ;
    prgopt [ 8 ] = 1 ;
    lsei ( w , & mdw , & me , & ma , & mg , & numunk , prgopt , sol , & rnorme , & rnorml , & mode , ws , ip ) ;
    memset ( fitimage , 0 , 33 * 65 * sizeof ( float ) ) ;
    for ( k = 0 ; k < na ; ++ k ) {
        if ( bound_flag [ k ] ) continue ;
        if ( k < ilxp || k > iuxp ) continue ;
        if ( efit [ k ] >= kpsi [ 0 ] && efit [ k ] <= kpsi [ nkpsi - 1 ] )
            fitimage [ k ] = basisv ( & numunk , & efit [ k ] , & T , sol ,
                kpsi , & nkpsi ) ;
        else fitimage [ k ] = 0.0 ;
         
    }
    memset ( errorimage , 0 , 33 * 65 * sizeof ( float ) ) ;
    for ( i = 0 ; i < numunk ; ++ i ) {
				if ( w [ i + i * mdw ] >= 0.0 ) sol [ i ] = __builtin_sqrt ( w [ i + i * mdw ] ) ;
				else sol [ i ] = 0.0 ;
    }
		T = 10.0 ;
    for ( k = 0 ; k < na ; ++ k ) {
        if ( bound_flag [ k ] ) continue ;
        if ( k < ilxp || k > iuxp ) continue ;
        if ( efit [ k ] >= kpsi [ 0 ] && efit [ k ] <= kpsi [ nkpsi - 1 ] )
            errorimage [ k ] = basisv ( & numunk , & efit [ k ] , & T , sol ,
                kpsi , & nkpsi ) ;
        else errorimage [ k ] = 0.0 ;
         
    }
     
    bolomclip ( shot , fitimage , chans , nchans ) ;
    bolomproj ( shot , fitimage , fitproj ) ;
    * fret = 0.0 ;
    for ( i = 0 ; i < nchans ; ++ i ) {
        * fret += __builtin_pow ( ( inproj [ chans [ i ] ] - fitproj [ chans [ i ] ] ) , 2.0 ) / __builtin_pow ( sigma [ chans [ i ] ] , 2.0 ) ;
    }
    if ( nchans > numunk )
        * fret = * fret / ( float ) ( nchans - numunk ) ;
     
    free ( cmat ) ;
    free ( imat ) ;
    free ( timat ) ;
    free ( amat ) ;
    free ( u ) ;
    free ( uz ) ;
    free ( sol ) ;
    free ( tz ) ;
    free ( work ) ;
    free ( g ) ;
    free ( h ) ;
    free ( w ) ;
    free ( ws ) ;
}
splcore ( )
{
    FILE * f ;
    int i ;
    image gimage ;
    struct sigvec vec , ovec ;
    printf ( "dumping unknowns in corecoefs.dat\n" ) ;
    f = fopen ( "corecoefs.dat" , "w" ) ;
    for ( i = 1 ; i <= nfitpts ; ++ i ) {
        fprintf ( f , "%g\n" , undump [ i ] ) ;
    }
    fflush ( f ) ;
    close ( f ) ;
    if ( prev_handler ) {
        vec . sv_handler = prev_handler ;
        vec . sv_mask = 0xffff ;
        vec . sv_flags = 0 ;
        vec . sv_flags = 0 ;
        sigvector ( 2 , & vec , 0 ) ;
        prev_handler ( ) ;
        longjmp ( sjbuf , 0 ) ;
    } else
        exit ( 0 ) ;
}
 
c ( )
{
}
