/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McXtrace <http://www.mcxtrace.org>
 * Instrument: Monochromator_testing.instr (Monochroamtor_testing)
 * Date:       Wed Oct 24 12:51:32 2018
 * File:       Monochromator_testing.c
 * Compile:    cc -o Monochroamtor_testing.out Monochromator_testing.c  -lgsl -lgslcblas
 * CFLAGS= -lgsl -lgslcblas
 */


#define MCCODE_STRING "McXtrace 1.4 - Jun. 09, 2017"
#define FLAVOR "mcxtrace"
#define FLAVOR_UPPER "MCXTRACE"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McXtrace 1.4
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */



/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McXtrace 1.4 - Jun. 09, 2017"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Jun. 09, 2017"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "1.4"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McXtrace"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "xray"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCXTRACE"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */

/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
#  define PI 3.14159265358979323846
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic inline double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic inline void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic inline void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 691 "Monochromator_testing.c"

#line 1 "mcxtrace-r.h"
/*******************************************************************************
*
* McXtrace, X-ray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcxtrace-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McXtrace X.Y
* Version: $Revision$
*
* Runtime system header for McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
*******************************************************************************/

#ifndef MCXTRACE_R_H
#define MCXTRACE_R_H "$Revision$"

/* Following part is only embedded when not redundant with mcstas.h ========= */

#ifndef MCCODE_H

#define CELE     1.602176487e-19   /* [C] Elementary charge CODATA 2006*/
#define M_C      299792458         /* [m/s] speed of light CODATA 2006*/
#define E2K      0.506773091264796 /* Convert k[1/AA] to E [keV] (CELE/(HBAR*M_C)*1e-10)*1e3 */
#define K2E      1.97326972808327  /*Convert E[keV] to k[1/AA] (1e10*M_C*HBAR/CELE)/1e3 */
#define RE       2.8179402894e-5   /*[AA] Thomson scattering length*/

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlkx, mcnlky, mcnlkz, \
    mcnlphi, mcnlt, mcnlEx,mcnlEy,mcnlEz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlkx, mcnlky, mcnlkz, \
    mcnlphi, mcnlt, mcnlEx,mcnlEy,mcnlEz, mcnlp); mcDEBUG_ABSORB(); goto mcabsorb;} while(0)

#define STORE_XRAY(index, x,y,z, kx,ky,kz, phi, t, Ex,Ey,Ez, p) \
  mcstore_xray(mccomp_storein,index, x,y,z, kx,ky,kz, phi, t, Ex,Ey,Ez, p);
#define RESTORE_XRAY(index, x,y,z, kx,ky,kz, phi, t, Ex,Ey,Ez, p) \
  mcrestore_xray(mccomp_storein,index, &x,&y,&z, &kx,&ky,&kz, &phi, &t, &Ex,&Ey,&Ez, &p);

/*magnet stuff is probably redundant*/
#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
    /* change coordinates from local system to magnet system */ \
    Rotation rotLM, rotTemp; \
    Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
    rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
    rot_mul(rotTemp, mcMagnetRot, rotLM); \
    mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
	   	       &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
  } while(0)

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/*An interrupt a'la mcMagnet should be inserted below if there's non-zero permeability*/
/*perhaps some kind of PROP_POL*/

#define mcPROP_DL(dl) \
  do { \
    MCNUM k=sqrt( scalar_prod(mcnlkx,mcnlky,mcnlkz,mcnlkx,mcnlky,mcnlkz));\
    mcnlx += (dl)*mcnlkx/k;\
    mcnly += (dl)*mcnlky/k;\
    mcnlz += (dl)*mcnlkz/k;\
    mcnlphi += 1e10*k*(dl);\
    mcnlt += (dl)/((double)M_C);\
  }while (0)

/*gravity not an issue with x-rays*/
/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration. */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)

/*adapted from PROP_DT(dt)*//*{{{*/
#define PROP_DL(dl) \
  do{ \
    if( dl <0 && mcallowbackprop == 0) { (mcAbsorbProp[INDEX_CURRENT_COMP])++; ABSORB; }; \
    mcPROP_DL(dl); \
    DISALLOW_BACKPROP;\
  } while (0)

#define PROP_DT(dt) \
  do { \
    if(dt < 0 ) { RESTORE=1; goto mcabsorbComp; };		    \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)/*}}}*/

#define PROP_Z0 \
  mcPROP_P0(z)

#define PROP_X0 \
  mcPROP_P0(x)

#define PROP_Y0 \
  mcPROP_P0(y)

#define mcPROP_P0(P) \
  do { \
    MCNUM mc_dl,mc_k; \
    if(mcnlk ## P == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_k=sqrt(scalar_prod(mcnlkx,mcnlky,mcnlkz,mcnlkx,mcnlky,mcnlkz));\
    mc_dl= -mcnl ## P * mc_k / mcnlk ## P;\
    if(mc_dl<0 && mcallowbackprop==0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; };\
    PROP_DL(mc_dl);\
  } while(0)

void mcsetstate(double x, double y, double z, double kx, double ky, double kz,
    double phi, double t, double Ex, double Ey, double Ez, double p);


#endif /* !MCCODE_H */


#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
      x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p);
#define mcDEBUG_SCATTER(x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
      x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p);

#else

#define mcDEBUG_STATE(x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p)
#define mcDEBUG_SCATTER(x,y,z,kx,ky,kz,phi,t,Ex,Ey,Ez,p)

#endif


#endif /* MCXTRACE_R_H */
/* End of file "mcxtrace-r.h". */

#line 886 "Monochromator_testing.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;
  
  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);
      
    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+"); 
    
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n", 
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m || !detector.filename)
    return(detector);
  
  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];
  
  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m) 
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin); 
      else x = 0;
      if (detector.n && detector.p) 
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin); 
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw) 
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */
      
      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }
      
      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
    
    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH, 
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }
  
  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  
  return(detector);
  
} /* mcdetector_statistics */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=abs(m); n=abs(n); p=abs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre, 
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
    if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
      if (mcinputtable[i].par == NULL)
        strncpy(Parameters, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
      else
        (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);

      fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
    }
  }
  fflush(f);
} /* mcruninfo_out */

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;
  
  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" : 
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f, 
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" : 
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre, 
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
    
  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  /* Write data set information to simulation description file. */
  MPI_MASTER(
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n", 
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);
      
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
  
}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
        fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      }
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2, 
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0, 
          outfile, detector.istransposed);
      }
      fclose(outfile);
      
      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",  
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=32; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;
  
  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);
  
  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);
  
    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }
  
  return(ret);
  
} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "r");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: mcsiminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];
  
  if (!f || mcdisable_output_files) return;
  
  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */ 
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */
      
      NXclosegroup(f); /* simulation */
    } /* NXsimulation */
    
    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  if (!f || !detector.m || mcdisable_output_files) return;
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
    
    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
    
      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" : 
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" : 
                 "xylimits", detector.limits);
      nxprintattr(f, "variables", 
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");
        
      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */
  
} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[32];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);
    
    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{
  
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;
  
  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);
  
  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;
  
  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */
  
  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",  
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }
  
  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  
  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
  
      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar, 
          1, detector.m, detector.xmin, detector.xmax);
          
        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar, 
          2, detector.n, detector.ymin, detector.ymax);
          
        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar, 
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */
      
      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);
      
      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may. 
*   Then the number of recv/send is not constant along nodes, and simulation stalls.  
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );
  
  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
	    
	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );   
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );
  
  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  
#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
{
  int exists=0;
  int index;
  
  /* check format */      
  if (!mcformat || !strlen(mcformat) 
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE") 
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }
  
  /* open the SIM file if not defined yet */
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));
    
} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));
  
} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];
  
  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (abs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (abs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));
  
} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;
  
  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,abs(m),1,abs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);
  
  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out("  ", stdout);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  printf("MCDISPLAY: magnify('%s')\n", what);
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic inline void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic inline double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir)) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  
#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif
  srandom(mcseed);

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 4818 "Monochromator_testing.c"

#line 1 "mcxtrace-r.c"
/*******************************************************************************
*
* McXtrace, X-ray tracing package
*           Copyright (C) 1997-2009, All rights reserved
*           Risoe National Laboratory, Roskilde, Denmark
*           Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcxtrace-r.c
*
* %Identification
* Edited by: EK
* Date:    May 29, 2009
* Release: McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McXtrace.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embedded in the c code whenever required.
*
*******************************************************************************/

#ifndef MCXTRACE_H

/*******************************************************************************
* mcstore_xray: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_xray(MCNUM *s, int index, double x, double y, double z,
               double kx, double ky, double kz, double phi, double t,
               double Ex, double Ey, double Ez, double p)
{
    double *dptr = &s[12*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = kx;
    *dptr++  = ky;
    *dptr++  = kz;
    *dptr++  = phi;
    *dptr++  = t;
    *dptr++  = Ex;
    *dptr++  = Ey;
    *dptr++  = Ez;
    *dptr    = p ;
}

/*******************************************************************************
* mcrestore_xray: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_xray(MCNUM *s, int index, double *x, double *y, double *z,
               double *kx, double *ky, double *kz, double *phi, double *t,
               double *Ex, double *Ey, double *Ez, double *p)
{
    double *dptr = &s[12*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *kx =  *dptr++;
    *ky =  *dptr++;
    *kz =  *dptr++;
    *phi=  *dptr++;
    *t  =  *dptr++;
    *Ex =  *dptr++;
    *Ey =  *dptr++;
    *Ez =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_xray */

/*******************************************************************************
* mcsetstate: transfer parameters into global McXtrace variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double kx, double ky, double kz,
           double phi, double t, double Ex, double Ey, double Ez, double p)
{
  extern double mcnx, mcny, mcnz, mcnkx, mcnky, mcnkz;
  extern double mcnphi, mcnt, mcnEx, mcnEy, mcnEz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnkx = kx;
  mcnky = ky;
  mcnkz = kz;
  mcnphi = phi;
  mcnt = t;
  mcnEx = Ex;
  mcnEy = Ey;
  mcnEz = Ez;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default xray parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute length intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting travelling lengths dl_in and dl_out
*******************************************************************************/
int box_intersect(double *dl_in, double *dl_out,
                  double x, double y, double z,
                  double kx, double ky, double kz,
                  double dx, double dy, double dz)
{

  double k, l,xf,yf,zf, l_[6],dx_2,dy_2,dz_2;
  double ab[2];
  unsigned int count=0;
  k=sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));
  dx_2=dx/2.0;dy_2=dy/2.0;dz_2=dz/2.0; 
  /*we really don't need to store the 6 intersects as only two are possible. i.e. should remove that.*/
  if (kx) {
    l=(-dx_2-x)/kx*k;
    yf=l*ky/k+y;zf=l*kz/k+z;
    if(yf > -dy_2 && yf<dy_2 && zf > -dz_2 && zf<dz_2){
      l_[0]=l;
      ab[count++]=l_[0];
    }else{
      l_[0]=0;
    }
    l=(dx_2-x)/kx*k;
    yf=l*ky/k+y;zf=l*kz/k+z;
    if(yf > -dy_2 && yf<dy_2 && zf > -dz_2 && zf<dz_2){
      l_[1]=l;
      ab[count++]=l_[1];
    }else{
      l_[1]=0;
    }
  }
  if (ky) {
    l=(-dy_2-y)/ky*k;
    xf=l*kx/k+x;zf=l*kz/k+z;
    if(xf > -dx_2 && xf<dx_2 && zf > -dz_2 && zf<dz_2){
      l_[2]=l;
      ab[count++]=l_[2];
    }else{
      l_[2]=0;
    } 
    l=(dy_2-y)/ky*k;
    xf=l*kx/k+x;zf=l*kz/k+z;
    if(xf > -dx_2 && xf<dx_2 && zf > -dz_2 && zf<dz_2){
      l_[3]=l;
      ab[count++]=l_[3];
    }else{
      l_[3]=0;
    }
  }
  if (kz) {
    l=(-dz_2-z)/kz*k;
    xf=l*kx/k+x; yf=l*ky/k+y;
    if(xf > -dx_2 && xf<dx_2 && yf > -dy_2 && yf<dy_2){
      l_[4]=l;
      ab[count++]=l_[4];
    }else{
      l_[4]=0;
    }
    l=(dz_2-z)/kz*k;
    xf=l*kx/k+x; yf=l*ky/k+y;
    if(xf > -dx_2 && xf<dx_2 && yf > -dy_2 && yf<dy_2){
      l_[5]=l;
      ab[count++]=l_[5];
    }else{
      l_[5]=0;
    }
  }
  /*check validity of intersects*/
  if (count>2){
    fprintf(stderr,"box_instersect: xray hitting box more than twice\n");
  }
  if (!count){
    *dl_in=0;*dl_out=0;
    return 0;
  }

  if (ab[0]<ab[1]){
    *dl_in=ab[0];*dl_out=ab[1];
    return 1;
  }else{
    *dl_in=ab[1];*dl_out=ab[0];
    return 1;
  }
} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 1/2/4/8/16 bits depending on intersection,
 *     and resulting times l0 and l1
 * Written by: EK 11.6.09 
 *******************************************************************************/
int
cylinder_intersect(double *l0, double *l1, double x, double y, double z,
                   double kx, double ky, double kz, double r, double h)
{
  double A,B,C,D,k2,k;
  double dl1p=0,dl0p=0,dl1c=0,dl0c=0,y0,y1;
  int ret=1,stat=0,plane_stat=0;
  enum {HIT_CYL=01,ENTER_TOP=02,ENTER_BOT=04,EXIT_TOP=010,EXIT_BOT=020,ENTER_MASK=06,EXIT_MASK=030};
  k2=(kx*kx + ky*ky + kz*kz);
  k=sqrt(k2);

  /*check for prop. vector 0*/
  if(!k2) return 0;

  A= (k2 - ky*ky);
  B= 2*(x*kx + z*kz);
  C=(x*x + z*z - r*r);
  D=B*B-4*A*C;
  if(D>=0){
    if (kx || kz){
      stat|=HIT_CYL;
    /*propagation not parallel to y-axis*/
    /*hit infinitely high cylinder?*/
      D=sqrt(D);
      dl0c=k*(-B-D)/(2*A);
      dl1c=k*(-B+D)/(2*A);
      y0=dl0c*ky/k+y;
      y1=dl1c*ky/k+y;
      if ( (y0<-h/2 && y1<-h/2) || (y0>h/2 && y1>h/2) ){
        /*ray passes above or below cylinder*/
        return 0;
      }
    }
    /*now check top and bottom planes*/
    if (ky){
      dl0p = k*(-h/2-y)/ky;
      dl1p = k*(h/2-y)/ky;
      /*switch solutions?*/
      if (dl0p<dl1p){
        plane_stat|=(ENTER_BOT|EXIT_TOP);
      }else{
        double tmp=dl1p;
        dl1p=dl0p;dl0p=tmp;
        plane_stat|=(ENTER_TOP|EXIT_BOT);
      }
    }
  }
  if (stat & HIT_CYL){
    if (ky && dl0p>dl0c){
      *l0=dl0p;/*1st top/bottom plane intersection happens after 1st cylinder intersect*/
      stat|= plane_stat & ENTER_MASK;
    } else
      *l0=dl0c;
    if(ky && dl1p<dl1c){
      *l1=dl1p;/*2nd top/bottom plane intersection happens before 2nd cylinder intersect*/
      stat|= plane_stat & EXIT_MASK;
    }else
      *l1=dl1c;
  }
  return stat;
} /* cylinder_intersect */

/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting lengths l0 and l1 
 *******************************************************************************/
int
sphere_intersect(double *l0, double *l1, double x, double y, double z,
                 double kx, double ky, double kz, double r)
{
  double B, C, D, k;

  k = kx*kx + ky*ky + kz*kz;
  B = (x*kx + y*ky + z*kz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - k*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *l0 = (-B - D) / sqrt(k);
  *l1 = (-B + D) / sqrt(k);
  return 1;
} /* sphere_intersect */

/******************************************************************************
 * ellipsoid_intersect: Calculate intersection between a line and an ellipsoid.
 * They ellisoid is fixed by a set of half-axis (a,b,c) and a matrix Q, with the
 * columns of Q being the (orthogonal) vectors along which the half-axis lie.
 * This allows for complete freedom in orienting th eellipsoid.
 * returns 0 when no intersection is found
 *      or 1 when they are found with resulting lemngths l0 and l1.
 *****************************************************************************/
int
ellipsoid_intersect(double *l0, double *l1, double x, double y, double z,
    double kx, double ky, double kz, double a, double b, double c,
    Rotation Q)
{
  Rotation A,Gamma,Q_t,Tmp;
  double u,v,w;

  Gamma[0][0]=Gamma[0][1]=Gamma[0][2]=0;
  Gamma[1][1]=Gamma[1][0]=Gamma[1][2]=0;
  Gamma[2][2]=Gamma[2][0]=Gamma[2][1]=0;
  /*now set diagonal to ellipsoid half axis if non-zero.
   * This way a zero value mean the sllipsoid extends infinitely along that axis,
   * which is useful for objects only curved in one direction*/ 
  if (a!=0){
    Gamma[0][0]=1/(a*a);
  }
  if (b!=0){
    Gamma[1][1]=1/(b*b);
  }
  if (c!=0){
    Gamma[2][2]=1/(c*c);
  }

  if (Q!=NULL){
    rot_transpose(Q,Q_t);
    rot_mul(Gamma,Q_t,Tmp);
    rot_mul(Q,Tmp,A);
  }else{
    rot_copy(A,Gamma);
  }

  /*to get the solutions as lengths in m use unit vector along k*/
  double ex,ey,ez,k;
  k=sqrt(kx*kx+ky*ky+kz*kz);
  ex=kx/k;
  ey=ky/k;
  ez=kz/k;

  u=ex*(A[0][0]*ex + A[1][0]*ey + A[2][0]*ez) + ey*( A[0][1]*ex + A[1][1]*ey + A[2][1]*ez) + ez*(A[0][2]*ex + A[1][2]*ey + A[2][2]*ez);
  v=x *(A[0][0]*ex + A[1][0]*ey + A[2][0]*ez) + ex*(A[0][0]*x + A[1][0]*y + A[2][0]*z) +
    y *(A[0][1]*ex + A[1][1]*ey + A[2][1]*ez) + ey*(A[0][1]*x + A[1][1]*y + A[2][1]*z) +
    z *(A[0][2]*ex + A[1][2]*ey + A[2][2]*ez) + ez*(A[0][2]*x + A[1][2]*y + A[2][2]*z);
  w=x*(A[0][0]*x + A[1][0]*y + A[2][0]*z) + y*(A[0][1]*x + A[1][1]*y + A[2][1]*z) + z*(A[0][2]*x + A[1][2]*y + A[2][2]*z);

  double D=v*v-4*u*w+4*u;
  if (D<0) return 0;

  D=sqrt(D);

  *l0=(-v-D) / (2*u);
  *l1=(-v+D) / (2*u);
  return 1;
}


/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane (with normal n including the point w)
 * and a line through x along the direction k.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection length is positive and negative, respectively
 *******************************************************************************/
int
plane_intersect(double *l, double x, double y, double z,
                 double kx, double ky, double kz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s,k2;
  k2=scalar_prod(kx,ky,kz,kx,ky,kz);
  s=scalar_prod(kx,ky,kz,nx,ny,nz);
  if (k2<FLT_EPSILON || fabs(s)<FLT_EPSILON) return 0;
  *l = - sqrt(k2)*scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*l<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCXTRACE_H */
/* End of file "mcxtrace-r.c". */

#line 5206 "Monochromator_testing.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Monochroamtor_testing";
char mcinstrument_source[] = "Monochromator_testing.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Source_flat'. */
#line 58 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 5052 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value shoud just be used as
     * if the table had been read. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we shoudl proceed with freeing the memory
     * associated with the table item - otherwise do nothing since there are more references that may need it.*/ 

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref=(t_Table *)calloc(1,sizeof(t_Table));
            /*copy the contents of the table handle*/
            *(tr->table_ref)= *((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found - no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found - move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            return (void *)0x1 ;/*item not found*/ 
    } 

}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && mcinstrument_source && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array+CHAR_BUF_LENGTH;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index < Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1,j);
  Y2 = Table_Index(Table,Index  ,j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* fisrt allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /* if t_Table array is not long enough, expand and realocate */
      if (block_number >= allocated-1) {
        allocated += 256;
        Table_Array = (t_Table *)realloc(Table_Array,
           allocated*sizeof(t_Table));
        if (!Table_Array) {
          fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
              allocated*sizeof(t_Table));
          *blocks = 0;
          return (NULL);
        }
      }
      /* store it into t_Table array */
      //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
      Table_Array[block_number-1] = Table;
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index=0;
    if (!Table) return;
    while (Table[index].data || Table[index].header){
            Table_Free(&Table[index]);
            index++;
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */

#line 6655 "Monochromator_testing.c"

/* Shared user declarations for all components 'Undulator'. */
#line 67 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

double mxundulator_Bsig_integrand(double x, void *params){
    double w_w1 = *((double *) params);
    double p = *((double *) params+1);
    double q = *((double *) params+2);
    double angle_term = *((double *) params+3);/* xi*gamma/K  horizontal angle term*/

    double f1 = angle_term - cos(x);
    double inner= w_w1*x - p*sin(x) + q*sin(2*x);
    double f2 =  cos(inner);

    return f1*f2;
}

double mxundulator_Bpi_integrand(double x, void *params){
    double w_w1 = *((double *) params);
    double p = *((double *) params+1);
    double q = *((double *) params+2);
    double angle_term = *((double *) params+3); /*phi*gamma/K vertical angle term*/

    double f1 = angle_term;
    double inner= w_w1*x - p*sin(x) + q*sin(2*x);
    double f2 =  cos(inner);

    return f1*f2;
}

double mxundulator_S_N(double w_w1, int N){
    return pow(sin(N*M_PI*w_w1)/(N*sin(M_PI*w_w1)),2.0);
}


#line 6695 "Monochromator_testing.c"

/* Instrument parameters. */
MCNUM mcipr_rho;
MCNUM mcipm;
MCNUM mcipE0;
MCNUM mcipdE;
MCNUM mcipblazed;
MCNUM mcipblazed_angle;
MCNUM mcipcff;
MCNUM mcipWanted_energy;
MCNUM mcipSourceChoice;

#define mcNUMIPAR 9
int mcnumipar = 9;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "r_rho", &mcipr_rho, instr_type_double, "800", 
  "m", &mcipm, instr_type_double, "6", 
  "E0", &mcipE0, instr_type_double, "0.09", 
  "dE", &mcipdE, instr_type_double, "0.03", 
  "blazed", &mcipblazed, instr_type_double, "1", 
  "blazed_angle", &mcipblazed_angle, instr_type_double, "2", 
  "cff", &mcipcff, instr_type_double, "2.25", 
  "Wanted_energy", &mcipWanted_energy, instr_type_double, "0.08", 
  "SourceChoice", &mcipSourceChoice, instr_type_double, "0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Monochroamtor_testing
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaMonochroamtor_testing coords_set(0,0,0)
#define r_rho mcipr_rho
#define m mcipm
#define E0 mcipE0
#define dE mcipdE
#define blazed mcipblazed
#define blazed_angle mcipblazed_angle
#define cff mcipcff
#define Wanted_energy mcipWanted_energy
#define SourceChoice mcipSourceChoice
#line 39 "Monochromator_testing.instr"
   double mirror2_angle, angle_grating;
   double Wanted_energy,m;
   double MCAngleVariation;
   double A,B,C,D,E,F,G,X; 
   double Lmax, Lmin, Emax, Emin;
   double E1st, Nper, undK;
#line 6744 "Monochromator_testing.c"
#undef SourceChoice
#undef Wanted_energy
#undef cff
#undef blazed_angle
#undef blazed
#undef dE
#undef E0
#undef m
#undef r_rho
#undef mcposaMonochroamtor_testing
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* xray state table at each component input (local coords) */
/* [x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p] */
MCNUM mccomp_storein[12*22];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[22];
Coords mccomp_posr[22];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[22];
MCNUM  mcPCounter[22];
MCNUM  mcP2Counter[22];
#define mcNUMCOMP 21 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[22];
/* Flag true when previous component acted on the xray (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when xray should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Definition parameters for component 'origin' [1]. */
#define mccorigin_profile 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'origin' [1]. */
MCNUM mccorigin_percent;
MCNUM mccorigin_flag_save;
MCNUM mccorigin_minutes;

/* Definition parameters for component 'source_flat' [2]. */
#define mccsource_flat_spectrum_file NULL /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'source_flat' [2]. */
MCNUM mccsource_flat_radius;
MCNUM mccsource_flat_yheight;
MCNUM mccsource_flat_xwidth;
MCNUM mccsource_flat_xmin;
MCNUM mccsource_flat_xmax;
MCNUM mccsource_flat_dist;
MCNUM mccsource_flat_focus_xw;
MCNUM mccsource_flat_focus_yh;
MCNUM mccsource_flat_E0;
MCNUM mccsource_flat_dE;
MCNUM mccsource_flat_lambda0;
MCNUM mccsource_flat_dlambda;
MCNUM mccsource_flat_flux;
MCNUM mccsource_flat_gauss;
MCNUM mccsource_flat_randomphase;
MCNUM mccsource_flat_phase;

/* Definition parameters for component 'dmu' [3]. */
#define mccdmu_verbose 1
/* Setting parameters for component 'dmu' [3]. */
MCNUM mccdmu_E0;
MCNUM mccdmu_dE;
MCNUM mccdmu_phase;
MCNUM mccdmu_randomphase;
MCNUM mccdmu_Ee;
MCNUM mccdmu_dEe;
MCNUM mccdmu_Ie;
MCNUM mccdmu_tbunch;
MCNUM mccdmu_t0;
MCNUM mccdmu_B;
MCNUM mccdmu_K;
MCNUM mccdmu_gap;
int mccdmu_Nper;
MCNUM mccdmu_lu;
MCNUM mccdmu_sigey;
MCNUM mccdmu_sigex;
MCNUM mccdmu_sigepx;
MCNUM mccdmu_sigepy;
MCNUM mccdmu_focus_xw;
MCNUM mccdmu_focus_yh;
MCNUM mccdmu_dist;
MCNUM mccdmu_gauss_t;
MCNUM mccdmu_quick_integ;
MCNUM mccdmu_E1st;

/* Definition parameters for component 'Mbefore_E_monitor' [4]. */
#define mccMbefore_E_monitor_nE 20
#define mccMbefore_E_monitor_filename "Mbefore_E_monitor" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Mbefore_E_monitor' [4]. */
MCNUM mccMbefore_E_monitor_xmin;
MCNUM mccMbefore_E_monitor_xmax;
MCNUM mccMbefore_E_monitor_ymin;
MCNUM mccMbefore_E_monitor_ymax;
MCNUM mccMbefore_E_monitor_xwidth;
MCNUM mccMbefore_E_monitor_yheight;
MCNUM mccMbefore_E_monitor_Emin;
MCNUM mccMbefore_E_monitor_Emax;
MCNUM mccMbefore_E_monitor_restore_xray;

/* Definition parameters for component 'Mbefore_wl_monitor' [5]. */
#define mccMbefore_wl_monitor_nL 20
#define mccMbefore_wl_monitor_filename "Mbefore_wl_monitor" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Mbefore_wl_monitor' [5]. */
MCNUM mccMbefore_wl_monitor_xmin;
MCNUM mccMbefore_wl_monitor_xmax;
MCNUM mccMbefore_wl_monitor_ymin;
MCNUM mccMbefore_wl_monitor_ymax;
MCNUM mccMbefore_wl_monitor_xwidth;
MCNUM mccMbefore_wl_monitor_yheight;
MCNUM mccMbefore_wl_monitor_Lmin;
MCNUM mccMbefore_wl_monitor_Lmax;
MCNUM mccMbefore_wl_monitor_restore_xray;

/* Definition parameters for component 'Mbefore_psd_monitor' [6]. */
#define mccMbefore_psd_monitor_nx 100
#define mccMbefore_psd_monitor_ny 10
#define mccMbefore_psd_monitor_nr 0
#define mccMbefore_psd_monitor_filename "Mbefore_psd_monitor" /* declared as a string. May produce warnings at compile */
#define mccMbefore_psd_monitor_restore_xray 1
/* Setting parameters for component 'Mbefore_psd_monitor' [6]. */
MCNUM mccMbefore_psd_monitor_xmin;
MCNUM mccMbefore_psd_monitor_xmax;
MCNUM mccMbefore_psd_monitor_ymin;
MCNUM mccMbefore_psd_monitor_ymax;
MCNUM mccMbefore_psd_monitor_xwidth;
MCNUM mccMbefore_psd_monitor_yheight;
MCNUM mccMbefore_psd_monitor_radius;

/* Definition parameters for component 'mirror2' [11]. */
#define mccmirror2_reflec "" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'mirror2' [11]. */
MCNUM mccmirror2_zdepth;
MCNUM mccmirror2_xwidth;
MCNUM mccmirror2_R0;

/* Definition parameters for component 'Plane_grating' [12]. */
#define mccPlane_grating_blazed mcipblazed
#define mccPlane_grating_display 0
#define mccPlane_grating_zdepth 0.136
#define mccPlane_grating_xwidth 0.15
/* Setting parameters for component 'Plane_grating' [12]. */
MCNUM mccPlane_grating_MCangleVariation;
MCNUM mccPlane_grating_R0;
MCNUM mccPlane_grating_r_rho;
MCNUM mccPlane_grating_b;
MCNUM mccPlane_grating_N_slits;
MCNUM mccPlane_grating_d;
MCNUM mccPlane_grating_blazed_angle;

/* Definition parameters for component 'Close_wl_monitor1' [14]. */
#define mccClose_wl_monitor1_nL 20
#define mccClose_wl_monitor1_filename "Closewl_monitor2" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Close_wl_monitor1' [14]. */
MCNUM mccClose_wl_monitor1_xmin;
MCNUM mccClose_wl_monitor1_xmax;
MCNUM mccClose_wl_monitor1_ymin;
MCNUM mccClose_wl_monitor1_ymax;
MCNUM mccClose_wl_monitor1_xwidth;
MCNUM mccClose_wl_monitor1_yheight;
MCNUM mccClose_wl_monitor1_Lmin;
MCNUM mccClose_wl_monitor1_Lmax;
MCNUM mccClose_wl_monitor1_restore_xray;

/* Definition parameters for component 'close_E_monitor1' [15]. */
#define mccclose_E_monitor1_nE 20
#define mccclose_E_monitor1_filename "CloseE_monitor" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'close_E_monitor1' [15]. */
MCNUM mccclose_E_monitor1_xmin;
MCNUM mccclose_E_monitor1_xmax;
MCNUM mccclose_E_monitor1_ymin;
MCNUM mccclose_E_monitor1_ymax;
MCNUM mccclose_E_monitor1_xwidth;
MCNUM mccclose_E_monitor1_yheight;
MCNUM mccclose_E_monitor1_Emin;
MCNUM mccclose_E_monitor1_Emax;
MCNUM mccclose_E_monitor1_restore_xray;

/* Definition parameters for component 'close_psd_monitor2' [16]. */
#define mccclose_psd_monitor2_nx 50
#define mccclose_psd_monitor2_ny 20
#define mccclose_psd_monitor2_nr 0
#define mccclose_psd_monitor2_filename "Close_psd_monitor" /* declared as a string. May produce warnings at compile */
#define mccclose_psd_monitor2_restore_xray 1
/* Setting parameters for component 'close_psd_monitor2' [16]. */
MCNUM mccclose_psd_monitor2_xmin;
MCNUM mccclose_psd_monitor2_xmax;
MCNUM mccclose_psd_monitor2_ymin;
MCNUM mccclose_psd_monitor2_ymax;
MCNUM mccclose_psd_monitor2_xwidth;
MCNUM mccclose_psd_monitor2_yheight;
MCNUM mccclose_psd_monitor2_radius;

/* Definition parameters for component 'Far_E_monitor2' [17]. */
#define mccFar_E_monitor2_nE 20
#define mccFar_E_monitor2_filename "Far_E_monitor" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Far_E_monitor2' [17]. */
MCNUM mccFar_E_monitor2_xmin;
MCNUM mccFar_E_monitor2_xmax;
MCNUM mccFar_E_monitor2_ymin;
MCNUM mccFar_E_monitor2_ymax;
MCNUM mccFar_E_monitor2_xwidth;
MCNUM mccFar_E_monitor2_yheight;
MCNUM mccFar_E_monitor2_Emin;
MCNUM mccFar_E_monitor2_Emax;
MCNUM mccFar_E_monitor2_restore_xray;

/* Definition parameters for component 'Fardivergence_monitor' [18]. */
#define mccFardivergence_monitor_nh 20
#define mccFardivergence_monitor_nv 20
#define mccFardivergence_monitor_rad 0
/* Setting parameters for component 'Fardivergence_monitor' [18]. */
char mccFardivergence_monitor_filename[16384];
MCNUM mccFardivergence_monitor_xmin;
MCNUM mccFardivergence_monitor_xmax;
MCNUM mccFardivergence_monitor_ymin;
MCNUM mccFardivergence_monitor_ymax;
MCNUM mccFardivergence_monitor_xwidth;
MCNUM mccFardivergence_monitor_yheight;
MCNUM mccFardivergence_monitor_maxdiv_h;
MCNUM mccFardivergence_monitor_maxdiv_v;
MCNUM mccFardivergence_monitor_restore_xray;
MCNUM mccFardivergence_monitor_nx;
MCNUM mccFardivergence_monitor_ny;
MCNUM mccFardivergence_monitor_nz;

/* Definition parameters for component 'Far_wl_monitor2' [19]. */
#define mccFar_wl_monitor2_nL 20
#define mccFar_wl_monitor2_filename "far_wl_monitor" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Far_wl_monitor2' [19]. */
MCNUM mccFar_wl_monitor2_xmin;
MCNUM mccFar_wl_monitor2_xmax;
MCNUM mccFar_wl_monitor2_ymin;
MCNUM mccFar_wl_monitor2_ymax;
MCNUM mccFar_wl_monitor2_xwidth;
MCNUM mccFar_wl_monitor2_yheight;
MCNUM mccFar_wl_monitor2_Lmin;
MCNUM mccFar_wl_monitor2_Lmax;
MCNUM mccFar_wl_monitor2_restore_xray;

/* Definition parameters for component 'Far_psd_monitor2' [20]. */
#define mccFar_psd_monitor2_nx 50
#define mccFar_psd_monitor2_ny 20
#define mccFar_psd_monitor2_nr 0
#define mccFar_psd_monitor2_filename "far_psd_monitor" /* declared as a string. May produce warnings at compile */
#define mccFar_psd_monitor2_restore_xray 1
/* Setting parameters for component 'Far_psd_monitor2' [20]. */
MCNUM mccFar_psd_monitor2_xmin;
MCNUM mccFar_psd_monitor2_xmax;
MCNUM mccFar_psd_monitor2_ymin;
MCNUM mccFar_psd_monitor2_ymax;
MCNUM mccFar_psd_monitor2_xwidth;
MCNUM mccFar_psd_monitor2_yheight;
MCNUM mccFar_psd_monitor2_radius;

/* User component declarations. */

/* User declarations for component 'origin' [1]. */
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define profile mccorigin_profile
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 50 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

  double IntermediateCnts=0;
  time_t StartTime       =0;
  time_t EndTime         =0;
  time_t CurrentTime     =0;
#line 7026 "Monochromator_testing.c"
#undef minutes
#undef flag_save
#undef percent
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef profile
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source_flat' [2]. */
#define mccompcurname  source_flat
#define mccompcurtype  Source_flat
#define mccompcurindex 2
#define spectrum_file mccsource_flat_spectrum_file
#define prms mccsource_flat_prms
#define square mccsource_flat_square
#define radius mccsource_flat_radius
#define yheight mccsource_flat_yheight
#define xwidth mccsource_flat_xwidth
#define xmin mccsource_flat_xmin
#define xmax mccsource_flat_xmax
#define dist mccsource_flat_dist
#define focus_xw mccsource_flat_focus_xw
#define focus_yh mccsource_flat_focus_yh
#define E0 mccsource_flat_E0
#define dE mccsource_flat_dE
#define lambda0 mccsource_flat_lambda0
#define dlambda mccsource_flat_dlambda
#define flux mccsource_flat_flux
#define gauss mccsource_flat_gauss
#define randomphase mccsource_flat_randomphase
#define phase mccsource_flat_phase
#line 63 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
  double pmul, srcArea;
  int square;
  struct {
    double l0,dl;
    double pmul,pint;
    t_Table T;
  } prms;
#line 7069 "Monochromator_testing.c"
#undef phase
#undef randomphase
#undef gauss
#undef flux
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_yh
#undef focus_xw
#undef dist
#undef xmax
#undef xmin
#undef xwidth
#undef yheight
#undef radius
#undef square
#undef prms
#undef spectrum_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'dmu' [3]. */
#define mccompcurname  dmu
#define mccompcurtype  Undulator
#define mccompcurindex 3
#define verbose mccdmu_verbose
#define prms mccdmu_prms
#define alpha mccdmu_alpha
#define MELE mccdmu_MELE
#define E0 mccdmu_E0
#define dE mccdmu_dE
#define phase mccdmu_phase
#define randomphase mccdmu_randomphase
#define Ee mccdmu_Ee
#define dEe mccdmu_dEe
#define Ie mccdmu_Ie
#define tbunch mccdmu_tbunch
#define t0 mccdmu_t0
#define B mccdmu_B
#define K mccdmu_K
#define gap mccdmu_gap
#define Nper mccdmu_Nper
#define lu mccdmu_lu
#define sigey mccdmu_sigey
#define sigex mccdmu_sigex
#define sigepx mccdmu_sigepx
#define sigepy mccdmu_sigepy
#define focus_xw mccdmu_focus_xw
#define focus_yh mccdmu_focus_yh
#define dist mccdmu_dist
#define gauss_t mccdmu_gauss_t
#define quick_integ mccdmu_quick_integ
#define E1st mccdmu_E1st
#line 107 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"
  struct {
    double gamma,gamma2,igamma;
    double s1x,s1y; /*beam's size at dist (convolution of sigex/sigey and igamma)*/
    double length; /*undulator magnetic length*/
    double kc; /*undulator kritical wavenumber*/
    double pmul; /*initial photon weight*/
    gsl_function Bsig,Bpi;
    gsl_integration_workspace *gsl_int_ws;
  } prms;

  /*fine structure constant from CODATA*/
  const double alpha=7.2973525698e-3;
  /*electron mass from CODATA in kg*/
  const double MELE=9.10938291e-31;
  //double besselj[nharm],besselh[nharm];


#line 7143 "Monochromator_testing.c"
#undef E1st
#undef quick_integ
#undef gauss_t
#undef dist
#undef focus_yh
#undef focus_xw
#undef sigepy
#undef sigepx
#undef sigex
#undef sigey
#undef lu
#undef Nper
#undef gap
#undef K
#undef B
#undef t0
#undef tbunch
#undef Ie
#undef dEe
#undef Ee
#undef randomphase
#undef phase
#undef dE
#undef E0
#undef MELE
#undef alpha
#undef prms
#undef verbose
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mbefore_E_monitor' [4]. */
#define mccompcurname  Mbefore_E_monitor
#define mccompcurtype  E_monitor
#define mccompcurindex 4
#define nE mccMbefore_E_monitor_nE
#define filename mccMbefore_E_monitor_filename
#define E_N mccMbefore_E_monitor_E_N
#define E_p mccMbefore_E_monitor_E_p
#define E_p2 mccMbefore_E_monitor_E_p2
#define xmin mccMbefore_E_monitor_xmin
#define xmax mccMbefore_E_monitor_xmax
#define ymin mccMbefore_E_monitor_ymin
#define ymax mccMbefore_E_monitor_ymax
#define xwidth mccMbefore_E_monitor_xwidth
#define yheight mccMbefore_E_monitor_yheight
#define Emin mccMbefore_E_monitor_Emin
#define Emax mccMbefore_E_monitor_Emax
#define restore_xray mccMbefore_E_monitor_restore_xray
#line 60 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
    double E_N[nE];
    double E_p[nE], E_p2[nE];
#line 7197 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mbefore_wl_monitor' [5]. */
#define mccompcurname  Mbefore_wl_monitor
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccMbefore_wl_monitor_nL
#define filename mccMbefore_wl_monitor_filename
#define L_N mccMbefore_wl_monitor_L_N
#define L_p mccMbefore_wl_monitor_L_p
#define L_p2 mccMbefore_wl_monitor_L_p2
#define xmin mccMbefore_wl_monitor_xmin
#define xmax mccMbefore_wl_monitor_xmax
#define ymin mccMbefore_wl_monitor_ymin
#define ymax mccMbefore_wl_monitor_ymax
#define xwidth mccMbefore_wl_monitor_xwidth
#define yheight mccMbefore_wl_monitor_yheight
#define Lmin mccMbefore_wl_monitor_Lmin
#define Lmax mccMbefore_wl_monitor_Lmax
#define restore_xray mccMbefore_wl_monitor_restore_xray
#line 59 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
    double L_N[nL];
    double L_p[nL], L_p2[nL];
#line 7237 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mbefore_psd_monitor' [6]. */
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
#define xmin mccMbefore_psd_monitor_xmin
#define xmax mccMbefore_psd_monitor_xmax
#define ymin mccMbefore_psd_monitor_ymin
#define ymax mccMbefore_psd_monitor_ymax
#define xwidth mccMbefore_psd_monitor_xwidth
#define yheight mccMbefore_psd_monitor_yheight
#define radius mccMbefore_psd_monitor_radius
#line 62 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
    double **PSD_N;
    double **PSD_p;
    double **PSD_p2;
#line 7279 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PG1_arm' [7]. */
#define mccompcurname  PG1_arm
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'M2_rotation_arm1' [8]. */
#define mccompcurname  M2_rotation_arm1
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'M2_rotation_arm2' [9]. */
#define mccompcurname  M2_rotation_arm2
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'M2_rotation_arm3' [10]. */
#define mccompcurname  M2_rotation_arm3
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mirror2' [11]. */
#define mccompcurname  mirror2
#define mccompcurtype  Mirror
#define mccompcurindex 11
#define reflec mccmirror2_reflec
#define reflec_table mccmirror2_reflec_table
#define prms mccmirror2_prms
#define zdepth mccmirror2_zdepth
#define xwidth mccmirror2_xwidth
#define R0 mccmirror2_R0
#line 48 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Mirror.comp"
    struct {
        double e_min,e_max,e_step,theta_min,theta_max,theta_step;
        int use_reflec_table;
    } prms;
    t_Table reflec_table;

#line 7348 "Monochromator_testing.c"
#undef R0
#undef xwidth
#undef zdepth
#undef prms
#undef reflec_table
#undef reflec
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Plane_grating' [12]. */
#define mccompcurname  Plane_grating
#define mccompcurtype  MultiPurpose_grating
#define mccompcurindex 12
#define blazed mccPlane_grating_blazed
#define display mccPlane_grating_display
#define zdepth mccPlane_grating_zdepth
#define xwidth mccPlane_grating_xwidth
#define MCangleVariation mccPlane_grating_MCangleVariation
#define R0 mccPlane_grating_R0
#define r_rho mccPlane_grating_r_rho
#define b mccPlane_grating_b
#define N_slits mccPlane_grating_N_slits
#define d mccPlane_grating_d
#define blazed_angle mccPlane_grating_blazed_angle
#line 83 "MultiPurpose_grating.comp"
#include <complex.h>

    double N_slits, b,d;
    double MCangleVariation,DirectionalSamplingFactor,R0;
    double blazed_angle, r_rho;
#line 7380 "Monochromator_testing.c"
#undef blazed_angle
#undef d
#undef N_slits
#undef b
#undef r_rho
#undef R0
#undef MCangleVariation
#undef xwidth
#undef zdepth
#undef display
#undef blazed
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Monochromator_Monitor_arm' [13]. */
#define mccompcurname  Monochromator_Monitor_arm
#define mccompcurtype  Arm
#define mccompcurindex 13
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Close_wl_monitor1' [14]. */
#define mccompcurname  Close_wl_monitor1
#define mccompcurtype  L_monitor
#define mccompcurindex 14
#define nL mccClose_wl_monitor1_nL
#define filename mccClose_wl_monitor1_filename
#define L_N mccClose_wl_monitor1_L_N
#define L_p mccClose_wl_monitor1_L_p
#define L_p2 mccClose_wl_monitor1_L_p2
#define xmin mccClose_wl_monitor1_xmin
#define xmax mccClose_wl_monitor1_xmax
#define ymin mccClose_wl_monitor1_ymin
#define ymax mccClose_wl_monitor1_ymax
#define xwidth mccClose_wl_monitor1_xwidth
#define yheight mccClose_wl_monitor1_yheight
#define Lmin mccClose_wl_monitor1_Lmin
#define Lmax mccClose_wl_monitor1_Lmax
#define restore_xray mccClose_wl_monitor1_restore_xray
#line 59 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
    double L_N[nL];
    double L_p[nL], L_p2[nL];
#line 7425 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'close_E_monitor1' [15]. */
#define mccompcurname  close_E_monitor1
#define mccompcurtype  E_monitor
#define mccompcurindex 15
#define nE mccclose_E_monitor1_nE
#define filename mccclose_E_monitor1_filename
#define E_N mccclose_E_monitor1_E_N
#define E_p mccclose_E_monitor1_E_p
#define E_p2 mccclose_E_monitor1_E_p2
#define xmin mccclose_E_monitor1_xmin
#define xmax mccclose_E_monitor1_xmax
#define ymin mccclose_E_monitor1_ymin
#define ymax mccclose_E_monitor1_ymax
#define xwidth mccclose_E_monitor1_xwidth
#define yheight mccclose_E_monitor1_yheight
#define Emin mccclose_E_monitor1_Emin
#define Emax mccclose_E_monitor1_Emax
#define restore_xray mccclose_E_monitor1_restore_xray
#line 60 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
    double E_N[nE];
    double E_p[nE], E_p2[nE];
#line 7465 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'close_psd_monitor2' [16]. */
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
#define xmin mccclose_psd_monitor2_xmin
#define xmax mccclose_psd_monitor2_xmax
#define ymin mccclose_psd_monitor2_ymin
#define ymax mccclose_psd_monitor2_ymax
#define xwidth mccclose_psd_monitor2_xwidth
#define yheight mccclose_psd_monitor2_yheight
#define radius mccclose_psd_monitor2_radius
#line 62 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
    double **PSD_N;
    double **PSD_p;
    double **PSD_p2;
#line 7507 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Far_E_monitor2' [17]. */
#define mccompcurname  Far_E_monitor2
#define mccompcurtype  E_monitor
#define mccompcurindex 17
#define nE mccFar_E_monitor2_nE
#define filename mccFar_E_monitor2_filename
#define E_N mccFar_E_monitor2_E_N
#define E_p mccFar_E_monitor2_E_p
#define E_p2 mccFar_E_monitor2_E_p2
#define xmin mccFar_E_monitor2_xmin
#define xmax mccFar_E_monitor2_xmax
#define ymin mccFar_E_monitor2_ymin
#define ymax mccFar_E_monitor2_ymax
#define xwidth mccFar_E_monitor2_xwidth
#define yheight mccFar_E_monitor2_yheight
#define Emin mccFar_E_monitor2_Emin
#define Emax mccFar_E_monitor2_Emax
#define restore_xray mccFar_E_monitor2_restore_xray
#line 60 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
    double E_N[nE];
    double E_p[nE], E_p2[nE];
#line 7548 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fardivergence_monitor' [18]. */
#define mccompcurname  Fardivergence_monitor
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 18
#define nh mccFardivergence_monitor_nh
#define nv mccFardivergence_monitor_nv
#define rad mccFardivergence_monitor_rad
#define Div_N mccFardivergence_monitor_Div_N
#define Div_p mccFardivergence_monitor_Div_p
#define Div_p2 mccFardivergence_monitor_Div_p2
#define filename mccFardivergence_monitor_filename
#define xmin mccFardivergence_monitor_xmin
#define xmax mccFardivergence_monitor_xmax
#define ymin mccFardivergence_monitor_ymin
#define ymax mccFardivergence_monitor_ymax
#define xwidth mccFardivergence_monitor_xwidth
#define yheight mccFardivergence_monitor_yheight
#define maxdiv_h mccFardivergence_monitor_maxdiv_h
#define maxdiv_v mccFardivergence_monitor_maxdiv_v
#define restore_xray mccFardivergence_monitor_restore_xray
#define nx mccFardivergence_monitor_nx
#define ny mccFardivergence_monitor_ny
#define nz mccFardivergence_monitor_nz
#line 64 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/Divergence_monitor.comp"
    double Div_N[nh][nv];
    double Div_p[nh][nv];
    double Div_p2[nh][nv];
#line 7594 "Monochromator_testing.c"
#undef nz
#undef ny
#undef nx
#undef restore_xray
#undef maxdiv_v
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef rad
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Far_wl_monitor2' [19]. */
#define mccompcurname  Far_wl_monitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mccFar_wl_monitor2_nL
#define filename mccFar_wl_monitor2_filename
#define L_N mccFar_wl_monitor2_L_N
#define L_p mccFar_wl_monitor2_L_p
#define L_p2 mccFar_wl_monitor2_L_p2
#define xmin mccFar_wl_monitor2_xmin
#define xmax mccFar_wl_monitor2_xmax
#define ymin mccFar_wl_monitor2_ymin
#define ymax mccFar_wl_monitor2_ymax
#define xwidth mccFar_wl_monitor2_xwidth
#define yheight mccFar_wl_monitor2_yheight
#define Lmin mccFar_wl_monitor2_Lmin
#define Lmax mccFar_wl_monitor2_Lmax
#define restore_xray mccFar_wl_monitor2_restore_xray
#line 59 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
    double L_N[nL];
    double L_p[nL], L_p2[nL];
#line 7639 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Far_psd_monitor2' [20]. */
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
#define xmin mccFar_psd_monitor2_xmin
#define xmax mccFar_psd_monitor2_xmax
#define ymin mccFar_psd_monitor2_ymin
#define ymax mccFar_psd_monitor2_ymax
#define xwidth mccFar_psd_monitor2_xwidth
#define yheight mccFar_psd_monitor2_yheight
#define radius mccFar_psd_monitor2_radius
#line 62 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
    double **PSD_N;
    double **PSD_p;
    double **PSD_p2;
#line 7681 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaorigin, mcposrorigin;
Rotation mcrotaorigin, mcrotrorigin;
Coords mcposasource_flat, mcposrsource_flat;
Rotation mcrotasource_flat, mcrotrsource_flat;
Coords mcposadmu, mcposrdmu;
Rotation mcrotadmu, mcrotrdmu;
Coords mcposaMbefore_E_monitor, mcposrMbefore_E_monitor;
Rotation mcrotaMbefore_E_monitor, mcrotrMbefore_E_monitor;
Coords mcposaMbefore_wl_monitor, mcposrMbefore_wl_monitor;
Rotation mcrotaMbefore_wl_monitor, mcrotrMbefore_wl_monitor;
Coords mcposaMbefore_psd_monitor, mcposrMbefore_psd_monitor;
Rotation mcrotaMbefore_psd_monitor, mcrotrMbefore_psd_monitor;
Coords mcposaPG1_arm, mcposrPG1_arm;
Rotation mcrotaPG1_arm, mcrotrPG1_arm;
Coords mcposaM2_rotation_arm1, mcposrM2_rotation_arm1;
Rotation mcrotaM2_rotation_arm1, mcrotrM2_rotation_arm1;
Coords mcposaM2_rotation_arm2, mcposrM2_rotation_arm2;
Rotation mcrotaM2_rotation_arm2, mcrotrM2_rotation_arm2;
Coords mcposaM2_rotation_arm3, mcposrM2_rotation_arm3;
Rotation mcrotaM2_rotation_arm3, mcrotrM2_rotation_arm3;
Coords mcposamirror2, mcposrmirror2;
Rotation mcrotamirror2, mcrotrmirror2;
Coords mcposaPlane_grating, mcposrPlane_grating;
Rotation mcrotaPlane_grating, mcrotrPlane_grating;
Coords mcposaMonochromator_Monitor_arm, mcposrMonochromator_Monitor_arm;
Rotation mcrotaMonochromator_Monitor_arm, mcrotrMonochromator_Monitor_arm;
Coords mcposaClose_wl_monitor1, mcposrClose_wl_monitor1;
Rotation mcrotaClose_wl_monitor1, mcrotrClose_wl_monitor1;
Coords mcposaclose_E_monitor1, mcposrclose_E_monitor1;
Rotation mcrotaclose_E_monitor1, mcrotrclose_E_monitor1;
Coords mcposaclose_psd_monitor2, mcposrclose_psd_monitor2;
Rotation mcrotaclose_psd_monitor2, mcrotrclose_psd_monitor2;
Coords mcposaFar_E_monitor2, mcposrFar_E_monitor2;
Rotation mcrotaFar_E_monitor2, mcrotrFar_E_monitor2;
Coords mcposaFardivergence_monitor, mcposrFardivergence_monitor;
Rotation mcrotaFardivergence_monitor, mcrotrFardivergence_monitor;
Coords mcposaFar_wl_monitor2, mcposrFar_wl_monitor2;
Rotation mcrotaFar_wl_monitor2, mcrotrFar_wl_monitor2;
Coords mcposaFar_psd_monitor2, mcposrFar_psd_monitor2;
Rotation mcrotaFar_psd_monitor2, mcrotrFar_psd_monitor2;

MCNUM mcnx, mcny, mcnz, mcnkx, mcnky, mcnkz, mcnphi, mcnt, mcnEx, mcnEy, mcnEz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Monochroamtor_testing
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaMonochroamtor_testing coords_set(0,0,0)
#define r_rho mcipr_rho
#define m mcipm
#define E0 mcipE0
#define dE mcipdE
#define blazed mcipblazed
#define blazed_angle mcipblazed_angle
#define cff mcipcff
#define Wanted_energy mcipWanted_energy
#define SourceChoice mcipSourceChoice
#line 48 "Monochromator_testing.instr"
{

/*******************************************************
Finding wanted energy if none is given
*******************************************************/
if (!Wanted_energy){
        /* If no wanted energy is given, it is assumed the wanted energy is E0. */
        printf("Warning: No wanted energy is given. Default is E0=%f keV \n",E0);
        Wanted_energy = E0;
}
/*******************************************************
Distances used in the monochromator:
*******************************************************/
A = 0.064819;
B = 0.000125;
C = 0.043821;
D = 0.020;
E = 0.620;
F = 0.042;
G = 0.140;
X = 0.072746;
/*******************************************************
Calculating the angle for the monochromator. For motivation, see Urpelainen, Samuli 2014.
*******************************************************/

if (cff && blazed){
            double cPGM_wl,cPGM_rho, cPGM_a,cPGM_b,cPGM_c,cPGM_beta,cPGM_alpha;
            cPGM_wl = (12.398/Wanted_energy)*pow(10,-10);
            cPGM_rho = 1/(r_rho*1000);
            cPGM_a = (1-pow(cff,2));
            cPGM_b = 2*pow(cff,2)*(m*cPGM_wl/cPGM_rho);
            cPGM_c = pow(cff,2)-1-pow(cff,2)*((pow(m,2)*pow(cPGM_wl,2))/pow(cPGM_rho,2)); 
            cPGM_beta = asin((-cPGM_b+sqrt(pow(cPGM_b,2)-4*cPGM_a*cPGM_c))/(2*cPGM_a));      
            cPGM_alpha =acos(cos(cPGM_beta)/cff);
            cPGM_beta = cPGM_beta*RAD2DEG;
            cPGM_alpha = cPGM_alpha*RAD2DEG;
            angle_grating = (cPGM_beta+90);   
            mirror2_angle = ((90+angle_grating-cPGM_alpha))/2;
            printf("a=%f, b=%f, c=%f, beta=%f, alpha=%f,PG=%f, M2%f \n",cPGM_a,cPGM_b,cPGM_c,cPGM_beta,cPGM_alpha,angle_grating,mirror2_angle);
} else {
            printf("Error. The cPGM testing instrument is used for blazed gratings only.");
            exit(-1);
}
printf("Input specs: \n     Wanted energy: %f keV. \n     Grating order: %f. \n     Wanted wavelength: %f AA. \n",Wanted_energy,m,12.398/Wanted_energy);
printf("\nMonochromator specs: \n     Angle of pre-mirror=%f deg. \n     Angle of grating=%f deg. \n",mirror2_angle,angle_grating);
//Here, in testing the monochromator the MC angle need to be "big enough". 
MCAngleVariation=angle_grating*2; 

Emin=(Wanted_energy-dE)/4;
Emax=(Wanted_energy+dE)*4;
Lmax=12.398/Emax;
Lmin=12.398/Emin;

if (Wanted_energy<0.025){
printf("Warning. Energy below recommended value.(%f)\n",Wanted_energy);
} else if (Wanted_energy>1) {
printf("Warning. Energy above recommended value.(%f)\n",Wanted_energy);
}
if (mirror2_angle>21){
printf("Warning. Angle of pre-mirror is higher than allowed value on Bloch.(%f)\n",mirror2_angle);
}
if (angle_grating>30){
printf("Warning. Angle of grating is higher than allowed value on Bloch.(%f)\n",angle_grating);
}

if (blazed_angle==2 && r_rho==800){
printf("Grating PG1 is used.\n") ;
}
else if(blazed_angle==6 && r_rho==2000){
printf("Grating PG2 is used.\n"); 
}
else {
printf("Grating not in use on Bloch");
}

E1st=1.0018*Wanted_energy/7;
undK=5.6;
Nper=187;
}
#line 7840 "Monochromator_testing.c"
#undef SourceChoice
#undef Wanted_energy
#undef cff
#undef blazed_angle
#undef blazed
#undef dE
#undef E0
#undef m
#undef r_rho
#undef mcposaMonochroamtor_testing
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component origin. */
  /* Setting parameters for component origin. */
  SIG_MESSAGE("origin (Init:SetPar)");
#line 44 "Monochromator_testing.instr"
  mccorigin_percent = 10;
#line 44 "Monochromator_testing.instr"
  mccorigin_flag_save = 0;
#line 44 "Monochromator_testing.instr"
  mccorigin_minutes = 0;
#line 7873 "Monochromator_testing.c"

  SIG_MESSAGE("origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaorigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7880 "Monochromator_testing.c"
  rot_copy(mcrotrorigin, mcrotaorigin);
  mcposaorigin = coords_set(
#line 130 "Monochromator_testing.instr"
    0,
#line 130 "Monochromator_testing.instr"
    0,
#line 130 "Monochromator_testing.instr"
    0);
#line 7889 "Monochromator_testing.c"
  mctc1 = coords_neg(mcposaorigin);
  mcposrorigin = rot_apply(mcrotaorigin, mctc1);
  mcDEBUG_COMPONENT("origin", mcposaorigin, mcrotaorigin)
  mccomp_posa[1] = mcposaorigin;
  mccomp_posr[1] = mcposrorigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source_flat. */
  /* Setting parameters for component source_flat. */
  SIG_MESSAGE("source_flat (Init:SetPar)");
#line 50 "Monochromator_testing.instr"
  mccsource_flat_radius = 0;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_yheight = 0.001;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_xwidth = 0.001;
#line 50 "Monochromator_testing.instr"
  mccsource_flat_xmin = 0;
#line 50 "Monochromator_testing.instr"
  mccsource_flat_xmax = 0;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_dist = 2;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_focus_xw = 0.0001;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_focus_yh = 0.0001;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_E0 = mcipE0;
#line 133 "Monochromator_testing.instr"
  mccsource_flat_dE = mcipdE;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_lambda0 = 0;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_dlambda = 0;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_flux = 0;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_gauss = 0;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_randomphase = 1;
#line 52 "Monochromator_testing.instr"
  mccsource_flat_phase = 0;
#line 7932 "Monochromator_testing.c"

  SIG_MESSAGE("source_flat (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7939 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaorigin, mcrotasource_flat);
  rot_transpose(mcrotaorigin, mctr1);
  rot_mul(mcrotasource_flat, mctr1, mcrotrsource_flat);
  mctc1 = coords_set(
#line 135 "Monochromator_testing.instr"
    0,
#line 135 "Monochromator_testing.instr"
    0,
#line 135 "Monochromator_testing.instr"
    0);
#line 7950 "Monochromator_testing.c"
  rot_transpose(mcrotaorigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource_flat = coords_add(mcposaorigin, mctc2);
  mctc1 = coords_sub(mcposaorigin, mcposasource_flat);
  mcposrsource_flat = rot_apply(mcrotasource_flat, mctc1);
  mcDEBUG_COMPONENT("source_flat", mcposasource_flat, mcrotasource_flat)
  mccomp_posa[2] = mcposasource_flat;
  mccomp_posr[2] = mcposrsource_flat;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component dmu. */
  /* Setting parameters for component dmu. */
  SIG_MESSAGE("dmu (Init:SetPar)");
#line 137 "Monochromator_testing.instr"
  mccdmu_E0 = mcipE0;
#line 137 "Monochromator_testing.instr"
  mccdmu_dE = mcipdE;
#line 58 "Monochromator_testing.instr"
  mccdmu_phase = 0;
#line 58 "Monochromator_testing.instr"
  mccdmu_randomphase = 1;
#line 137 "Monochromator_testing.instr"
  mccdmu_Ee = 1.5;
#line 137 "Monochromator_testing.instr"
  mccdmu_dEe = ( ( 6e-9 ) * ( 60e-12 ) ) / 1.5;
#line 137 "Monochromator_testing.instr"
  mccdmu_Ie = 0.5;
#line 137 "Monochromator_testing.instr"
  mccdmu_tbunch = 43;
#line 58 "Monochromator_testing.instr"
  mccdmu_t0 = 0;
#line 58 "Monochromator_testing.instr"
  mccdmu_B = 0;
#line 137 "Monochromator_testing.instr"
  mccdmu_K = undK;
#line 137 "Monochromator_testing.instr"
  mccdmu_gap = 14e-3;
#line 137 "Monochromator_testing.instr"
  mccdmu_Nper = Nper;
#line 138 "Monochromator_testing.instr"
  mccdmu_lu = 84e-3;
#line 138 "Monochromator_testing.instr"
  mccdmu_sigey = 1.3e-5;
#line 138 "Monochromator_testing.instr"
  mccdmu_sigex = 185e-5;
#line 138 "Monochromator_testing.instr"
  mccdmu_sigepx = 32e-6;
#line 138 "Monochromator_testing.instr"
  mccdmu_sigepy = 4.6e-6;
#line 138 "Monochromator_testing.instr"
  mccdmu_focus_xw = 1.1e-4;
#line 138 "Monochromator_testing.instr"
  mccdmu_focus_yh = 1.1e-4;
#line 138 "Monochromator_testing.instr"
  mccdmu_dist = 2;
#line 59 "Monochromator_testing.instr"
  mccdmu_gauss_t = 0;
#line 59 "Monochromator_testing.instr"
  mccdmu_quick_integ = 0;
#line 138 "Monochromator_testing.instr"
  mccdmu_E1st = E1st;
#line 8012 "Monochromator_testing.c"

  SIG_MESSAGE("dmu (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8019 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaorigin, mcrotadmu);
  rot_transpose(mcrotasource_flat, mctr1);
  rot_mul(mcrotadmu, mctr1, mcrotrdmu);
  mctc1 = coords_set(
#line 140 "Monochromator_testing.instr"
    0,
#line 140 "Monochromator_testing.instr"
    0,
#line 140 "Monochromator_testing.instr"
    0);
#line 8030 "Monochromator_testing.c"
  rot_transpose(mcrotaorigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadmu = coords_add(mcposaorigin, mctc2);
  mctc1 = coords_sub(mcposasource_flat, mcposadmu);
  mcposrdmu = rot_apply(mcrotadmu, mctc1);
  mcDEBUG_COMPONENT("dmu", mcposadmu, mcrotadmu)
  mccomp_posa[3] = mcposadmu;
  mccomp_posr[3] = mcposrdmu;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component Mbefore_E_monitor. */
  /* Setting parameters for component Mbefore_E_monitor. */
  SIG_MESSAGE("Mbefore_E_monitor (Init:SetPar)");
#line 54 "Monochromator_testing.instr"
  mccMbefore_E_monitor_xmin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccMbefore_E_monitor_xmax = 0.05;
#line 54 "Monochromator_testing.instr"
  mccMbefore_E_monitor_ymin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccMbefore_E_monitor_ymax = 0.05;
#line 145 "Monochromator_testing.instr"
  mccMbefore_E_monitor_xwidth = 0.06;
#line 145 "Monochromator_testing.instr"
  mccMbefore_E_monitor_yheight = 0.1;
#line 145 "Monochromator_testing.instr"
  mccMbefore_E_monitor_Emin = Emin;
#line 145 "Monochromator_testing.instr"
  mccMbefore_E_monitor_Emax = Emax;
#line 145 "Monochromator_testing.instr"
  mccMbefore_E_monitor_restore_xray = 1;
#line 8062 "Monochromator_testing.c"

  SIG_MESSAGE("Mbefore_E_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8069 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotadmu, mcrotaMbefore_E_monitor);
  rot_transpose(mcrotadmu, mctr1);
  rot_mul(mcrotaMbefore_E_monitor, mctr1, mcrotrMbefore_E_monitor);
  mctc1 = coords_set(
#line 146 "Monochromator_testing.instr"
    0,
#line 146 "Monochromator_testing.instr"
    0,
#line 146 "Monochromator_testing.instr"
    0.5);
#line 8080 "Monochromator_testing.c"
  rot_transpose(mcrotadmu, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMbefore_E_monitor = coords_add(mcposadmu, mctc2);
  mctc1 = coords_sub(mcposadmu, mcposaMbefore_E_monitor);
  mcposrMbefore_E_monitor = rot_apply(mcrotaMbefore_E_monitor, mctc1);
  mcDEBUG_COMPONENT("Mbefore_E_monitor", mcposaMbefore_E_monitor, mcrotaMbefore_E_monitor)
  mccomp_posa[4] = mcposaMbefore_E_monitor;
  mccomp_posr[4] = mcposrMbefore_E_monitor;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component Mbefore_wl_monitor. */
  /* Setting parameters for component Mbefore_wl_monitor. */
  SIG_MESSAGE("Mbefore_wl_monitor (Init:SetPar)");
#line 53 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_xmin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_xmax = 0.05;
#line 53 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_ymin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_ymax = 0.05;
#line 148 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_xwidth = 0.06;
#line 148 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_yheight = 0.1;
#line 148 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_Lmin = Lmin;
#line 148 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_Lmax = Lmax;
#line 148 "Monochromator_testing.instr"
  mccMbefore_wl_monitor_restore_xray = 1;
#line 8112 "Monochromator_testing.c"

  SIG_MESSAGE("Mbefore_wl_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8119 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMbefore_E_monitor, mcrotaMbefore_wl_monitor);
  rot_transpose(mcrotaMbefore_E_monitor, mctr1);
  rot_mul(mcrotaMbefore_wl_monitor, mctr1, mcrotrMbefore_wl_monitor);
  mctc1 = coords_set(
#line 149 "Monochromator_testing.instr"
    0,
#line 149 "Monochromator_testing.instr"
    0,
#line 149 "Monochromator_testing.instr"
    0);
#line 8130 "Monochromator_testing.c"
  rot_transpose(mcrotaMbefore_E_monitor, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMbefore_wl_monitor = coords_add(mcposaMbefore_E_monitor, mctc2);
  mctc1 = coords_sub(mcposaMbefore_E_monitor, mcposaMbefore_wl_monitor);
  mcposrMbefore_wl_monitor = rot_apply(mcrotaMbefore_wl_monitor, mctc1);
  mcDEBUG_COMPONENT("Mbefore_wl_monitor", mcposaMbefore_wl_monitor, mcrotaMbefore_wl_monitor)
  mccomp_posa[5] = mcposaMbefore_wl_monitor;
  mccomp_posr[5] = mcposrMbefore_wl_monitor;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Mbefore_psd_monitor. */
  /* Setting parameters for component Mbefore_psd_monitor. */
  SIG_MESSAGE("Mbefore_psd_monitor (Init:SetPar)");
#line 56 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_xmin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_xmax = 0.05;
#line 56 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_ymin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_ymax = 0.05;
#line 151 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_xwidth = 0.1;
#line 151 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_yheight = 0.1;
#line 56 "Monochromator_testing.instr"
  mccMbefore_psd_monitor_radius = 0;
#line 8158 "Monochromator_testing.c"

  SIG_MESSAGE("Mbefore_psd_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8165 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMbefore_wl_monitor, mcrotaMbefore_psd_monitor);
  rot_transpose(mcrotaMbefore_wl_monitor, mctr1);
  rot_mul(mcrotaMbefore_psd_monitor, mctr1, mcrotrMbefore_psd_monitor);
  mctc1 = coords_set(
#line 152 "Monochromator_testing.instr"
    0,
#line 152 "Monochromator_testing.instr"
    0,
#line 152 "Monochromator_testing.instr"
    0);
#line 8176 "Monochromator_testing.c"
  rot_transpose(mcrotaMbefore_wl_monitor, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMbefore_psd_monitor = coords_add(mcposaMbefore_wl_monitor, mctc2);
  mctc1 = coords_sub(mcposaMbefore_wl_monitor, mcposaMbefore_psd_monitor);
  mcposrMbefore_psd_monitor = rot_apply(mcrotaMbefore_psd_monitor, mctc1);
  mcDEBUG_COMPONENT("Mbefore_psd_monitor", mcposaMbefore_psd_monitor, mcrotaMbefore_psd_monitor)
  mccomp_posa[6] = mcposaMbefore_psd_monitor;
  mccomp_posr[6] = mcposrMbefore_psd_monitor;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component PG1_arm. */
  /* Setting parameters for component PG1_arm. */
  SIG_MESSAGE("PG1_arm (Init:SetPar)");

  SIG_MESSAGE("PG1_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 158 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 158 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 158 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8199 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaorigin, mcrotaPG1_arm);
  rot_transpose(mcrotaMbefore_psd_monitor, mctr1);
  rot_mul(mcrotaPG1_arm, mctr1, mcrotrPG1_arm);
  mctc1 = coords_set(
#line 157 "Monochromator_testing.instr"
    0,
#line 157 "Monochromator_testing.instr"
    F,
#line 157 "Monochromator_testing.instr"
    2);
#line 8210 "Monochromator_testing.c"
  rot_transpose(mcrotaorigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPG1_arm = coords_add(mcposaorigin, mctc2);
  mctc1 = coords_sub(mcposaMbefore_psd_monitor, mcposaPG1_arm);
  mcposrPG1_arm = rot_apply(mcrotaPG1_arm, mctc1);
  mcDEBUG_COMPONENT("PG1_arm", mcposaPG1_arm, mcrotaPG1_arm)
  mccomp_posa[7] = mcposaPG1_arm;
  mccomp_posr[7] = mcposrPG1_arm;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component M2_rotation_arm1. */
  /* Setting parameters for component M2_rotation_arm1. */
  SIG_MESSAGE("M2_rotation_arm1 (Init:SetPar)");

  SIG_MESSAGE("M2_rotation_arm1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 162 "Monochromator_testing.instr"
    (- mirror2_angle)*DEG2RAD,
#line 162 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 162 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8233 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaPG1_arm, mcrotaM2_rotation_arm1);
  rot_transpose(mcrotaPG1_arm, mctr1);
  rot_mul(mcrotaM2_rotation_arm1, mctr1, mcrotrM2_rotation_arm1);
  mctc1 = coords_set(
#line 161 "Monochromator_testing.instr"
    0,
#line 161 "Monochromator_testing.instr"
    A - F,
#line 161 "Monochromator_testing.instr"
    B);
#line 8244 "Monochromator_testing.c"
  rot_transpose(mcrotaPG1_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaM2_rotation_arm1 = coords_add(mcposaPG1_arm, mctc2);
  mctc1 = coords_sub(mcposaPG1_arm, mcposaM2_rotation_arm1);
  mcposrM2_rotation_arm1 = rot_apply(mcrotaM2_rotation_arm1, mctc1);
  mcDEBUG_COMPONENT("M2_rotation_arm1", mcposaM2_rotation_arm1, mcrotaM2_rotation_arm1)
  mccomp_posa[8] = mcposaM2_rotation_arm1;
  mccomp_posr[8] = mcposrM2_rotation_arm1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component M2_rotation_arm2. */
  /* Setting parameters for component M2_rotation_arm2. */
  SIG_MESSAGE("M2_rotation_arm2 (Init:SetPar)");

  SIG_MESSAGE("M2_rotation_arm2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 166 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 166 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 166 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8267 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaM2_rotation_arm1, mcrotaM2_rotation_arm2);
  rot_transpose(mcrotaM2_rotation_arm1, mctr1);
  rot_mul(mcrotaM2_rotation_arm2, mctr1, mcrotrM2_rotation_arm2);
  mctc1 = coords_set(
#line 165 "Monochromator_testing.instr"
    0,
#line 165 "Monochromator_testing.instr"
    0,
#line 165 "Monochromator_testing.instr"
    - D - ( E / 2 ));
#line 8278 "Monochromator_testing.c"
  rot_transpose(mcrotaM2_rotation_arm1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaM2_rotation_arm2 = coords_add(mcposaM2_rotation_arm1, mctc2);
  mctc1 = coords_sub(mcposaM2_rotation_arm1, mcposaM2_rotation_arm2);
  mcposrM2_rotation_arm2 = rot_apply(mcrotaM2_rotation_arm2, mctc1);
  mcDEBUG_COMPONENT("M2_rotation_arm2", mcposaM2_rotation_arm2, mcrotaM2_rotation_arm2)
  mccomp_posa[9] = mcposaM2_rotation_arm2;
  mccomp_posr[9] = mcposrM2_rotation_arm2;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component M2_rotation_arm3. */
  /* Setting parameters for component M2_rotation_arm3. */
  SIG_MESSAGE("M2_rotation_arm3 (Init:SetPar)");

  SIG_MESSAGE("M2_rotation_arm3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 170 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 170 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 170 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8301 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaM2_rotation_arm2, mcrotaM2_rotation_arm3);
  rot_transpose(mcrotaM2_rotation_arm2, mctr1);
  rot_mul(mcrotaM2_rotation_arm3, mctr1, mcrotrM2_rotation_arm3);
  mctc1 = coords_set(
#line 169 "Monochromator_testing.instr"
    0,
#line 169 "Monochromator_testing.instr"
    - C,
#line 169 "Monochromator_testing.instr"
    0);
#line 8312 "Monochromator_testing.c"
  rot_transpose(mcrotaM2_rotation_arm2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaM2_rotation_arm3 = coords_add(mcposaM2_rotation_arm2, mctc2);
  mctc1 = coords_sub(mcposaM2_rotation_arm2, mcposaM2_rotation_arm3);
  mcposrM2_rotation_arm3 = rot_apply(mcrotaM2_rotation_arm3, mctc1);
  mcDEBUG_COMPONENT("M2_rotation_arm3", mcposaM2_rotation_arm3, mcrotaM2_rotation_arm3)
  mccomp_posa[10] = mcposaM2_rotation_arm3;
  mccomp_posr[10] = mcposrM2_rotation_arm3;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component mirror2. */
  /* Setting parameters for component mirror2. */
  SIG_MESSAGE("mirror2 (Init:SetPar)");
#line 173 "Monochromator_testing.instr"
  mccmirror2_zdepth = 0.57;
#line 173 "Monochromator_testing.instr"
  mccmirror2_xwidth = 0.15;
#line 173 "Monochromator_testing.instr"
  mccmirror2_R0 = 1;
#line 8332 "Monochromator_testing.c"

  SIG_MESSAGE("mirror2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 175 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 175 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 175 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8342 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaM2_rotation_arm3, mcrotamirror2);
  rot_transpose(mcrotaM2_rotation_arm3, mctr1);
  rot_mul(mcrotamirror2, mctr1, mcrotrmirror2);
  mctc1 = coords_set(
#line 174 "Monochromator_testing.instr"
    0,
#line 174 "Monochromator_testing.instr"
    0,
#line 174 "Monochromator_testing.instr"
    0);
#line 8353 "Monochromator_testing.c"
  rot_transpose(mcrotaM2_rotation_arm3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamirror2 = coords_add(mcposaM2_rotation_arm3, mctc2);
  mctc1 = coords_sub(mcposaM2_rotation_arm3, mcposamirror2);
  mcposrmirror2 = rot_apply(mcrotamirror2, mctc1);
  mcDEBUG_COMPONENT("mirror2", mcposamirror2, mcrotamirror2)
  mccomp_posa[11] = mcposamirror2;
  mccomp_posr[11] = mcposrmirror2;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component Plane_grating. */
  /* Setting parameters for component Plane_grating. */
  SIG_MESSAGE("Plane_grating (Init:SetPar)");
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_MCangleVariation = MCAngleVariation;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_R0 = 1;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_r_rho = mcipr_rho;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_b = 0;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_N_slits = 0;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_d = 0;
#line 178 "Monochromator_testing.instr"
  mccPlane_grating_blazed_angle = mcipblazed_angle;
#line 8381 "Monochromator_testing.c"

  SIG_MESSAGE("Plane_grating (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 180 "Monochromator_testing.instr"
    (- angle_grating)*DEG2RAD,
#line 180 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 180 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8391 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaPG1_arm, mcrotaPlane_grating);
  rot_transpose(mcrotamirror2, mctr1);
  rot_mul(mcrotaPlane_grating, mctr1, mcrotrPlane_grating);
  mctc1 = coords_set(
#line 179 "Monochromator_testing.instr"
    0,
#line 179 "Monochromator_testing.instr"
    0,
#line 179 "Monochromator_testing.instr"
    0);
#line 8402 "Monochromator_testing.c"
  rot_transpose(mcrotaPG1_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPlane_grating = coords_add(mcposaPG1_arm, mctc2);
  mctc1 = coords_sub(mcposamirror2, mcposaPlane_grating);
  mcposrPlane_grating = rot_apply(mcrotaPlane_grating, mctc1);
  mcDEBUG_COMPONENT("Plane_grating", mcposaPlane_grating, mcrotaPlane_grating)
  mccomp_posa[12] = mcposaPlane_grating;
  mccomp_posr[12] = mcposrPlane_grating;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component Monochromator_Monitor_arm. */
  /* Setting parameters for component Monochromator_Monitor_arm. */
  SIG_MESSAGE("Monochromator_Monitor_arm (Init:SetPar)");

  SIG_MESSAGE("Monochromator_Monitor_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 186 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 186 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 186 "Monochromator_testing.instr"
    (0)*DEG2RAD);
#line 8425 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaPG1_arm, mcrotaMonochromator_Monitor_arm);
  rot_transpose(mcrotaPlane_grating, mctr1);
  rot_mul(mcrotaMonochromator_Monitor_arm, mctr1, mcrotrMonochromator_Monitor_arm);
  mctc1 = coords_set(
#line 185 "Monochromator_testing.instr"
    0,
#line 185 "Monochromator_testing.instr"
    0,
#line 185 "Monochromator_testing.instr"
    0);
#line 8436 "Monochromator_testing.c"
  rot_transpose(mcrotaPG1_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMonochromator_Monitor_arm = coords_add(mcposaPG1_arm, mctc2);
  mctc1 = coords_sub(mcposaPlane_grating, mcposaMonochromator_Monitor_arm);
  mcposrMonochromator_Monitor_arm = rot_apply(mcrotaMonochromator_Monitor_arm, mctc1);
  mcDEBUG_COMPONENT("Monochromator_Monitor_arm", mcposaMonochromator_Monitor_arm, mcrotaMonochromator_Monitor_arm)
  mccomp_posa[13] = mcposaMonochromator_Monitor_arm;
  mccomp_posr[13] = mcposrMonochromator_Monitor_arm;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component Close_wl_monitor1. */
  /* Setting parameters for component Close_wl_monitor1. */
  SIG_MESSAGE("Close_wl_monitor1 (Init:SetPar)");
#line 53 "Monochromator_testing.instr"
  mccClose_wl_monitor1_xmin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccClose_wl_monitor1_xmax = 0.05;
#line 53 "Monochromator_testing.instr"
  mccClose_wl_monitor1_ymin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccClose_wl_monitor1_ymax = 0.05;
#line 189 "Monochromator_testing.instr"
  mccClose_wl_monitor1_xwidth = 0.8;
#line 189 "Monochromator_testing.instr"
  mccClose_wl_monitor1_yheight = 0.8;
#line 189 "Monochromator_testing.instr"
  mccClose_wl_monitor1_Lmin = Lmin;
#line 189 "Monochromator_testing.instr"
  mccClose_wl_monitor1_Lmax = Lmax;
#line 189 "Monochromator_testing.instr"
  mccClose_wl_monitor1_restore_xray = 1;
#line 8468 "Monochromator_testing.c"

  SIG_MESSAGE("Close_wl_monitor1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 191 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 191 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 191 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8478 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaClose_wl_monitor1);
  rot_transpose(mcrotaMonochromator_Monitor_arm, mctr1);
  rot_mul(mcrotaClose_wl_monitor1, mctr1, mcrotrClose_wl_monitor1);
  mctc1 = coords_set(
#line 190 "Monochromator_testing.instr"
    0,
#line 190 "Monochromator_testing.instr"
    0,
#line 190 "Monochromator_testing.instr"
    0.4);
#line 8489 "Monochromator_testing.c"
  rot_transpose(mcrotaMonochromator_Monitor_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaClose_wl_monitor1 = coords_add(mcposaMonochromator_Monitor_arm, mctc2);
  mctc1 = coords_sub(mcposaMonochromator_Monitor_arm, mcposaClose_wl_monitor1);
  mcposrClose_wl_monitor1 = rot_apply(mcrotaClose_wl_monitor1, mctc1);
  mcDEBUG_COMPONENT("Close_wl_monitor1", mcposaClose_wl_monitor1, mcrotaClose_wl_monitor1)
  mccomp_posa[14] = mcposaClose_wl_monitor1;
  mccomp_posr[14] = mcposrClose_wl_monitor1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component close_E_monitor1. */
  /* Setting parameters for component close_E_monitor1. */
  SIG_MESSAGE("close_E_monitor1 (Init:SetPar)");
#line 54 "Monochromator_testing.instr"
  mccclose_E_monitor1_xmin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccclose_E_monitor1_xmax = 0.05;
#line 54 "Monochromator_testing.instr"
  mccclose_E_monitor1_ymin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccclose_E_monitor1_ymax = 0.05;
#line 193 "Monochromator_testing.instr"
  mccclose_E_monitor1_xwidth = 0.8;
#line 193 "Monochromator_testing.instr"
  mccclose_E_monitor1_yheight = 0.8;
#line 193 "Monochromator_testing.instr"
  mccclose_E_monitor1_Emin = Emin;
#line 193 "Monochromator_testing.instr"
  mccclose_E_monitor1_Emax = Emax;
#line 193 "Monochromator_testing.instr"
  mccclose_E_monitor1_restore_xray = 1;
#line 8521 "Monochromator_testing.c"

  SIG_MESSAGE("close_E_monitor1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 195 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 195 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 195 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8531 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaclose_E_monitor1);
  rot_transpose(mcrotaClose_wl_monitor1, mctr1);
  rot_mul(mcrotaclose_E_monitor1, mctr1, mcrotrclose_E_monitor1);
  mctc1 = coords_set(
#line 194 "Monochromator_testing.instr"
    0,
#line 194 "Monochromator_testing.instr"
    0,
#line 194 "Monochromator_testing.instr"
    0);
#line 8542 "Monochromator_testing.c"
  rot_transpose(mcrotaClose_wl_monitor1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaclose_E_monitor1 = coords_add(mcposaClose_wl_monitor1, mctc2);
  mctc1 = coords_sub(mcposaClose_wl_monitor1, mcposaclose_E_monitor1);
  mcposrclose_E_monitor1 = rot_apply(mcrotaclose_E_monitor1, mctc1);
  mcDEBUG_COMPONENT("close_E_monitor1", mcposaclose_E_monitor1, mcrotaclose_E_monitor1)
  mccomp_posa[15] = mcposaclose_E_monitor1;
  mccomp_posr[15] = mcposrclose_E_monitor1;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component close_psd_monitor2. */
  /* Setting parameters for component close_psd_monitor2. */
  SIG_MESSAGE("close_psd_monitor2 (Init:SetPar)");
#line 56 "Monochromator_testing.instr"
  mccclose_psd_monitor2_xmin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccclose_psd_monitor2_xmax = 0.05;
#line 56 "Monochromator_testing.instr"
  mccclose_psd_monitor2_ymin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccclose_psd_monitor2_ymax = 0.05;
#line 197 "Monochromator_testing.instr"
  mccclose_psd_monitor2_xwidth = 0.8;
#line 197 "Monochromator_testing.instr"
  mccclose_psd_monitor2_yheight = 0.8;
#line 56 "Monochromator_testing.instr"
  mccclose_psd_monitor2_radius = 0;
#line 8570 "Monochromator_testing.c"

  SIG_MESSAGE("close_psd_monitor2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 199 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 199 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 199 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8580 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaclose_psd_monitor2);
  rot_transpose(mcrotaclose_E_monitor1, mctr1);
  rot_mul(mcrotaclose_psd_monitor2, mctr1, mcrotrclose_psd_monitor2);
  mctc1 = coords_set(
#line 198 "Monochromator_testing.instr"
    0,
#line 198 "Monochromator_testing.instr"
    0,
#line 198 "Monochromator_testing.instr"
    0);
#line 8591 "Monochromator_testing.c"
  rot_transpose(mcrotaclose_E_monitor1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaclose_psd_monitor2 = coords_add(mcposaclose_E_monitor1, mctc2);
  mctc1 = coords_sub(mcposaclose_E_monitor1, mcposaclose_psd_monitor2);
  mcposrclose_psd_monitor2 = rot_apply(mcrotaclose_psd_monitor2, mctc1);
  mcDEBUG_COMPONENT("close_psd_monitor2", mcposaclose_psd_monitor2, mcrotaclose_psd_monitor2)
  mccomp_posa[16] = mcposaclose_psd_monitor2;
  mccomp_posr[16] = mcposrclose_psd_monitor2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component Far_E_monitor2. */
  /* Setting parameters for component Far_E_monitor2. */
  SIG_MESSAGE("Far_E_monitor2 (Init:SetPar)");
#line 54 "Monochromator_testing.instr"
  mccFar_E_monitor2_xmin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccFar_E_monitor2_xmax = 0.05;
#line 54 "Monochromator_testing.instr"
  mccFar_E_monitor2_ymin = -0.05;
#line 54 "Monochromator_testing.instr"
  mccFar_E_monitor2_ymax = 0.05;
#line 202 "Monochromator_testing.instr"
  mccFar_E_monitor2_xwidth = 0.1;
#line 202 "Monochromator_testing.instr"
  mccFar_E_monitor2_yheight = 0.1;
#line 202 "Monochromator_testing.instr"
  mccFar_E_monitor2_Emin = Emin;
#line 202 "Monochromator_testing.instr"
  mccFar_E_monitor2_Emax = Emax;
#line 202 "Monochromator_testing.instr"
  mccFar_E_monitor2_restore_xray = 1;
#line 8623 "Monochromator_testing.c"

  SIG_MESSAGE("Far_E_monitor2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 204 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 204 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 204 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8633 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaFar_E_monitor2);
  rot_transpose(mcrotaclose_psd_monitor2, mctr1);
  rot_mul(mcrotaFar_E_monitor2, mctr1, mcrotrFar_E_monitor2);
  mctc1 = coords_set(
#line 203 "Monochromator_testing.instr"
    0,
#line 203 "Monochromator_testing.instr"
    0,
#line 203 "Monochromator_testing.instr"
    15);
#line 8644 "Monochromator_testing.c"
  rot_transpose(mcrotaMonochromator_Monitor_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFar_E_monitor2 = coords_add(mcposaMonochromator_Monitor_arm, mctc2);
  mctc1 = coords_sub(mcposaclose_psd_monitor2, mcposaFar_E_monitor2);
  mcposrFar_E_monitor2 = rot_apply(mcrotaFar_E_monitor2, mctc1);
  mcDEBUG_COMPONENT("Far_E_monitor2", mcposaFar_E_monitor2, mcrotaFar_E_monitor2)
  mccomp_posa[17] = mcposaFar_E_monitor2;
  mccomp_posr[17] = mcposrFar_E_monitor2;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component Fardivergence_monitor. */
  /* Setting parameters for component Fardivergence_monitor. */
  SIG_MESSAGE("Fardivergence_monitor (Init:SetPar)");
#line 207 "Monochromator_testing.instr"
  if("FarDivergence") strncpy(mccFardivergence_monitor_filename, "FarDivergence" ? "FarDivergence" : "", 16384); else mccFardivergence_monitor_filename[0]='\0';
#line 58 "Monochromator_testing.instr"
  mccFardivergence_monitor_xmin = -0.05;
#line 58 "Monochromator_testing.instr"
  mccFardivergence_monitor_xmax = 0.05;
#line 58 "Monochromator_testing.instr"
  mccFardivergence_monitor_ymin = -0.05;
#line 58 "Monochromator_testing.instr"
  mccFardivergence_monitor_ymax = 0.05;
#line 208 "Monochromator_testing.instr"
  mccFardivergence_monitor_xwidth = 0.1;
#line 209 "Monochromator_testing.instr"
  mccFardivergence_monitor_yheight = 0.1;
#line 59 "Monochromator_testing.instr"
  mccFardivergence_monitor_maxdiv_h = 1;
#line 59 "Monochromator_testing.instr"
  mccFardivergence_monitor_maxdiv_v = 1;
#line 210 "Monochromator_testing.instr"
  mccFardivergence_monitor_restore_xray = 1;
#line 59 "Monochromator_testing.instr"
  mccFardivergence_monitor_nx = 0;
#line 59 "Monochromator_testing.instr"
  mccFardivergence_monitor_ny = 0;
#line 59 "Monochromator_testing.instr"
  mccFardivergence_monitor_nz = 1;
#line 8684 "Monochromator_testing.c"

  SIG_MESSAGE("Fardivergence_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 212 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 212 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 212 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8694 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaFardivergence_monitor);
  rot_transpose(mcrotaFar_E_monitor2, mctr1);
  rot_mul(mcrotaFardivergence_monitor, mctr1, mcrotrFardivergence_monitor);
  mctc1 = coords_set(
#line 211 "Monochromator_testing.instr"
    0,
#line 211 "Monochromator_testing.instr"
    0,
#line 211 "Monochromator_testing.instr"
    0);
#line 8705 "Monochromator_testing.c"
  rot_transpose(mcrotaFar_E_monitor2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFardivergence_monitor = coords_add(mcposaFar_E_monitor2, mctc2);
  mctc1 = coords_sub(mcposaFar_E_monitor2, mcposaFardivergence_monitor);
  mcposrFardivergence_monitor = rot_apply(mcrotaFardivergence_monitor, mctc1);
  mcDEBUG_COMPONENT("Fardivergence_monitor", mcposaFardivergence_monitor, mcrotaFardivergence_monitor)
  mccomp_posa[18] = mcposaFardivergence_monitor;
  mccomp_posr[18] = mcposrFardivergence_monitor;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component Far_wl_monitor2. */
  /* Setting parameters for component Far_wl_monitor2. */
  SIG_MESSAGE("Far_wl_monitor2 (Init:SetPar)");
#line 53 "Monochromator_testing.instr"
  mccFar_wl_monitor2_xmin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccFar_wl_monitor2_xmax = 0.05;
#line 53 "Monochromator_testing.instr"
  mccFar_wl_monitor2_ymin = -0.05;
#line 53 "Monochromator_testing.instr"
  mccFar_wl_monitor2_ymax = 0.05;
#line 214 "Monochromator_testing.instr"
  mccFar_wl_monitor2_xwidth = 0.1;
#line 214 "Monochromator_testing.instr"
  mccFar_wl_monitor2_yheight = 0.1;
#line 214 "Monochromator_testing.instr"
  mccFar_wl_monitor2_Lmin = Lmin;
#line 214 "Monochromator_testing.instr"
  mccFar_wl_monitor2_Lmax = Lmax;
#line 214 "Monochromator_testing.instr"
  mccFar_wl_monitor2_restore_xray = 1;
#line 8737 "Monochromator_testing.c"

  SIG_MESSAGE("Far_wl_monitor2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 216 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 216 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 216 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8747 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaFar_wl_monitor2);
  rot_transpose(mcrotaFardivergence_monitor, mctr1);
  rot_mul(mcrotaFar_wl_monitor2, mctr1, mcrotrFar_wl_monitor2);
  mctc1 = coords_set(
#line 215 "Monochromator_testing.instr"
    0,
#line 215 "Monochromator_testing.instr"
    0,
#line 215 "Monochromator_testing.instr"
    0);
#line 8758 "Monochromator_testing.c"
  rot_transpose(mcrotaFardivergence_monitor, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFar_wl_monitor2 = coords_add(mcposaFardivergence_monitor, mctc2);
  mctc1 = coords_sub(mcposaFardivergence_monitor, mcposaFar_wl_monitor2);
  mcposrFar_wl_monitor2 = rot_apply(mcrotaFar_wl_monitor2, mctc1);
  mcDEBUG_COMPONENT("Far_wl_monitor2", mcposaFar_wl_monitor2, mcrotaFar_wl_monitor2)
  mccomp_posa[19] = mcposaFar_wl_monitor2;
  mccomp_posr[19] = mcposrFar_wl_monitor2;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component Far_psd_monitor2. */
  /* Setting parameters for component Far_psd_monitor2. */
  SIG_MESSAGE("Far_psd_monitor2 (Init:SetPar)");
#line 56 "Monochromator_testing.instr"
  mccFar_psd_monitor2_xmin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccFar_psd_monitor2_xmax = 0.05;
#line 56 "Monochromator_testing.instr"
  mccFar_psd_monitor2_ymin = -0.05;
#line 56 "Monochromator_testing.instr"
  mccFar_psd_monitor2_ymax = 0.05;
#line 219 "Monochromator_testing.instr"
  mccFar_psd_monitor2_xwidth = 0.1;
#line 219 "Monochromator_testing.instr"
  mccFar_psd_monitor2_yheight = 0.1;
#line 56 "Monochromator_testing.instr"
  mccFar_psd_monitor2_radius = 0;
#line 8786 "Monochromator_testing.c"

  SIG_MESSAGE("Far_psd_monitor2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 221 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 221 "Monochromator_testing.instr"
    (0)*DEG2RAD,
#line 221 "Monochromator_testing.instr"
    (90)*DEG2RAD);
#line 8796 "Monochromator_testing.c"
  rot_mul(mctr1, mcrotaMonochromator_Monitor_arm, mcrotaFar_psd_monitor2);
  rot_transpose(mcrotaFar_wl_monitor2, mctr1);
  rot_mul(mcrotaFar_psd_monitor2, mctr1, mcrotrFar_psd_monitor2);
  mctc1 = coords_set(
#line 220 "Monochromator_testing.instr"
    0,
#line 220 "Monochromator_testing.instr"
    0,
#line 220 "Monochromator_testing.instr"
    0);
#line 8807 "Monochromator_testing.c"
  rot_transpose(mcrotaFar_wl_monitor2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFar_psd_monitor2 = coords_add(mcposaFar_wl_monitor2, mctc2);
  mctc1 = coords_sub(mcposaFar_wl_monitor2, mcposaFar_psd_monitor2);
  mcposrFar_psd_monitor2 = rot_apply(mcrotaFar_psd_monitor2, mctc1);
  mcDEBUG_COMPONENT("Far_psd_monitor2", mcposaFar_psd_monitor2, mcrotaFar_psd_monitor2)
  mccomp_posa[20] = mcposaFar_psd_monitor2;
  mccomp_posr[20] = mcposrFar_psd_monitor2;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
  /* Component initializations. */
  /* Initializations for component origin. */
  SIG_MESSAGE("origin (Init)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define profile mccorigin_profile
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 63 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../misc/Progress_bar.comp"
{
  fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 8838 "Monochromator_testing.c"
#undef minutes
#undef flag_save
#undef percent
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef profile
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component source_flat. */
  SIG_MESSAGE("source_flat (Init)");
#define mccompcurname  source_flat
#define mccompcurtype  Source_flat
#define mccompcurindex 2
#define spectrum_file mccsource_flat_spectrum_file
#define prms mccsource_flat_prms
#define square mccsource_flat_square
#define radius mccsource_flat_radius
#define yheight mccsource_flat_yheight
#define xwidth mccsource_flat_xwidth
#define xmin mccsource_flat_xmin
#define xmax mccsource_flat_xmax
#define dist mccsource_flat_dist
#define focus_xw mccsource_flat_focus_xw
#define focus_yh mccsource_flat_focus_yh
#define E0 mccsource_flat_E0
#define dE mccsource_flat_dE
#define lambda0 mccsource_flat_lambda0
#define dlambda mccsource_flat_dlambda
#define flux mccsource_flat_flux
#define gauss mccsource_flat_gauss
#define randomphase mccsource_flat_randomphase
#define phase mccsource_flat_phase
#line 73 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
{
  square = 0;srcArea=0;
  /* Determine source area */

  if (xmin && xmax && !xwidth){
    xwidth=xmax-xmin;
  }

  if (radius && !yheight && !xwidth ) {
    square = 0;
    srcArea = PI*radius*radius;
  } else if(yheight && xwidth) {
    square = 1;
    srcArea = xwidth * yheight;
  }

  
  if (srcArea <= 0) {
    printf("Source_flat: %s: Source area is <= 0 !\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  if (dist <= 0 || focus_xw <= 0 || focus_yh <= 0) {
    printf("Source_flat: %s: Target area unmeaningful! (negative dist / focus_xw / focus_yh)\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  
  if (spectrum_file){
    /*read spectrum from file*/
    int status=0;
    if ( (status=Table_Read(&(prms.T),spectrum_file,0))==-1){
      fprintf(stderr,"Source_flat(%s) Error: Could not parse file \"%s\"\n",NAME_CURRENT_COMP,spectrum_file?spectrum_file:"");
      exit(-1);
    }
    /*data is now in table t*/
    /*integrate to get total flux, assuming raw numbers have been corrected for measuring aperture*/
    int i;
    prms.pint=0;
    t_Table *T=&(prms.T);
    for (i=0;i<prms.T.rows-1;i++){
      prms.pint+=((T->data[i*T->columns+1]+T->data[(i+1)*T->columns+1])/2.0)*fabs(T->data[(i+1)*T->columns]-T->data[i*T->columns]); 
    }
    printf("Source_flat(%s) Integrated intensity radiated is %g pht/s\n",NAME_CURRENT_COMP,prms.pint);
    if(E0) printf("Source_flat(%s) E0!=0 -> assuming intensity spectrum is parametrized by energy [keV]\n",NAME_CURRENT_COMP);
  }else if ( !E0 && !lambda0){
    fprintf(stderr,"Source_flat(%s): Error: Must specify either wavelength or energy distribution\n",NAME_CURRENT_COMP); 
    exit(0);  
  }  
  /*calculate the X-ray weight from the flux*/
  if (flux){
    prms.pmul=flux;
  }else{
    prms.pmul=1;
  }
  prms.pmul*=1.0/((double) mcget_ncount());
}
#line 8932 "Monochromator_testing.c"
#undef phase
#undef randomphase
#undef gauss
#undef flux
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_yh
#undef focus_xw
#undef dist
#undef xmax
#undef xmin
#undef xwidth
#undef yheight
#undef radius
#undef square
#undef prms
#undef spectrum_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component dmu. */
  SIG_MESSAGE("dmu (Init)");
#define mccompcurname  dmu
#define mccompcurtype  Undulator
#define mccompcurindex 3
#define verbose mccdmu_verbose
#define prms mccdmu_prms
#define alpha mccdmu_alpha
#define MELE mccdmu_MELE
#define E0 mccdmu_E0
#define dE mccdmu_dE
#define phase mccdmu_phase
#define randomphase mccdmu_randomphase
#define Ee mccdmu_Ee
#define dEe mccdmu_dEe
#define Ie mccdmu_Ie
#define tbunch mccdmu_tbunch
#define t0 mccdmu_t0
#define B mccdmu_B
#define K mccdmu_K
#define gap mccdmu_gap
#define Nper mccdmu_Nper
#define lu mccdmu_lu
#define sigey mccdmu_sigey
#define sigex mccdmu_sigex
#define sigepx mccdmu_sigepx
#define sigepy mccdmu_sigepy
#define focus_xw mccdmu_focus_xw
#define focus_yh mccdmu_focus_yh
#define dist mccdmu_dist
#define gauss_t mccdmu_gauss_t
#define quick_integ mccdmu_quick_integ
#define E1st mccdmu_E1st
#line 128 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"
{
  fprintf(stderr,"Warning (%s): Undulator is an experimental component - testing is ongoing\n",NAME_CURRENT_COMP);
  
  prms.length=lu*Nper;

  if( (!E1st && !K && B<=0) || Ee<=0 || Ie<=0 ){
    fprintf(stderr, "Error (%s): E1st, K, B, Ee, and Ie do not have a sane set of values. Found (%g %g %g %g %g). Aborting.\n",NAME_CURRENT_COMP,E1st,K,B,Ee,Ie);
    exit(1);
  }

  /*compute gamma*/
  prms.gamma=(Ee*1e9)/(MELE/CELE*M_C*M_C);/*the extra CELE is to convert to eV*/
  prms.gamma2=prms.gamma*prms.gamma;
  prms.igamma=1.0/prms.gamma;

  if(E1st){
      /*compute K and B from desired target energy*/
      K=sqrt(2*(4*M_PI*prms.gamma2/(E2K*E1st*lu*1e10) -1));
      B=2*M_PI*MELE*M_C*K/CELE/lu;
      if (verbose) printf("Undulator (%s) (K,B)=f(E1st=%g)= (%g , %g)\n",NAME_CURRENT_COMP,E1st,K,B);
  }else if(B>0){
    K=CELE*B*lu/(2*M_PI*MELE*M_C);
    if (verbose) printf("Undulator (%s) K=f(B=%g)= %g\n",NAME_CURRENT_COMP,B,K);
  }else if (K){
    B=2*M_PI*MELE*M_C*K/CELE/lu;
    if (verbose) printf("Undulator (%s) B=f(K=%g)= %g\n",NAME_CURRENT_COMP,K,B);
  }

  if (sigex <0 || sigey<0){
   fprintf(stderr, "Error (%s): sigex and sigey must be >= 0. Negative beam size isn't meaningful. Aborting.\n",NAME_CURRENT_COMP,sigex,sigey);
    exit(1);
  }
  if (dist<=0 || focus_xw < 0 || focus_yh < 0){
    fprintf(stderr,"Error (%s): Target undefined.\n",NAME_CURRENT_COMP);
    exit(1);
  }


  //printf("Undulator (%s): gamma=%g, divergence is 1/gamma=%g rad.\n",NAME_CURRENT_COMP,prms.gamma,prms.igamma);
  /*compute characteristic energy in keV*/
  double Ec=0.665*Ee*Ee*B;
  //double Ec=1.5*prms.gamma2*HBAR*CELE*B/MELE *1e-3; /*check units on this one. The 1e-3 factor is because energy is assumed to be in keV*/
  /*We normally do computations in k so use that for transfer*/
  prms.kc=E2K*Ec;

  /*allocate an integration workspace*/
  if(!quick_integ){
      prms.gsl_int_ws = gsl_integration_workspace_alloc (1000);
  }
  prms.Bsig.function= &mxundulator_Bsig_integrand;
  prms.Bpi.function = &mxundulator_Bpi_integrand;
  gsl_set_error_handler_off();

  /*correct for number of rays*/
  prms.pmul=1.0/( (double) mcget_ncount());

  /*correct for finite energy interval*/
  if(dE){
      prms.pmul*=dE*2.0;
  }

}
#line 9052 "Monochromator_testing.c"
#undef E1st
#undef quick_integ
#undef gauss_t
#undef dist
#undef focus_yh
#undef focus_xw
#undef sigepy
#undef sigepx
#undef sigex
#undef sigey
#undef lu
#undef Nper
#undef gap
#undef K
#undef B
#undef t0
#undef tbunch
#undef Ie
#undef dEe
#undef Ee
#undef randomphase
#undef phase
#undef dE
#undef E0
#undef MELE
#undef alpha
#undef prms
#undef verbose
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mbefore_E_monitor. */
  SIG_MESSAGE("Mbefore_E_monitor (Init)");
#define mccompcurname  Mbefore_E_monitor
#define mccompcurtype  E_monitor
#define mccompcurindex 4
#define nE mccMbefore_E_monitor_nE
#define filename mccMbefore_E_monitor_filename
#define E_N mccMbefore_E_monitor_E_N
#define E_p mccMbefore_E_monitor_E_p
#define E_p2 mccMbefore_E_monitor_E_p2
#define xmin mccMbefore_E_monitor_xmin
#define xmax mccMbefore_E_monitor_xmax
#define ymin mccMbefore_E_monitor_ymin
#define ymax mccMbefore_E_monitor_ymax
#define xwidth mccMbefore_E_monitor_xwidth
#define yheight mccMbefore_E_monitor_yheight
#define Emin mccMbefore_E_monitor_Emin
#define Emax mccMbefore_E_monitor_Emax
#define restore_xray mccMbefore_E_monitor_restore_xray
#line 64 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
}
#line 9125 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mbefore_wl_monitor. */
  SIG_MESSAGE("Mbefore_wl_monitor (Init)");
#define mccompcurname  Mbefore_wl_monitor
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccMbefore_wl_monitor_nL
#define filename mccMbefore_wl_monitor_filename
#define L_N mccMbefore_wl_monitor_L_N
#define L_p mccMbefore_wl_monitor_L_p
#define L_p2 mccMbefore_wl_monitor_L_p2
#define xmin mccMbefore_wl_monitor_xmin
#define xmax mccMbefore_wl_monitor_xmax
#define ymin mccMbefore_wl_monitor_ymin
#define ymax mccMbefore_wl_monitor_ymax
#define xwidth mccMbefore_wl_monitor_xwidth
#define yheight mccMbefore_wl_monitor_yheight
#define Lmin mccMbefore_wl_monitor_Lmin
#define Lmax mccMbefore_wl_monitor_Lmax
#define restore_xray mccMbefore_wl_monitor_restore_xray
#line 63 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 9184 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mbefore_psd_monitor. */
  SIG_MESSAGE("Mbefore_psd_monitor (Init)");
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
#define xmin mccMbefore_psd_monitor_xmin
#define xmax mccMbefore_psd_monitor_xmax
#define ymin mccMbefore_psd_monitor_ymin
#define ymax mccMbefore_psd_monitor_ymax
#define xwidth mccMbefore_psd_monitor_xwidth
#define yheight mccMbefore_psd_monitor_yheight
#define radius mccMbefore_psd_monitor_radius
#line 67 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;
    double *p1,*p2,*p3;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ( ((xmin >= xmax) || (ymin >= ymax)) && !radius ) {
            printf("PSD_monitor: %s: Null detection area !\n"
                "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax,radius). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
    if(!radius){
      p1=calloc(nx*ny,sizeof(double));
      p2=calloc(nx*ny,sizeof(double));
      p3=calloc(nx*ny,sizeof(double));
      
      PSD_N=calloc(nx,sizeof(double *));
      PSD_p=calloc(nx,sizeof(double *));
      PSD_p2=calloc(nx,sizeof(double *));

      for (i=0; i<nx; i++){
        PSD_N[i]=&(p1[i*ny]);//calloc(ny,sizeof(double));
        PSD_p[i]=&(p2[i*ny]);//calloc(ny,sizeof(double));
        PSD_p2[i]=&(p3[i*ny]);//calloc(ny,sizeof(double));
      }
    }else{
      PSD_N=calloc(1,sizeof(double *));
      PSD_p=calloc(1,sizeof(double *));
      PSD_p2=calloc(1,sizeof(double *));
      *PSD_N=calloc(nr,sizeof(double));
      *PSD_p=calloc(nr,sizeof(double));
      *PSD_p2=calloc(nr,sizeof(double));
    }

}
#line 9261 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PG1_arm. */
  SIG_MESSAGE("PG1_arm (Init)");

  /* Initializations for component M2_rotation_arm1. */
  SIG_MESSAGE("M2_rotation_arm1 (Init)");

  /* Initializations for component M2_rotation_arm2. */
  SIG_MESSAGE("M2_rotation_arm2 (Init)");

  /* Initializations for component M2_rotation_arm3. */
  SIG_MESSAGE("M2_rotation_arm3 (Init)");

  /* Initializations for component mirror2. */
  SIG_MESSAGE("mirror2 (Init)");
#define mccompcurname  mirror2
#define mccompcurtype  Mirror
#define mccompcurindex 11
#define reflec mccmirror2_reflec
#define reflec_table mccmirror2_reflec_table
#define prms mccmirror2_prms
#define zdepth mccmirror2_zdepth
#define xwidth mccmirror2_xwidth
#define R0 mccmirror2_R0
#line 57 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Mirror.comp"
{
    if (reflec && strlen(reflec)) {
        char **header_parsed;
        t_Table *tp=&reflec_table;
        /* read 1st block data from file into tp */
        if (Table_Read(tp, reflec, 1) <= 0)
        {
            exit(fprintf(stderr,"Error: %s: cannot read file %s\n",NAME_CURRENT_COMP, reflec));
        }
        header_parsed = Table_ParseHeader(tp->header,
                "e_min=","e_max=","e_step=","theta_min=","theta_max=","theta_step=",NULL);
        if (header_parsed[0] && header_parsed[1] && header_parsed[2] &&
                header_parsed[3] && header_parsed[4] && header_parsed[5])
        {
            prms.e_min=strtod(header_parsed[0],NULL);
            prms.e_max=strtod(header_parsed[1],NULL);
            prms.e_step=strtod(header_parsed[2],NULL);
            prms.theta_min=strtod(header_parsed[3],NULL);
            prms.theta_max=strtod(header_parsed[4],NULL);
            prms.theta_step=strtod(header_parsed[5],NULL);
        } else {
            exit(fprintf(stderr,"Error: %s: wrong/missing header line(s) in file %s\n", NAME_CURRENT_COMP, reflec));
        }
        if ((prms.e_max-prms.e_min) != ((tp->rows-1)*prms.e_step) )
        {
            exit(fprintf(stderr,"Error: %s: e_step does not match e_min and e_max in file %s\n",NAME_CURRENT_COMP, reflec));
        }
        if ((prms.theta_max-prms.theta_min) != ((tp->columns-1)*prms.theta_step) )
        {
            exit(fprintf(stderr,"Error: %s: theta_step does not match theta_min and theta_max in file %s\n",NAME_CURRENT_COMP, reflec));
        }
        prms.use_reflec_table=1;
    }else{
        prms.use_reflec_table=0;
    }
}
#line 9341 "Monochromator_testing.c"
#undef R0
#undef xwidth
#undef zdepth
#undef prms
#undef reflec_table
#undef reflec
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Plane_grating. */
  SIG_MESSAGE("Plane_grating (Init)");
#define mccompcurname  Plane_grating
#define mccompcurtype  MultiPurpose_grating
#define mccompcurindex 12
#define blazed mccPlane_grating_blazed
#define display mccPlane_grating_display
#define zdepth mccPlane_grating_zdepth
#define xwidth mccPlane_grating_xwidth
#define MCangleVariation mccPlane_grating_MCangleVariation
#define R0 mccPlane_grating_R0
#define r_rho mccPlane_grating_r_rho
#define b mccPlane_grating_b
#define N_slits mccPlane_grating_N_slits
#define d mccPlane_grating_d
#define blazed_angle mccPlane_grating_blazed_angle
#line 96 "MultiPurpose_grating.comp"
{
    int status;
if (MCangleVariation>15){
       if(display){      
          printf("\n Warning(%s) : MCangleVariation of:%f is high.\n",NAME_CURRENT_COMP,MCangleVariation); 
       } 
}
if (N_slits<0 || b<0 || d<0 || r_rho<0 || MCangleVariation<0 || blazed<0 || blazed_angle<0 || zdepth<0 || xwidth<0 || blazed<0 || blazed_angle<0){
       printf("Error(%s) : Negative input parameter given.\n",NAME_CURRENT_COMP);
       exit(-1);
}
if (!blazed && blazed_angle){
       printf("Error(%s) : A blazed angle is given, but not a blazing mode. \n",NAME_CURRENT_COMP);
       exit(-1);
}
if (!r_rho){
       printf("Error(%s) : Need line density [l/mm] to define grating. \n",NAME_CURRENT_COMP);
       exit(-1);
}
if (R0<0 || R0>1){
      printf("Error(%s) reflectivity (%f) is specified but is not in [0:1]. \n",NAME_CURRENT_COMP,R0);
      exit(-1);
}

      if(!d){
      // Angstrom/line
      d=10000000/r_rho; 
      }
      if(!N_slits){
      // Number of slits 
      N_slits=(zdepth*1000)*r_rho; 
      } 
      if(!b){
      // Approximated slit-width in Angstrom
      b=d/3;
      }
      if(display){      
      printf("\n Line width, d=%f [AA]. Number of slits,  N=%f. Slit width, b=%f [AA].\n \n",d,N_slits,b);  
      }   

DirectionalSamplingFactor=((2*(MCangleVariation*DEG2RAD))/(4*PI));       
}
#line 9411 "Monochromator_testing.c"
#undef blazed_angle
#undef d
#undef N_slits
#undef b
#undef r_rho
#undef R0
#undef MCangleVariation
#undef xwidth
#undef zdepth
#undef display
#undef blazed
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Monochromator_Monitor_arm. */
  SIG_MESSAGE("Monochromator_Monitor_arm (Init)");

  /* Initializations for component Close_wl_monitor1. */
  SIG_MESSAGE("Close_wl_monitor1 (Init)");
#define mccompcurname  Close_wl_monitor1
#define mccompcurtype  L_monitor
#define mccompcurindex 14
#define nL mccClose_wl_monitor1_nL
#define filename mccClose_wl_monitor1_filename
#define L_N mccClose_wl_monitor1_L_N
#define L_p mccClose_wl_monitor1_L_p
#define L_p2 mccClose_wl_monitor1_L_p2
#define xmin mccClose_wl_monitor1_xmin
#define xmax mccClose_wl_monitor1_xmax
#define ymin mccClose_wl_monitor1_ymin
#define ymax mccClose_wl_monitor1_ymax
#define xwidth mccClose_wl_monitor1_xwidth
#define yheight mccClose_wl_monitor1_yheight
#define Lmin mccClose_wl_monitor1_Lmin
#define Lmax mccClose_wl_monitor1_Lmax
#define restore_xray mccClose_wl_monitor1_restore_xray
#line 63 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 9470 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component close_E_monitor1. */
  SIG_MESSAGE("close_E_monitor1 (Init)");
#define mccompcurname  close_E_monitor1
#define mccompcurtype  E_monitor
#define mccompcurindex 15
#define nE mccclose_E_monitor1_nE
#define filename mccclose_E_monitor1_filename
#define E_N mccclose_E_monitor1_E_N
#define E_p mccclose_E_monitor1_E_p
#define E_p2 mccclose_E_monitor1_E_p2
#define xmin mccclose_E_monitor1_xmin
#define xmax mccclose_E_monitor1_xmax
#define ymin mccclose_E_monitor1_ymin
#define ymax mccclose_E_monitor1_ymax
#define xwidth mccclose_E_monitor1_xwidth
#define yheight mccclose_E_monitor1_yheight
#define Emin mccclose_E_monitor1_Emin
#define Emax mccclose_E_monitor1_Emax
#define restore_xray mccclose_E_monitor1_restore_xray
#line 64 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
}
#line 9529 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component close_psd_monitor2. */
  SIG_MESSAGE("close_psd_monitor2 (Init)");
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
#define xmin mccclose_psd_monitor2_xmin
#define xmax mccclose_psd_monitor2_xmax
#define ymin mccclose_psd_monitor2_ymin
#define ymax mccclose_psd_monitor2_ymax
#define xwidth mccclose_psd_monitor2_xwidth
#define yheight mccclose_psd_monitor2_yheight
#define radius mccclose_psd_monitor2_radius
#line 67 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;
    double *p1,*p2,*p3;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ( ((xmin >= xmax) || (ymin >= ymax)) && !radius ) {
            printf("PSD_monitor: %s: Null detection area !\n"
                "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax,radius). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
    if(!radius){
      p1=calloc(nx*ny,sizeof(double));
      p2=calloc(nx*ny,sizeof(double));
      p3=calloc(nx*ny,sizeof(double));
      
      PSD_N=calloc(nx,sizeof(double *));
      PSD_p=calloc(nx,sizeof(double *));
      PSD_p2=calloc(nx,sizeof(double *));

      for (i=0; i<nx; i++){
        PSD_N[i]=&(p1[i*ny]);//calloc(ny,sizeof(double));
        PSD_p[i]=&(p2[i*ny]);//calloc(ny,sizeof(double));
        PSD_p2[i]=&(p3[i*ny]);//calloc(ny,sizeof(double));
      }
    }else{
      PSD_N=calloc(1,sizeof(double *));
      PSD_p=calloc(1,sizeof(double *));
      PSD_p2=calloc(1,sizeof(double *));
      *PSD_N=calloc(nr,sizeof(double));
      *PSD_p=calloc(nr,sizeof(double));
      *PSD_p2=calloc(nr,sizeof(double));
    }

}
#line 9606 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Far_E_monitor2. */
  SIG_MESSAGE("Far_E_monitor2 (Init)");
#define mccompcurname  Far_E_monitor2
#define mccompcurtype  E_monitor
#define mccompcurindex 17
#define nE mccFar_E_monitor2_nE
#define filename mccFar_E_monitor2_filename
#define E_N mccFar_E_monitor2_E_N
#define E_p mccFar_E_monitor2_E_p
#define E_p2 mccFar_E_monitor2_E_p2
#define xmin mccFar_E_monitor2_xmin
#define xmax mccFar_E_monitor2_xmax
#define ymin mccFar_E_monitor2_ymin
#define ymax mccFar_E_monitor2_ymax
#define xwidth mccFar_E_monitor2_xwidth
#define yheight mccFar_E_monitor2_yheight
#define Emin mccFar_E_monitor2_Emin
#define Emax mccFar_E_monitor2_Emax
#define restore_xray mccFar_E_monitor2_restore_xray
#line 64 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
}
#line 9666 "Monochromator_testing.c"
#undef restore_xray
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fardivergence_monitor. */
  SIG_MESSAGE("Fardivergence_monitor (Init)");
#define mccompcurname  Fardivergence_monitor
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 18
#define nh mccFardivergence_monitor_nh
#define nv mccFardivergence_monitor_nv
#define rad mccFardivergence_monitor_rad
#define Div_N mccFardivergence_monitor_Div_N
#define Div_p mccFardivergence_monitor_Div_p
#define Div_p2 mccFardivergence_monitor_Div_p2
#define filename mccFardivergence_monitor_filename
#define xmin mccFardivergence_monitor_xmin
#define xmax mccFardivergence_monitor_xmax
#define ymin mccFardivergence_monitor_ymin
#define ymax mccFardivergence_monitor_ymax
#define xwidth mccFardivergence_monitor_xwidth
#define yheight mccFardivergence_monitor_yheight
#define maxdiv_h mccFardivergence_monitor_maxdiv_h
#define maxdiv_v mccFardivergence_monitor_maxdiv_v
#define restore_xray mccFardivergence_monitor_restore_xray
#define nx mccFardivergence_monitor_nx
#define ny mccFardivergence_monitor_ny
#define nz mccFardivergence_monitor_nz
#line 69 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
    int i,j;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 9732 "Monochromator_testing.c"
#undef nz
#undef ny
#undef nx
#undef restore_xray
#undef maxdiv_v
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef rad
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Far_wl_monitor2. */
  SIG_MESSAGE("Far_wl_monitor2 (Init)");
#define mccompcurname  Far_wl_monitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mccFar_wl_monitor2_nL
#define filename mccFar_wl_monitor2_filename
#define L_N mccFar_wl_monitor2_L_N
#define L_p mccFar_wl_monitor2_L_p
#define L_p2 mccFar_wl_monitor2_L_p2
#define xmin mccFar_wl_monitor2_xmin
#define xmax mccFar_wl_monitor2_xmax
#define ymin mccFar_wl_monitor2_ymin
#define ymax mccFar_wl_monitor2_ymax
#define xwidth mccFar_wl_monitor2_xwidth
#define yheight mccFar_wl_monitor2_yheight
#define Lmin mccFar_wl_monitor2_Lmin
#define Lmax mccFar_wl_monitor2_Lmax
#define restore_xray mccFar_wl_monitor2_restore_xray
#line 63 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 9796 "Monochromator_testing.c"
#undef restore_xray
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Far_psd_monitor2. */
  SIG_MESSAGE("Far_psd_monitor2 (Init)");
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
#define xmin mccFar_psd_monitor2_xmin
#define xmax mccFar_psd_monitor2_xmax
#define ymin mccFar_psd_monitor2_ymin
#define ymax mccFar_psd_monitor2_ymax
#define xwidth mccFar_psd_monitor2_xwidth
#define yheight mccFar_psd_monitor2_yheight
#define radius mccFar_psd_monitor2_radius
#line 67 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;
    double *p1,*p2,*p3;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ( ((xmin >= xmax) || (ymin >= ymax)) && !radius ) {
            printf("PSD_monitor: %s: Null detection area !\n"
                "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax,radius). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
    if(!radius){
      p1=calloc(nx*ny,sizeof(double));
      p2=calloc(nx*ny,sizeof(double));
      p3=calloc(nx*ny,sizeof(double));
      
      PSD_N=calloc(nx,sizeof(double *));
      PSD_p=calloc(nx,sizeof(double *));
      PSD_p2=calloc(nx,sizeof(double *));

      for (i=0; i<nx; i++){
        PSD_N[i]=&(p1[i*ny]);//calloc(ny,sizeof(double));
        PSD_p[i]=&(p2[i*ny]);//calloc(ny,sizeof(double));
        PSD_p2[i]=&(p3[i*ny]);//calloc(ny,sizeof(double));
      }
    }else{
      PSD_N=calloc(1,sizeof(double *));
      PSD_p=calloc(1,sizeof(double *));
      PSD_p2=calloc(1,sizeof(double *));
      *PSD_N=calloc(nr,sizeof(double));
      *PSD_p=calloc(nr,sizeof(double));
      *PSD_p2=calloc(nr,sizeof(double));
    }

}
#line 9873 "Monochromator_testing.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Copy xray state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlkx = mcnkx;
  MCNUM mcnlky = mcnky;
  MCNUM mcnlkz = mcnkz;
  MCNUM mcnlphi = mcnphi;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlEx = mcnEx;
  MCNUM mcnlEy = mcnEy;
  MCNUM mcnlEz = mcnEz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component origin [1] */
  mccoordschange(mcposrorigin, mcrotrorigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component origin (without coords transformations) */
  mcJumpTrace_origin:
  SIG_MESSAGE("origin (Trace)");
  mcDEBUG_COMP("origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbComporigin
  STORE_XRAY(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define profile mccorigin_profile
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 71 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../misc/Progress_bar.comp"
{
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  if (EndTime &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;
    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    if (IntermediateCnts >= mcget_ncount()) fprintf(stdout, "\n");
    if (flag_save) mcsave(NULL);
  }
}
#line 10038 "Monochromator_testing.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef profile
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbComporigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component source_flat [2] */
  mccoordschange(mcposrsource_flat, mcrotrsource_flat,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component source_flat (without coords transformations) */
  mcJumpTrace_source_flat:
  SIG_MESSAGE("source_flat (Trace)");
  mcDEBUG_COMP("source_flat")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource_flat
  STORE_XRAY(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  source_flat
#define mccompcurtype  Source_flat
#define mccompcurindex 2
#define spectrum_file mccsource_flat_spectrum_file
#define prms mccsource_flat_prms
#define square mccsource_flat_square
{   /* Declarations of source_flat=Source_flat() SETTING parameters. */
MCNUM radius = mccsource_flat_radius;
MCNUM yheight = mccsource_flat_yheight;
MCNUM xwidth = mccsource_flat_xwidth;
MCNUM xmin = mccsource_flat_xmin;
MCNUM xmax = mccsource_flat_xmax;
MCNUM dist = mccsource_flat_dist;
MCNUM focus_xw = mccsource_flat_focus_xw;
MCNUM focus_yh = mccsource_flat_focus_yh;
MCNUM E0 = mccsource_flat_E0;
MCNUM dE = mccsource_flat_dE;
MCNUM lambda0 = mccsource_flat_lambda0;
MCNUM dlambda = mccsource_flat_dlambda;
MCNUM flux = mccsource_flat_flux;
MCNUM gauss = mccsource_flat_gauss;
MCNUM randomphase = mccsource_flat_randomphase;
MCNUM phase = mccsource_flat_phase;
/* 'source_flat=Source_flat()' component instance has conditional execution */
if (( ! mcipSourceChoice ))

#line 133 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
{
  double chi,e,k,l,r, xf, yf, rf, dx, dy, pdir;

  phi=0;
  z=0;
  if (square == 1) {
    if (xmin && xmax){
      x = xmin + xwidth*rand01();
    }else{
      x = xwidth * (rand01() - 0.5);
    }
    y = yheight * (rand01() - 0.5);
  } else {
    chi=2*PI*rand01();                          /* Choose point on source */
    r=sqrt(rand01())*radius;                    /* with uniform distribution. */
    x=r*cos(chi);
    y=r*sin(chi);
  }
  randvec_target_rect_real(&xf, &yf, &rf, &pdir,
      0, 0, dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

  dx = xf-x;
  dy = yf-y;
  rf = sqrt(dx*dx+dy*dy+dist*dist);
  /*pdir contains the unnormalized solid angle weighting */
  p = prms.pmul*pdir/(4*M_PI);

  if (spectrum_file){
    double pp=0;
    //while (pp<=0){ 
    l=prms.T.data[0]+ (prms.T.data[(prms.T.rows-1)*prms.T.columns] -prms.T.data[0])*rand01();
    pp=Table_Value(prms.T,l,1);
    //}
    p*=pp;
    /*if E0!=0 the tabled value is assumed to energy in keV*/
    if (E0!=0){
      k=E2K*l;
    }else{
      k=(2*M_PI/l);
    }
  }else if (E0){
    if(!dE){
      e=E0;
    }else if (gauss){
      e=E0+dE*randnorm();
    }else{
      e=randpm1()*dE + E0;
    }
    k=E2K*e;
  }else if (lambda0){
    if (!dlambda){
      l=lambda0;
    }else if (gauss){
      l=lambda0+dlambda*randnorm();
    }else{
      l=randpm1()*dlambda*0.5 + lambda0;
    }
    k=(2*M_PI/l);
  }
  ky=k*dy/rf;
  kx=k*dx/rf;
  kz=k*dist/rf;

  /*randomly pick phase or set to something real*/
  if (randomphase){
    phi=rand01()*2*M_PI;
  }else{
    phi=phase;
  }

  /*set polarization vector*/
  Ex=0;Ey=0;Ez=0;
}
#line 10249 "Monochromator_testing.c"
}   /* End of source_flat=Source_flat() SETTING parameter declarations. */
#undef square
#undef prms
#undef spectrum_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompsource_flat:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component dmu [3] */
  mccoordschange(mcposrdmu, mcrotrdmu,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component dmu (without coords transformations) */
  mcJumpTrace_dmu:
  SIG_MESSAGE("dmu (Trace)");
  mcDEBUG_COMP("dmu")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdmu
  STORE_XRAY(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  dmu
#define mccompcurtype  Undulator
#define mccompcurindex 3
#define verbose mccdmu_verbose
#define prms mccdmu_prms
#define alpha mccdmu_alpha
#define MELE mccdmu_MELE
{   /* Declarations of dmu=Undulator() SETTING parameters. */
MCNUM E0 = mccdmu_E0;
MCNUM dE = mccdmu_dE;
MCNUM phase = mccdmu_phase;
MCNUM randomphase = mccdmu_randomphase;
MCNUM Ee = mccdmu_Ee;
MCNUM dEe = mccdmu_dEe;
MCNUM Ie = mccdmu_Ie;
MCNUM tbunch = mccdmu_tbunch;
MCNUM t0 = mccdmu_t0;
MCNUM B = mccdmu_B;
MCNUM K = mccdmu_K;
MCNUM gap = mccdmu_gap;
int Nper = mccdmu_Nper;
MCNUM lu = mccdmu_lu;
MCNUM sigey = mccdmu_sigey;
MCNUM sigex = mccdmu_sigex;
MCNUM sigepx = mccdmu_sigepx;
MCNUM sigepy = mccdmu_sigepy;
MCNUM focus_xw = mccdmu_focus_xw;
MCNUM focus_yh = mccdmu_focus_yh;
MCNUM dist = mccdmu_dist;
MCNUM gauss_t = mccdmu_gauss_t;
MCNUM quick_integ = mccdmu_quick_integ;
MCNUM E1st = mccdmu_E1st;
/* 'dmu=Undulator()' component instance has conditional execution */
if (( mcipSourceChoice ))

#line 193 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"
{

    double k,e,l,w_u,r,w;
    double xo,yo,zo,xi,psi,theta2,Omega;
    double bsigma_integral, bpi_integral, s_n2;
    double bsigma_error, bpi_error;

    /* pick an energy in the given interval */
    e=E0+randpm1()*dE;

    /* add electron beam energy spread to gamma parameters (if necessary).*/
    if( dEe){
        double deltaEe=(randnorm()*dEe*Ee)+ Ee;
        prms.gamma=(deltaEe*1e9)/(MELE/CELE*M_C*M_C);/*the extra CELE is to convert to eV*/
        prms.gamma2=prms.gamma*prms.gamma;
        prms.igamma=1.0/prms.gamma;
    }

    /* the undulator's fundamental wavelength*/
    w_u=2*M_PI*M_C/lu;

    /* now pick angles (\xi and \psi => \theta) within the focus_window 
     * (point on window + point in electron beam)
     * ... we can now weight properly according to ang. flux density*/
    if(focus_xw && focus_yh){
        randvec_target_rect(&xo,&yo,&zo, &Omega, 0,0,dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP);
    }else if (focus_yh!=0){
        xi=0;xo=0;yo=0;
        yo=randpm1()*0.5*focus_yh;
        Omega=1;
    }else if (focus_xw!=0){
        psi=0;xo=0;yo=0;
        xo=randpm1()*0.5*focus_xw;
        Omega=1;
    }else{
        xi=0;psi=0;xo=0;yo=0;
        Omega=1;
    }
    p=prms.pmul*Omega/(4*M_PI);

    /*add emittance effects - note that doing it this way will shoot some rays outside the focus-window*/
    x=y=z=0;
    if(sigex){
        x=randnorm()*sigex;
    }
    if(sigey){
        y=randnorm()*sigey;
    }

    xi=fabs(tan(xo/dist));
    psi=fabs(tan(yo/dist));
    /* This has to be after (xi,psi), else it will be convoluted into the weight calculation.*/
    if(sigepx){
        xo+=randnorm()*sigepx*dist;
    }
    if(sigepy){
        yo+=randnorm()*sigepy*dist;
    }

    theta2=xi*xi+psi*psi;

    k=E2K*e;
    /*anguar frequency w:*/
    w=M_C*k*1e10;

    kx=xo;ky=yo;kz=dist;
    NORM(kx,ky,kz);
    kx*=k;
    ky*=k;
    kz*=k;


    double w1theta=2*prms.gamma2/(1+K*K/2.0 + prms.gamma2 *theta2) * w_u; 
    double w10=2*prms.gamma2/(1+K*K/2.0) * w_u; 
    double w_w1=w/w1theta;
    double w_w10=w/w10;

    double Bsig_integ_prms[4],Bpi_integ_prms[4];
    Bsig_integ_prms[0] = Bpi_integ_prms[0]=w_w1; /*relative frequency*/
    Bsig_integ_prms[1] = Bpi_integ_prms[1] = 2*w_w10*xi*prms.gamma*K/(1+K*K/2.0); /*p*/
    Bsig_integ_prms[2] = Bpi_integ_prms[2] = 0.25 * w_w10*K*K/(1+K*K/2.0); /*q*/
    Bsig_integ_prms[3] = xi *prms.gamma/K; /*angle terms*/
    Bpi_integ_prms[3]  = psi*prms.gamma/K; 

    /*pass the parameters to the integrand*/
    prms.Bsig.params=Bsig_integ_prms;
    prms.Bpi.params =Bpi_integ_prms;
    if (!quick_integ){
        gsl_integration_qags (&(prms.Bsig),    0, M_PI,    0, 1e-6, 1000, prms.gsl_int_ws, &bsigma_integral, &bsigma_error); 
        gsl_integration_qags (&(prms.Bpi),     0, M_PI,    0, 1e-6, 1000, prms.gsl_int_ws, &bpi_integral, &bpi_error); 
    }else{
        size_t neval;
        gsl_integration_qng (&(prms.Bsig),    0, M_PI,    0, 1e-6, &bsigma_integral, &bsigma_error,&neval);
        gsl_integration_qng (&(prms.Bpi),     0, M_PI,    0, 1e-6, &bpi_integral, &bpi_error,&neval); 
    }
    bsigma_integral*=M_2_PI;/*correct for only integrating half interval and normalize by pi.*/
    bpi_integral*=M_2_PI;

    s_n2=mxundulator_S_N(w_w1, Nper);
    
    double prefactor=alpha*Ie/CELE*pow(K*prms.gamma/(1+K*K/2.0),2.0)*Nper*Nper*w_w10*w_w10;

    p*=prefactor* ( pow(bsigma_integral,2.0)+pow(bpi_integral,2.0) ) *s_n2;

    /*randomly pick phase*/
    if (randomphase){
        phi=rand01()*2*M_PI;
    }else{
        phi=phase;
    }

    /*Set polarization vector. TODO: Do this right.*/
    Ex=0;Ey=0;Ez=0;
}
#line 10509 "Monochromator_testing.c"
}   /* End of dmu=Undulator() SETTING parameter declarations. */
#undef MELE
#undef alpha
#undef prms
#undef verbose
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompdmu:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Mbefore_E_monitor [4] */
  mccoordschange(mcposrMbefore_E_monitor, mcrotrMbefore_E_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Mbefore_E_monitor (without coords transformations) */
  mcJumpTrace_Mbefore_E_monitor:
  SIG_MESSAGE("Mbefore_E_monitor (Trace)");
  mcDEBUG_COMP("Mbefore_E_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMbefore_E_monitor
  STORE_XRAY(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  Mbefore_E_monitor
#define mccompcurtype  E_monitor
#define mccompcurindex 4
#define nE mccMbefore_E_monitor_nE
#define filename mccMbefore_E_monitor_filename
#define E_N mccMbefore_E_monitor_E_N
#define E_p mccMbefore_E_monitor_E_p
#define E_p2 mccMbefore_E_monitor_E_p2
{   /* Declarations of Mbefore_E_monitor=E_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_E_monitor_xmin;
MCNUM xmax = mccMbefore_E_monitor_xmax;
MCNUM ymin = mccMbefore_E_monitor_ymin;
MCNUM ymax = mccMbefore_E_monitor_ymax;
MCNUM xwidth = mccMbefore_E_monitor_xwidth;
MCNUM yheight = mccMbefore_E_monitor_yheight;
MCNUM Emin = mccMbefore_E_monitor_Emin;
MCNUM Emax = mccMbefore_E_monitor_Emax;
MCNUM restore_xray = mccMbefore_E_monitor_restore_xray;
#line 85 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;
    double E;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      E = K2E*sqrt(kx*kx + ky*ky + kz*kz);

      i = floor((E-Emin)*nE/(Emax-Emin));
      if(i >= 0 && i < nE)
      {
        E_N[i]++;
        E_p[i] += p;
        E_p2[i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 10662 "Monochromator_testing.c"
}   /* End of Mbefore_E_monitor=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompMbefore_E_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Mbefore_wl_monitor [5] */
  mccoordschange(mcposrMbefore_wl_monitor, mcrotrMbefore_wl_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Mbefore_wl_monitor (without coords transformations) */
  mcJumpTrace_Mbefore_wl_monitor:
  SIG_MESSAGE("Mbefore_wl_monitor (Trace)");
  mcDEBUG_COMP("Mbefore_wl_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMbefore_wl_monitor
  STORE_XRAY(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  Mbefore_wl_monitor
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccMbefore_wl_monitor_nL
#define filename mccMbefore_wl_monitor_filename
#define L_N mccMbefore_wl_monitor_L_N
#define L_p mccMbefore_wl_monitor_L_p
#define L_p2 mccMbefore_wl_monitor_L_p2
{   /* Declarations of Mbefore_wl_monitor=L_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_wl_monitor_xmin;
MCNUM xmax = mccMbefore_wl_monitor_xmax;
MCNUM ymin = mccMbefore_wl_monitor_ymin;
MCNUM ymax = mccMbefore_wl_monitor_ymax;
MCNUM xwidth = mccMbefore_wl_monitor_xwidth;
MCNUM yheight = mccMbefore_wl_monitor_yheight;
MCNUM Lmin = mccMbefore_wl_monitor_Lmin;
MCNUM Lmax = mccMbefore_wl_monitor_Lmax;
MCNUM restore_xray = mccMbefore_wl_monitor_restore_xray;
#line 84 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = 2*PI/sqrt(kx*kx + ky*ky + kz*kz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
      if(i >= 0 && i < nL)
      {
        L_N[i]++;
        L_p[i] += p;
        L_p2[i] += p*p;
        SCATTER;
      }
    } 
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 10815 "Monochromator_testing.c"
}   /* End of Mbefore_wl_monitor=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompMbefore_wl_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Mbefore_psd_monitor [6] */
  mccoordschange(mcposrMbefore_psd_monitor, mcrotrMbefore_psd_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Mbefore_psd_monitor (without coords transformations) */
  mcJumpTrace_Mbefore_psd_monitor:
  SIG_MESSAGE("Mbefore_psd_monitor (Trace)");
  mcDEBUG_COMP("Mbefore_psd_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMbefore_psd_monitor
  STORE_XRAY(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
{   /* Declarations of Mbefore_psd_monitor=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_psd_monitor_xmin;
MCNUM xmax = mccMbefore_psd_monitor_xmax;
MCNUM ymin = mccMbefore_psd_monitor_ymin;
MCNUM ymax = mccMbefore_psd_monitor_ymax;
MCNUM xwidth = mccMbefore_psd_monitor_xwidth;
MCNUM yheight = mccMbefore_psd_monitor_yheight;
MCNUM radius = mccMbefore_psd_monitor_radius;
#line 105 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;

    PROP_Z0;
    if (!radius){
      if (x>xmin && x<xmax && y>ymin && y<ymax)
      {
        i = floor((x - xmin)*nx/(xmax - xmin));
        j = floor((y - ymin)*ny/(ymax - ymin));
        PSD_N[i][j]++;
        PSD_p[i][j] += p;
        PSD_p2[i][j] += p*p;
        SCATTER;
      }
    }else{
      double r=sqrt(x*x+y*y);
      if (r<radius){
        i = floor(r*nr/radius);
        PSD_N[0][i]++;
        PSD_p[0][i] += p;
        PSD_p2[0][i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 10976 "Monochromator_testing.c"
}   /* End of Mbefore_psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompMbefore_psd_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component PG1_arm [7] */
  mccoordschange(mcposrPG1_arm, mcrotrPG1_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component PG1_arm (without coords transformations) */
  mcJumpTrace_PG1_arm:
  SIG_MESSAGE("PG1_arm (Trace)");
  mcDEBUG_COMP("PG1_arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPG1_arm
  STORE_XRAY(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  PG1_arm
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompPG1_arm:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component M2_rotation_arm1 [8] */
  mccoordschange(mcposrM2_rotation_arm1, mcrotrM2_rotation_arm1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component M2_rotation_arm1 (without coords transformations) */
  mcJumpTrace_M2_rotation_arm1:
  SIG_MESSAGE("M2_rotation_arm1 (Trace)");
  mcDEBUG_COMP("M2_rotation_arm1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompM2_rotation_arm1
  STORE_XRAY(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  M2_rotation_arm1
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompM2_rotation_arm1:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component M2_rotation_arm2 [9] */
  mccoordschange(mcposrM2_rotation_arm2, mcrotrM2_rotation_arm2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component M2_rotation_arm2 (without coords transformations) */
  mcJumpTrace_M2_rotation_arm2:
  SIG_MESSAGE("M2_rotation_arm2 (Trace)");
  mcDEBUG_COMP("M2_rotation_arm2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompM2_rotation_arm2
  STORE_XRAY(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  M2_rotation_arm2
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompM2_rotation_arm2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component M2_rotation_arm3 [10] */
  mccoordschange(mcposrM2_rotation_arm3, mcrotrM2_rotation_arm3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component M2_rotation_arm3 (without coords transformations) */
  mcJumpTrace_M2_rotation_arm3:
  SIG_MESSAGE("M2_rotation_arm3 (Trace)");
  mcDEBUG_COMP("M2_rotation_arm3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompM2_rotation_arm3
  STORE_XRAY(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  M2_rotation_arm3
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompM2_rotation_arm3:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component mirror2 [11] */
  mccoordschange(mcposrmirror2, mcrotrmirror2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component mirror2 (without coords transformations) */
  mcJumpTrace_mirror2:
  SIG_MESSAGE("mirror2 (Trace)");
  mcDEBUG_COMP("mirror2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompmirror2
  STORE_XRAY(11,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[11]++;
  mcPCounter[11] += p;
  mcP2Counter[11] += p*p;
#define mccompcurname  mirror2
#define mccompcurtype  Mirror
#define mccompcurindex 11
#define reflec mccmirror2_reflec
#define reflec_table mccmirror2_reflec_table
#define prms mccmirror2_prms
{   /* Declarations of mirror2=Mirror() SETTING parameters. */
MCNUM zdepth = mccmirror2_zdepth;
MCNUM xwidth = mccmirror2_xwidth;
MCNUM R0 = mccmirror2_R0;
#line 95 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Mirror.comp"
{
    int status;
    double l0,l1,l2,l3,tx,ty,tz;
    
    PROP_Y0;
    if(x<-xwidth/2.0|| x>xwidth/2.0 || z<-zdepth/2.0 || z>zdepth/2.0){
        RESTORE_XRAY(INDEX_CURRENT_COMP, x,y,z, kx,ky,kz, phi,t, Ex,Ey,Ez, p);
    }else{
        double nx,ny,nz;
        nx=0;
        ny=1;
        nz=0;

        double s=scalar_prod(kx,ky,kz,nx,ny,nz);
        double k=sqrt(scalar_prod(kx,ky,kz,kx,ky,kz)); 
         
        kx=kx-2*s*nx;
        ky=ky-2*s*ny;
        kz=kz-2*s*nz;

        /*find energy, and glancing angle*/ 
        if (prms.use_reflec_table){
            /*the get ref. by call to Table_value2d*/
            double R;
            double theta=RAD2DEG*(acos(s/k)-M_PI_2);
            double e=K2E*k;
            R=Table_Value2d(reflec_table,(e-prms.e_min)/prms.e_step, ((theta-prms.theta_min)/prms.theta_step));
            p*=R;
        }else{
            p*=R0;
        }
    }
}
#line 11572 "Monochromator_testing.c"
}   /* End of mirror2=Mirror() SETTING parameter declarations. */
#undef prms
#undef reflec_table
#undef reflec
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompmirror2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(11,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Plane_grating [12] */
  mccoordschange(mcposrPlane_grating, mcrotrPlane_grating,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Plane_grating (without coords transformations) */
  mcJumpTrace_Plane_grating:
  SIG_MESSAGE("Plane_grating (Trace)");
  mcDEBUG_COMP("Plane_grating")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPlane_grating
  STORE_XRAY(12,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[12]++;
  mcPCounter[12] += p;
  mcP2Counter[12] += p*p;
#define mccompcurname  Plane_grating
#define mccompcurtype  MultiPurpose_grating
#define mccompcurindex 12
#define blazed mccPlane_grating_blazed
#define display mccPlane_grating_display
#define zdepth mccPlane_grating_zdepth
#define xwidth mccPlane_grating_xwidth
{   /* Declarations of Plane_grating=MultiPurpose_grating() SETTING parameters. */
MCNUM MCangleVariation = mccPlane_grating_MCangleVariation;
MCNUM R0 = mccPlane_grating_R0;
MCNUM r_rho = mccPlane_grating_r_rho;
MCNUM b = mccPlane_grating_b;
MCNUM N_slits = mccPlane_grating_N_slits;
MCNUM d = mccPlane_grating_d;
MCNUM blazed_angle = mccPlane_grating_blazed_angle;
#line 140 "MultiPurpose_grating.comp"
{
/********************************************************************************************
Initializing by simulating a perfectly reflecting mirror:
********************************************************************************************/
    /*Placing the grating on the Y axis. */
     PROP_Y0;
    /*If the photon is passing the grating, it should not be reflected. Instead, for book-keeping, it is restored to previously state before the component using RESTORE_XRAY */
if(x<-xwidth/2.0|| x>xwidth/2.0 || z<-zdepth/2.0 || z>zdepth/2.0){
        RESTORE_XRAY(INDEX_CURRENT_COMP, x,y,z, kx,ky,kz, phi,t, Ex,Ey,Ez, p);
}else{
        double nx,ny,nz,theta_rand,Gamma,gamma,placeHolderIn,placeHolderOut,comp_angle,nlenght,DiffractionInterferenceFactor;
        /* Normal vector for the grating*/
        nx=0;
        ny=1;
        nz=0;
        // scalar-product, s and length of k.
        double s=scalar_prod(kx,ky,kz,nx,ny,nz);
        double k=sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));
        nlenght=sqrt(scalar_prod(nx,ny,nz,nx,ny,nz)); // Note, nlength should be 1.

        // outgoing vector for a perfectly reflecting mirror.
        kx=kx-2*s*nx;
        ky=ky-2*s*ny;
        kz=kz-2*s*nz;
        // Outgoing grazing angle for a perfectly reflecting mirror using angle between two vectors.
        comp_angle = acos(s/(k*nlenght))-M_PI_2;  
        // MC picked angle going from +- MCangleVariation
        theta_rand =rand01()*(2*MCangleVariation*DEG2RAD)-(MCangleVariation*DEG2RAD);
/********************************************************************************************
If the grating is blazed:
********************************************************************************************/
     if(blazed && blazed_angle){
              /*
               If there is a blazed angle. This will be used to find the outgoing angle based on
               angle_in + angle_out = 2*blazed_angle.
               For simplicity in the code, this is done in steps of two. 
               - Finding outgoing kx,ky,kz.
               */
               kx=kx;
               ky=(ky*cos(-2*blazed_angle*DEG2RAD)-kz*sin(-2*blazed_angle*DEG2RAD));
               kz=(ky*sin(-2*blazed_angle*DEG2RAD)+kz*cos(-2*blazed_angle*DEG2RAD));
               kx=kx;
               /* Adding the random angle. */
               ky=(ky*cos(theta_rand)-kz*sin(theta_rand));
               kz=(ky*sin(theta_rand)+kz*cos(theta_rand));
               k=sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));                             
     }
/********************************************************************************************
If the grating is lamellar:
********************************************************************************************/        
     if(!blazed){
        /* If grating is not blazed, lamellar grating is used.
         Thus finding outgoing angle only from ingoing angle.
        using old k and n and rotation matrix.*/
        kx=kx;
        ky=(ky*cos(theta_rand)-kz*sin(theta_rand));
        kz=(ky*sin(theta_rand)+kz*cos(theta_rand));
        k=sqrt(scalar_prod(kx,ky,kz,kx,ky,kz)); 
              if(display){      
              printf("Incidence grazing angle: %f deg. \n Reflecting grazing angle: %f deg \n",comp_angle*RAD2DEG,(theta_rand+comp_angle)*RAD2DEG);
              }  
     }
/********************************************************************************************
If photons are reflected into the grating:
********************************************************************************************/  
if(   ( ( (acos (s/ (k*nlenght))-M_PI_2)*RAD2DEG)  +(theta_rand*RAD2DEG))<0  ){
RESTORE_XRAY(INDEX_CURRENT_COMP, x,y,z, kx,ky,kz, phi,t, Ex,Ey,Ez, p);
}

/************************************************************************************************************
Finding the weight using diffraction theory. 
*************************************************************************************************************/               
     if (blazed && blazed_angle){
            placeHolderIn = sin(comp_angle+2*blazed_angle*DEG2RAD); 
            placeHolderOut = sin(comp_angle+2*blazed_angle*DEG2RAD+theta_rand); // asin angle out 
     } 
     if (!blazed){
            placeHolderIn = sin(comp_angle); 
            placeHolderOut = sin(comp_angle+theta_rand); // asin angle out
     }        
     /* Phase for interference pattern:  k = wave vector, d = lines/Angstrom, b=width in Aangsgrom */
     gamma = k*d*(placeHolderOut-placeHolderIn);
     /* Phase for diffraction pattern */
     Gamma = k*b*(placeHolderOut-placeHolderIn);
              if(display){      
              printf("\n Original weight= %e    \n with MC angle= %f , gamma = %f and  Gamma=%f. And DirectionalSamplingFactor=%f and DiffractionInterferenceFactor=%e,  the weight is updated. \n",p,theta_rand*RAD2DEG,gamma,Gamma,DirectionalSamplingFactor,DiffractionInterferenceFactor);
              }           
     DiffractionInterferenceFactor=((sin(Gamma/2)/(Gamma/2))*(sin(Gamma/2)/(Gamma/2)))*((sin(N_slits*gamma/2)/sin(gamma/2))*(sin(N_slits*gamma/2)/sin(gamma/2)));
     p=p*DirectionalSamplingFactor*(DiffractionInterferenceFactor/(pow(N_slits,2)));
              if(display){      
              printf("New weight = %e . \n",p);
              }  		
}

}
#line 11794 "Monochromator_testing.c"
}   /* End of Plane_grating=MultiPurpose_grating() SETTING parameter declarations. */
#undef xwidth
#undef zdepth
#undef display
#undef blazed
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompPlane_grating:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(12,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Monochromator_Monitor_arm [13] */
  mccoordschange(mcposrMonochromator_Monitor_arm, mcrotrMonochromator_Monitor_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Monochromator_Monitor_arm (without coords transformations) */
  mcJumpTrace_Monochromator_Monitor_arm:
  SIG_MESSAGE("Monochromator_Monitor_arm (Trace)");
  mcDEBUG_COMP("Monochromator_Monitor_arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMonochromator_Monitor_arm
  STORE_XRAY(13,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[13]++;
  mcPCounter[13] += p;
  mcP2Counter[13] += p*p;
#define mccompcurname  Monochromator_Monitor_arm
#define mccompcurtype  Arm
#define mccompcurindex 13
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompMonochromator_Monitor_arm:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(13,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Close_wl_monitor1 [14] */
  mccoordschange(mcposrClose_wl_monitor1, mcrotrClose_wl_monitor1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Close_wl_monitor1 (without coords transformations) */
  mcJumpTrace_Close_wl_monitor1:
  SIG_MESSAGE("Close_wl_monitor1 (Trace)");
  mcDEBUG_COMP("Close_wl_monitor1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompClose_wl_monitor1
  STORE_XRAY(14,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[14]++;
  mcPCounter[14] += p;
  mcP2Counter[14] += p*p;
#define mccompcurname  Close_wl_monitor1
#define mccompcurtype  L_monitor
#define mccompcurindex 14
#define nL mccClose_wl_monitor1_nL
#define filename mccClose_wl_monitor1_filename
#define L_N mccClose_wl_monitor1_L_N
#define L_p mccClose_wl_monitor1_L_p
#define L_p2 mccClose_wl_monitor1_L_p2
{   /* Declarations of Close_wl_monitor1=L_monitor() SETTING parameters. */
MCNUM xmin = mccClose_wl_monitor1_xmin;
MCNUM xmax = mccClose_wl_monitor1_xmax;
MCNUM ymin = mccClose_wl_monitor1_ymin;
MCNUM ymax = mccClose_wl_monitor1_ymax;
MCNUM xwidth = mccClose_wl_monitor1_xwidth;
MCNUM yheight = mccClose_wl_monitor1_yheight;
MCNUM Lmin = mccClose_wl_monitor1_Lmin;
MCNUM Lmax = mccClose_wl_monitor1_Lmax;
MCNUM restore_xray = mccClose_wl_monitor1_restore_xray;
#line 84 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = 2*PI/sqrt(kx*kx + ky*ky + kz*kz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
      if(i >= 0 && i < nL)
      {
        L_N[i]++;
        L_p[i] += p;
        L_p2[i] += p*p;
        SCATTER;
      }
    } 
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12055 "Monochromator_testing.c"
}   /* End of Close_wl_monitor1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompClose_wl_monitor1:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(14,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component close_E_monitor1 [15] */
  mccoordschange(mcposrclose_E_monitor1, mcrotrclose_E_monitor1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component close_E_monitor1 (without coords transformations) */
  mcJumpTrace_close_E_monitor1:
  SIG_MESSAGE("close_E_monitor1 (Trace)");
  mcDEBUG_COMP("close_E_monitor1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompclose_E_monitor1
  STORE_XRAY(15,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[15]++;
  mcPCounter[15] += p;
  mcP2Counter[15] += p*p;
#define mccompcurname  close_E_monitor1
#define mccompcurtype  E_monitor
#define mccompcurindex 15
#define nE mccclose_E_monitor1_nE
#define filename mccclose_E_monitor1_filename
#define E_N mccclose_E_monitor1_E_N
#define E_p mccclose_E_monitor1_E_p
#define E_p2 mccclose_E_monitor1_E_p2
{   /* Declarations of close_E_monitor1=E_monitor() SETTING parameters. */
MCNUM xmin = mccclose_E_monitor1_xmin;
MCNUM xmax = mccclose_E_monitor1_xmax;
MCNUM ymin = mccclose_E_monitor1_ymin;
MCNUM ymax = mccclose_E_monitor1_ymax;
MCNUM xwidth = mccclose_E_monitor1_xwidth;
MCNUM yheight = mccclose_E_monitor1_yheight;
MCNUM Emin = mccclose_E_monitor1_Emin;
MCNUM Emax = mccclose_E_monitor1_Emax;
MCNUM restore_xray = mccclose_E_monitor1_restore_xray;
#line 85 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;
    double E;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      E = K2E*sqrt(kx*kx + ky*ky + kz*kz);

      i = floor((E-Emin)*nE/(Emax-Emin));
      if(i >= 0 && i < nE)
      {
        E_N[i]++;
        E_p[i] += p;
        E_p2[i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12209 "Monochromator_testing.c"
}   /* End of close_E_monitor1=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompclose_E_monitor1:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(15,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component close_psd_monitor2 [16] */
  mccoordschange(mcposrclose_psd_monitor2, mcrotrclose_psd_monitor2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component close_psd_monitor2 (without coords transformations) */
  mcJumpTrace_close_psd_monitor2:
  SIG_MESSAGE("close_psd_monitor2 (Trace)");
  mcDEBUG_COMP("close_psd_monitor2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompclose_psd_monitor2
  STORE_XRAY(16,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[16]++;
  mcPCounter[16] += p;
  mcP2Counter[16] += p*p;
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
{   /* Declarations of close_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccclose_psd_monitor2_xmin;
MCNUM xmax = mccclose_psd_monitor2_xmax;
MCNUM ymin = mccclose_psd_monitor2_ymin;
MCNUM ymax = mccclose_psd_monitor2_ymax;
MCNUM xwidth = mccclose_psd_monitor2_xwidth;
MCNUM yheight = mccclose_psd_monitor2_yheight;
MCNUM radius = mccclose_psd_monitor2_radius;
#line 105 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;

    PROP_Z0;
    if (!radius){
      if (x>xmin && x<xmax && y>ymin && y<ymax)
      {
        i = floor((x - xmin)*nx/(xmax - xmin));
        j = floor((y - ymin)*ny/(ymax - ymin));
        PSD_N[i][j]++;
        PSD_p[i][j] += p;
        PSD_p2[i][j] += p*p;
        SCATTER;
      }
    }else{
      double r=sqrt(x*x+y*y);
      if (r<radius){
        i = floor(r*nr/radius);
        PSD_N[0][i]++;
        PSD_p[0][i] += p;
        PSD_p2[0][i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12370 "Monochromator_testing.c"
}   /* End of close_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompclose_psd_monitor2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(16,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Far_E_monitor2 [17] */
  mccoordschange(mcposrFar_E_monitor2, mcrotrFar_E_monitor2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Far_E_monitor2 (without coords transformations) */
  mcJumpTrace_Far_E_monitor2:
  SIG_MESSAGE("Far_E_monitor2 (Trace)");
  mcDEBUG_COMP("Far_E_monitor2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFar_E_monitor2
  STORE_XRAY(17,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[17]++;
  mcPCounter[17] += p;
  mcP2Counter[17] += p*p;
#define mccompcurname  Far_E_monitor2
#define mccompcurtype  E_monitor
#define mccompcurindex 17
#define nE mccFar_E_monitor2_nE
#define filename mccFar_E_monitor2_filename
#define E_N mccFar_E_monitor2_E_N
#define E_p mccFar_E_monitor2_E_p
#define E_p2 mccFar_E_monitor2_E_p2
{   /* Declarations of Far_E_monitor2=E_monitor() SETTING parameters. */
MCNUM xmin = mccFar_E_monitor2_xmin;
MCNUM xmax = mccFar_E_monitor2_xmax;
MCNUM ymin = mccFar_E_monitor2_ymin;
MCNUM ymax = mccFar_E_monitor2_ymax;
MCNUM xwidth = mccFar_E_monitor2_xwidth;
MCNUM yheight = mccFar_E_monitor2_yheight;
MCNUM Emin = mccFar_E_monitor2_Emin;
MCNUM Emax = mccFar_E_monitor2_Emax;
MCNUM restore_xray = mccFar_E_monitor2_restore_xray;
#line 85 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;
    double E;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      E = K2E*sqrt(kx*kx + ky*ky + kz*kz);

      i = floor((E-Emin)*nE/(Emax-Emin));
      if(i >= 0 && i < nE)
      {
        E_N[i]++;
        E_p[i] += p;
        E_p2[i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12527 "Monochromator_testing.c"
}   /* End of Far_E_monitor2=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompFar_E_monitor2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(17,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Fardivergence_monitor [18] */
  mccoordschange(mcposrFardivergence_monitor, mcrotrFardivergence_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Fardivergence_monitor (without coords transformations) */
  mcJumpTrace_Fardivergence_monitor:
  SIG_MESSAGE("Fardivergence_monitor (Trace)");
  mcDEBUG_COMP("Fardivergence_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFardivergence_monitor
  STORE_XRAY(18,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[18]++;
  mcPCounter[18] += p;
  mcP2Counter[18] += p*p;
#define mccompcurname  Fardivergence_monitor
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 18
#define nh mccFardivergence_monitor_nh
#define nv mccFardivergence_monitor_nv
#define rad mccFardivergence_monitor_rad
#define Div_N mccFardivergence_monitor_Div_N
#define Div_p mccFardivergence_monitor_Div_p
#define Div_p2 mccFardivergence_monitor_Div_p2
{   /* Declarations of Fardivergence_monitor=Divergence_monitor() SETTING parameters. */
char* filename = mccFardivergence_monitor_filename;
MCNUM xmin = mccFardivergence_monitor_xmin;
MCNUM xmax = mccFardivergence_monitor_xmax;
MCNUM ymin = mccFardivergence_monitor_ymin;
MCNUM ymax = mccFardivergence_monitor_ymax;
MCNUM xwidth = mccFardivergence_monitor_xwidth;
MCNUM yheight = mccFardivergence_monitor_yheight;
MCNUM maxdiv_h = mccFardivergence_monitor_maxdiv_h;
MCNUM maxdiv_v = mccFardivergence_monitor_maxdiv_v;
MCNUM restore_xray = mccFardivergence_monitor_restore_xray;
MCNUM nx = mccFardivergence_monitor_nx;
MCNUM ny = mccFardivergence_monitor_ny;
MCNUM nz = mccFardivergence_monitor_nz;
#line 92 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
    int i,j;
    double h_div, v_div;
    double k, kn;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      /* Find length of projection onto the [nx ny nz] axis */
      kn = scalar_prod(kx, ky, kz, nx, ny, nz);
      if (rad){
          h_div = atan2(kx,kn);
          v_div = atan2(ky,kn);
      }else{
          h_div = RAD2DEG*atan2(kx,kn);
          v_div = RAD2DEG*atan2(ky,kn);
      }
      if (h_div < maxdiv_h && h_div > -maxdiv_h &&
          v_div < maxdiv_v && v_div > -maxdiv_v)
      {
        i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
        j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
        Div_N[i][j]++;
        Div_p[i][j] += p;
        Div_p2[i][j] += p*p;
        SCATTER;
      }
    }   
    if (restore_xray) {
        RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12696 "Monochromator_testing.c"
}   /* End of Fardivergence_monitor=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef rad
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompFardivergence_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(18,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Far_wl_monitor2 [19] */
  mccoordschange(mcposrFar_wl_monitor2, mcrotrFar_wl_monitor2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Far_wl_monitor2 (without coords transformations) */
  mcJumpTrace_Far_wl_monitor2:
  SIG_MESSAGE("Far_wl_monitor2 (Trace)");
  mcDEBUG_COMP("Far_wl_monitor2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFar_wl_monitor2
  STORE_XRAY(19,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[19]++;
  mcPCounter[19] += p;
  mcP2Counter[19] += p*p;
#define mccompcurname  Far_wl_monitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mccFar_wl_monitor2_nL
#define filename mccFar_wl_monitor2_filename
#define L_N mccFar_wl_monitor2_L_N
#define L_p mccFar_wl_monitor2_L_p
#define L_p2 mccFar_wl_monitor2_L_p2
{   /* Declarations of Far_wl_monitor2=L_monitor() SETTING parameters. */
MCNUM xmin = mccFar_wl_monitor2_xmin;
MCNUM xmax = mccFar_wl_monitor2_xmax;
MCNUM ymin = mccFar_wl_monitor2_ymin;
MCNUM ymax = mccFar_wl_monitor2_ymax;
MCNUM xwidth = mccFar_wl_monitor2_xwidth;
MCNUM yheight = mccFar_wl_monitor2_yheight;
MCNUM Lmin = mccFar_wl_monitor2_Lmin;
MCNUM Lmax = mccFar_wl_monitor2_Lmax;
MCNUM restore_xray = mccFar_wl_monitor2_restore_xray;
#line 84 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = 2*PI/sqrt(kx*kx + ky*ky + kz*kz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
      if(i >= 0 && i < nL)
      {
        L_N[i]++;
        L_p[i] += p;
        L_p2[i] += p*p;
        SCATTER;
      }
    } 
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 12850 "Monochromator_testing.c"
}   /* End of Far_wl_monitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompFar_wl_monitor2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(19,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  /* TRACE Component Far_psd_monitor2 [20] */
  mccoordschange(mcposrFar_psd_monitor2, mcrotrFar_psd_monitor2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlkx,
    &mcnlky,
    &mcnlkz,
    &mcnlEx,
    &mcnlEy,
    &mcnlEz);
  /* define label inside component Far_psd_monitor2 (without coords transformations) */
  mcJumpTrace_Far_psd_monitor2:
  SIG_MESSAGE("Far_psd_monitor2 (Trace)");
  mcDEBUG_COMP("Far_psd_monitor2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define kx mcnlkx
#define ky mcnlky
#define kz mcnlkz
#define phi mcnlphi
#define t mcnlt
#define Ex mcnlEx
#define Ey mcnlEy
#define Ez mcnlEz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFar_psd_monitor2
  STORE_XRAY(20,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlkx,
    mcnlky,
    mcnlkz,
    mcnlphi,
    mcnlt,
    mcnlEx,
    mcnlEy,
    mcnlEz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[20]++;
  mcPCounter[20] += p;
  mcP2Counter[20] += p*p;
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
{   /* Declarations of Far_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccFar_psd_monitor2_xmin;
MCNUM xmax = mccFar_psd_monitor2_xmax;
MCNUM ymin = mccFar_psd_monitor2_ymin;
MCNUM ymax = mccFar_psd_monitor2_ymax;
MCNUM xwidth = mccFar_psd_monitor2_xwidth;
MCNUM yheight = mccFar_psd_monitor2_yheight;
MCNUM radius = mccFar_psd_monitor2_radius;
#line 105 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    int i,j;

    PROP_Z0;
    if (!radius){
      if (x>xmin && x<xmax && y>ymin && y<ymax)
      {
        i = floor((x - xmin)*nx/(xmax - xmin));
        j = floor((y - ymin)*ny/(ymax - ymin));
        PSD_N[i][j]++;
        PSD_p[i][j] += p;
        PSD_p2[i][j] += p*p;
        SCATTER;
      }
    }else{
      double r=sqrt(x*x+y*y);
      if (r<radius){
        i = floor(r*nr/radius);
        PSD_N[0][i]++;
        PSD_p[0][i] += p;
        PSD_p2[0][i] += p*p;
        SCATTER;
      }
    }
    if (restore_xray) {
      RESTORE_XRAY(INDEX_CURRENT_COMP, x, y, z, kx, ky, kz, phi, t, Ex, Ey, Ez, p);
    }
}
#line 13011 "Monochromator_testing.c"
}   /* End of Far_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  xray */
  mcabsorbCompFar_psd_monitor2:
  if (RESTORE) /* restore if needed */
  { RESTORE_XRAY(20,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlkx,
      mcnlky,
      mcnlkz,
      mcnlphi,
      mcnlt,
      mcnlEx,
      mcnlEy,
      mcnlEz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef Ez
#undef Ey
#undef Ex
#undef t
#undef phi
#undef kz
#undef ky
#undef kx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlkx,
mcnlky,
mcnlkz,
mcnlphi,
mcnlt,
mcnlEx,
mcnlEy,
mcnlEz,
mcnlp)
  /* Copy xray state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnkx = mcnlkx;
  mcnky = mcnlky;
  mcnkz = mcnlkz;
  mcnphi = mcnlphi;
  mcnt = mcnlt;
  mcnEx = mcnlEx;
  mcnEy = mcnlEy;
  mcnEz = mcnlEz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'origin'. */
  SIG_MESSAGE("origin (Save)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define profile mccorigin_profile
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 110 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile)) {
    char filename[256];
    if (!strlen(profile)) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 13132 "Monochromator_testing.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef profile
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Mbefore_E_monitor'. */
  SIG_MESSAGE("Mbefore_E_monitor (Save)");
#define mccompcurname  Mbefore_E_monitor
#define mccompcurtype  E_monitor
#define mccompcurindex 4
#define nE mccMbefore_E_monitor_nE
#define filename mccMbefore_E_monitor_filename
#define E_N mccMbefore_E_monitor_E_N
#define E_p mccMbefore_E_monitor_E_p
#define E_p2 mccMbefore_E_monitor_E_p2
{   /* Declarations of Mbefore_E_monitor=E_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_E_monitor_xmin;
MCNUM xmax = mccMbefore_E_monitor_xmax;
MCNUM ymin = mccMbefore_E_monitor_ymin;
MCNUM ymax = mccMbefore_E_monitor_ymax;
MCNUM xwidth = mccMbefore_E_monitor_xwidth;
MCNUM yheight = mccMbefore_E_monitor_yheight;
MCNUM Emin = mccMbefore_E_monitor_Emin;
MCNUM Emax = mccMbefore_E_monitor_Emax;
MCNUM restore_xray = mccMbefore_E_monitor_restore_xray;
#line 108 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [keV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
}
#line 13172 "Monochromator_testing.c"
}   /* End of Mbefore_E_monitor=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Mbefore_wl_monitor'. */
  SIG_MESSAGE("Mbefore_wl_monitor (Save)");
#define mccompcurname  Mbefore_wl_monitor
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccMbefore_wl_monitor_nL
#define filename mccMbefore_wl_monitor_filename
#define L_N mccMbefore_wl_monitor_L_N
#define L_p mccMbefore_wl_monitor_L_p
#define L_p2 mccMbefore_wl_monitor_L_p2
{   /* Declarations of Mbefore_wl_monitor=L_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_wl_monitor_xmin;
MCNUM xmax = mccMbefore_wl_monitor_xmax;
MCNUM ymin = mccMbefore_wl_monitor_ymin;
MCNUM ymax = mccMbefore_wl_monitor_ymax;
MCNUM xwidth = mccMbefore_wl_monitor_xwidth;
MCNUM yheight = mccMbefore_wl_monitor_yheight;
MCNUM Lmin = mccMbefore_wl_monitor_Lmin;
MCNUM Lmax = mccMbefore_wl_monitor_Lmax;
MCNUM restore_xray = mccMbefore_wl_monitor_restore_xray;
#line 106 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
}
#line 13213 "Monochromator_testing.c"
}   /* End of Mbefore_wl_monitor=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Mbefore_psd_monitor'. */
  SIG_MESSAGE("Mbefore_psd_monitor (Save)");
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
{   /* Declarations of Mbefore_psd_monitor=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_psd_monitor_xmin;
MCNUM xmax = mccMbefore_psd_monitor_xmax;
MCNUM ymin = mccMbefore_psd_monitor_ymin;
MCNUM ymax = mccMbefore_psd_monitor_ymax;
MCNUM xwidth = mccMbefore_psd_monitor_xwidth;
MCNUM yheight = mccMbefore_psd_monitor_yheight;
MCNUM radius = mccMbefore_psd_monitor_radius;
#line 134 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    if(!radius){
        if(nx==1 && ny==1){
            DETECTOR_OUT_0D("Intensity monitor " NAME_CURRENT_COMP, (double) PSD_N[0][0], PSD_p[0][0], PSD_p2[0][0]);
        }else if(nx==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","Y Position[m]", "Intensity", "Y",
                    ymin,ymax,ny,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else if (ny==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","X Position[m]", "Intensity", "X",
                    xmin,xmax,nx,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else{
            DETECTOR_OUT_2D(
                    "PSD monitor",
                    "X position [m]",
                    "Y position [m]",
                    xmin, xmax, ymin, ymax,
                    nx, ny,
                    *PSD_N,*PSD_p,*PSD_p2,
                    filename);
        }
    }else{
      DETECTOR_OUT_1D(
          "PSD_monitor","Radial Position[m]", "Intensity", "R",
          0,radius,nr,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
    }

}
#line 13275 "Monochromator_testing.c"
}   /* End of Mbefore_psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Close_wl_monitor1'. */
  SIG_MESSAGE("Close_wl_monitor1 (Save)");
#define mccompcurname  Close_wl_monitor1
#define mccompcurtype  L_monitor
#define mccompcurindex 14
#define nL mccClose_wl_monitor1_nL
#define filename mccClose_wl_monitor1_filename
#define L_N mccClose_wl_monitor1_L_N
#define L_p mccClose_wl_monitor1_L_p
#define L_p2 mccClose_wl_monitor1_L_p2
{   /* Declarations of Close_wl_monitor1=L_monitor() SETTING parameters. */
MCNUM xmin = mccClose_wl_monitor1_xmin;
MCNUM xmax = mccClose_wl_monitor1_xmax;
MCNUM ymin = mccClose_wl_monitor1_ymin;
MCNUM ymax = mccClose_wl_monitor1_ymax;
MCNUM xwidth = mccClose_wl_monitor1_xwidth;
MCNUM yheight = mccClose_wl_monitor1_yheight;
MCNUM Lmin = mccClose_wl_monitor1_Lmin;
MCNUM Lmax = mccClose_wl_monitor1_Lmax;
MCNUM restore_xray = mccClose_wl_monitor1_restore_xray;
#line 106 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
}
#line 13319 "Monochromator_testing.c"
}   /* End of Close_wl_monitor1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'close_E_monitor1'. */
  SIG_MESSAGE("close_E_monitor1 (Save)");
#define mccompcurname  close_E_monitor1
#define mccompcurtype  E_monitor
#define mccompcurindex 15
#define nE mccclose_E_monitor1_nE
#define filename mccclose_E_monitor1_filename
#define E_N mccclose_E_monitor1_E_N
#define E_p mccclose_E_monitor1_E_p
#define E_p2 mccclose_E_monitor1_E_p2
{   /* Declarations of close_E_monitor1=E_monitor() SETTING parameters. */
MCNUM xmin = mccclose_E_monitor1_xmin;
MCNUM xmax = mccclose_E_monitor1_xmax;
MCNUM ymin = mccclose_E_monitor1_ymin;
MCNUM ymax = mccclose_E_monitor1_ymax;
MCNUM xwidth = mccclose_E_monitor1_xwidth;
MCNUM yheight = mccclose_E_monitor1_yheight;
MCNUM Emin = mccclose_E_monitor1_Emin;
MCNUM Emax = mccclose_E_monitor1_Emax;
MCNUM restore_xray = mccclose_E_monitor1_restore_xray;
#line 108 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [keV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
}
#line 13360 "Monochromator_testing.c"
}   /* End of close_E_monitor1=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'close_psd_monitor2'. */
  SIG_MESSAGE("close_psd_monitor2 (Save)");
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
{   /* Declarations of close_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccclose_psd_monitor2_xmin;
MCNUM xmax = mccclose_psd_monitor2_xmax;
MCNUM ymin = mccclose_psd_monitor2_ymin;
MCNUM ymax = mccclose_psd_monitor2_ymax;
MCNUM xwidth = mccclose_psd_monitor2_xwidth;
MCNUM yheight = mccclose_psd_monitor2_yheight;
MCNUM radius = mccclose_psd_monitor2_radius;
#line 134 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    if(!radius){
        if(nx==1 && ny==1){
            DETECTOR_OUT_0D("Intensity monitor " NAME_CURRENT_COMP, (double) PSD_N[0][0], PSD_p[0][0], PSD_p2[0][0]);
        }else if(nx==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","Y Position[m]", "Intensity", "Y",
                    ymin,ymax,ny,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else if (ny==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","X Position[m]", "Intensity", "X",
                    xmin,xmax,nx,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else{
            DETECTOR_OUT_2D(
                    "PSD monitor",
                    "X position [m]",
                    "Y position [m]",
                    xmin, xmax, ymin, ymax,
                    nx, ny,
                    *PSD_N,*PSD_p,*PSD_p2,
                    filename);
        }
    }else{
      DETECTOR_OUT_1D(
          "PSD_monitor","Radial Position[m]", "Intensity", "R",
          0,radius,nr,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
    }

}
#line 13422 "Monochromator_testing.c"
}   /* End of close_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Far_E_monitor2'. */
  SIG_MESSAGE("Far_E_monitor2 (Save)");
#define mccompcurname  Far_E_monitor2
#define mccompcurtype  E_monitor
#define mccompcurindex 17
#define nE mccFar_E_monitor2_nE
#define filename mccFar_E_monitor2_filename
#define E_N mccFar_E_monitor2_E_N
#define E_p mccFar_E_monitor2_E_p
#define E_p2 mccFar_E_monitor2_E_p2
{   /* Declarations of Far_E_monitor2=E_monitor() SETTING parameters. */
MCNUM xmin = mccFar_E_monitor2_xmin;
MCNUM xmax = mccFar_E_monitor2_xmax;
MCNUM ymin = mccFar_E_monitor2_ymin;
MCNUM ymax = mccFar_E_monitor2_ymax;
MCNUM xwidth = mccFar_E_monitor2_xwidth;
MCNUM yheight = mccFar_E_monitor2_yheight;
MCNUM Emin = mccFar_E_monitor2_Emin;
MCNUM Emax = mccFar_E_monitor2_Emax;
MCNUM restore_xray = mccFar_E_monitor2_restore_xray;
#line 108 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [keV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
}
#line 13466 "Monochromator_testing.c"
}   /* End of Far_E_monitor2=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Fardivergence_monitor'. */
  SIG_MESSAGE("Fardivergence_monitor (Save)");
#define mccompcurname  Fardivergence_monitor
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 18
#define nh mccFardivergence_monitor_nh
#define nv mccFardivergence_monitor_nv
#define rad mccFardivergence_monitor_rad
#define Div_N mccFardivergence_monitor_Div_N
#define Div_p mccFardivergence_monitor_Div_p
#define Div_p2 mccFardivergence_monitor_Div_p2
{   /* Declarations of Fardivergence_monitor=Divergence_monitor() SETTING parameters. */
char* filename = mccFardivergence_monitor_filename;
MCNUM xmin = mccFardivergence_monitor_xmin;
MCNUM xmax = mccFardivergence_monitor_xmax;
MCNUM ymin = mccFardivergence_monitor_ymin;
MCNUM ymax = mccFardivergence_monitor_ymax;
MCNUM xwidth = mccFardivergence_monitor_xwidth;
MCNUM yheight = mccFardivergence_monitor_yheight;
MCNUM maxdiv_h = mccFardivergence_monitor_maxdiv_h;
MCNUM maxdiv_v = mccFardivergence_monitor_maxdiv_v;
MCNUM restore_xray = mccFardivergence_monitor_restore_xray;
MCNUM nx = mccFardivergence_monitor_nx;
MCNUM ny = mccFardivergence_monitor_ny;
MCNUM nz = mccFardivergence_monitor_nz;
#line 125 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
      if (rad){
          DETECTOR_OUT_2D(
                  "Divergence monitor",
                  "X divergence [rad]",
                  "Y divergence [rad]",
                  -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
                  nh, nv,
                  &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
                  filename);
      } else {
          DETECTOR_OUT_2D(
                  "Divergence monitor",
                  "X divergence [deg]",
                  "Y divergence [deg]",
                  -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
                  nh, nv,
                  &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
                  filename);
      }
}
#line 13524 "Monochromator_testing.c"
}   /* End of Fardivergence_monitor=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef rad
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Far_wl_monitor2'. */
  SIG_MESSAGE("Far_wl_monitor2 (Save)");
#define mccompcurname  Far_wl_monitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mccFar_wl_monitor2_nL
#define filename mccFar_wl_monitor2_filename
#define L_N mccFar_wl_monitor2_L_N
#define L_p mccFar_wl_monitor2_L_p
#define L_p2 mccFar_wl_monitor2_L_p2
{   /* Declarations of Far_wl_monitor2=L_monitor() SETTING parameters. */
MCNUM xmin = mccFar_wl_monitor2_xmin;
MCNUM xmax = mccFar_wl_monitor2_xmax;
MCNUM ymin = mccFar_wl_monitor2_ymin;
MCNUM ymax = mccFar_wl_monitor2_ymax;
MCNUM xwidth = mccFar_wl_monitor2_xwidth;
MCNUM yheight = mccFar_wl_monitor2_yheight;
MCNUM Lmin = mccFar_wl_monitor2_Lmin;
MCNUM Lmax = mccFar_wl_monitor2_Lmax;
MCNUM restore_xray = mccFar_wl_monitor2_restore_xray;
#line 106 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
}
#line 13566 "Monochromator_testing.c"
}   /* End of Far_wl_monitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Far_psd_monitor2'. */
  SIG_MESSAGE("Far_psd_monitor2 (Save)");
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
{   /* Declarations of Far_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccFar_psd_monitor2_xmin;
MCNUM xmax = mccFar_psd_monitor2_xmax;
MCNUM ymin = mccFar_psd_monitor2_ymin;
MCNUM ymax = mccFar_psd_monitor2_ymax;
MCNUM xwidth = mccFar_psd_monitor2_xwidth;
MCNUM yheight = mccFar_psd_monitor2_yheight;
MCNUM radius = mccFar_psd_monitor2_radius;
#line 134 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
    if(!radius){
        if(nx==1 && ny==1){
            DETECTOR_OUT_0D("Intensity monitor " NAME_CURRENT_COMP, (double) PSD_N[0][0], PSD_p[0][0], PSD_p2[0][0]);
        }else if(nx==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","Y Position[m]", "Intensity", "Y",
                    ymin,ymax,ny,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else if (ny==1){
            DETECTOR_OUT_1D(
                    "PSD_monitor","X Position[m]", "Intensity", "X",
                    xmin,xmax,nx,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
        }else{
            DETECTOR_OUT_2D(
                    "PSD monitor",
                    "X position [m]",
                    "Y position [m]",
                    xmin, xmax, ymin, ymax,
                    nx, ny,
                    *PSD_N,*PSD_p,*PSD_p2,
                    filename);
        }
    }else{
      DETECTOR_OUT_1D(
          "PSD_monitor","Radial Position[m]", "Intensity", "R",
          0,radius,nr,&PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],filename);
    }

}
#line 13628 "Monochromator_testing.c"
}   /* End of Far_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'origin'. */
  SIG_MESSAGE("origin (Finally)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define profile mccorigin_profile
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 128 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 13675 "Monochromator_testing.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef profile
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No xray could reach Component[1] origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'source_flat'. */
  SIG_MESSAGE("source_flat (Finally)");
#define mccompcurname  source_flat
#define mccompcurtype  Source_flat
#define mccompcurindex 2
#define spectrum_file mccsource_flat_spectrum_file
#define prms mccsource_flat_prms
#define square mccsource_flat_square
{   /* Declarations of source_flat=Source_flat() SETTING parameters. */
MCNUM radius = mccsource_flat_radius;
MCNUM yheight = mccsource_flat_yheight;
MCNUM xwidth = mccsource_flat_xwidth;
MCNUM xmin = mccsource_flat_xmin;
MCNUM xmax = mccsource_flat_xmax;
MCNUM dist = mccsource_flat_dist;
MCNUM focus_xw = mccsource_flat_focus_xw;
MCNUM focus_yh = mccsource_flat_focus_yh;
MCNUM E0 = mccsource_flat_E0;
MCNUM dE = mccsource_flat_dE;
MCNUM lambda0 = mccsource_flat_lambda0;
MCNUM dlambda = mccsource_flat_dlambda;
MCNUM flux = mccsource_flat_flux;
MCNUM gauss = mccsource_flat_gauss;
MCNUM randomphase = mccsource_flat_randomphase;
MCNUM phase = mccsource_flat_phase;
#line 208 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
{
  Table_Free(&(prms.T));
}
#line 13716 "Monochromator_testing.c"
}   /* End of source_flat=Source_flat() SETTING parameter declarations. */
#undef square
#undef prms
#undef spectrum_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No xray could reach Component[2] source_flat\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source_flat=Source_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'dmu'. */
  SIG_MESSAGE("dmu (Finally)");
#define mccompcurname  dmu
#define mccompcurtype  Undulator
#define mccompcurindex 3
#define verbose mccdmu_verbose
#define prms mccdmu_prms
#define alpha mccdmu_alpha
#define MELE mccdmu_MELE
{   /* Declarations of dmu=Undulator() SETTING parameters. */
MCNUM E0 = mccdmu_E0;
MCNUM dE = mccdmu_dE;
MCNUM phase = mccdmu_phase;
MCNUM randomphase = mccdmu_randomphase;
MCNUM Ee = mccdmu_Ee;
MCNUM dEe = mccdmu_dEe;
MCNUM Ie = mccdmu_Ie;
MCNUM tbunch = mccdmu_tbunch;
MCNUM t0 = mccdmu_t0;
MCNUM B = mccdmu_B;
MCNUM K = mccdmu_K;
MCNUM gap = mccdmu_gap;
int Nper = mccdmu_Nper;
MCNUM lu = mccdmu_lu;
MCNUM sigey = mccdmu_sigey;
MCNUM sigex = mccdmu_sigex;
MCNUM sigepx = mccdmu_sigepx;
MCNUM sigepy = mccdmu_sigepy;
MCNUM focus_xw = mccdmu_focus_xw;
MCNUM focus_yh = mccdmu_focus_yh;
MCNUM dist = mccdmu_dist;
MCNUM gauss_t = mccdmu_gauss_t;
MCNUM quick_integ = mccdmu_quick_integ;
MCNUM E1st = mccdmu_E1st;
#line 308 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"
{
    if(!quick_integ){
        gsl_integration_workspace_free (prms.gsl_int_ws);
    }
}
#line 13767 "Monochromator_testing.c"
}   /* End of dmu=Undulator() SETTING parameter declarations. */
#undef MELE
#undef alpha
#undef prms
#undef verbose
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No xray could reach Component[3] dmu\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] dmu=Undulator()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No xray could reach Component[4] Mbefore_E_monitor\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] Mbefore_E_monitor=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No xray could reach Component[5] Mbefore_wl_monitor\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] Mbefore_wl_monitor=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'Mbefore_psd_monitor'. */
  SIG_MESSAGE("Mbefore_psd_monitor (Finally)");
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
{   /* Declarations of Mbefore_psd_monitor=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_psd_monitor_xmin;
MCNUM xmax = mccMbefore_psd_monitor_xmax;
MCNUM ymin = mccMbefore_psd_monitor_ymin;
MCNUM ymax = mccMbefore_psd_monitor_ymax;
MCNUM xwidth = mccMbefore_psd_monitor_xwidth;
MCNUM yheight = mccMbefore_psd_monitor_yheight;
MCNUM radius = mccMbefore_psd_monitor_radius;
#line 165 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
   free(PSD_N[0]);
   free(PSD_N);
   free(PSD_p[0]);
   free(PSD_p);
   free(PSD_p2[0]);
   free(PSD_p2);
}
#line 13813 "Monochromator_testing.c"
}   /* End of Mbefore_psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No xray could reach Component[6] Mbefore_psd_monitor\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Mbefore_psd_monitor=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No xray could reach Component[7] PG1_arm\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] PG1_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No xray could reach Component[8] M2_rotation_arm1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] M2_rotation_arm1=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No xray could reach Component[9] M2_rotation_arm2\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] M2_rotation_arm2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No xray could reach Component[10] M2_rotation_arm3\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] M2_rotation_arm3=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No xray could reach Component[11] mirror2\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] mirror2=Mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No xray could reach Component[12] Plane_grating\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] Plane_grating=MultiPurpose_grating()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No xray could reach Component[13] Monochromator_Monitor_arm\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] Monochromator_Monitor_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No xray could reach Component[14] Close_wl_monitor1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] Close_wl_monitor1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No xray could reach Component[15] close_E_monitor1\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] close_E_monitor1=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
  /* User FINALLY code for component 'close_psd_monitor2'. */
  SIG_MESSAGE("close_psd_monitor2 (Finally)");
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
{   /* Declarations of close_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccclose_psd_monitor2_xmin;
MCNUM xmax = mccclose_psd_monitor2_xmax;
MCNUM ymin = mccclose_psd_monitor2_ymin;
MCNUM ymax = mccclose_psd_monitor2_ymax;
MCNUM xwidth = mccclose_psd_monitor2_xwidth;
MCNUM yheight = mccclose_psd_monitor2_yheight;
MCNUM radius = mccclose_psd_monitor2_radius;
#line 165 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
   free(PSD_N[0]);
   free(PSD_N);
   free(PSD_p[0]);
   free(PSD_p);
   free(PSD_p2[0]);
   free(PSD_p2);
}
#line 13877 "Monochromator_testing.c"
}   /* End of close_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[16]) fprintf(stderr, "Warning: No xray could reach Component[16] close_psd_monitor2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] close_psd_monitor2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No xray could reach Component[17] Far_E_monitor2\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] Far_E_monitor2=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No xray could reach Component[18] Fardivergence_monitor\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] Fardivergence_monitor=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No xray could reach Component[19] Far_wl_monitor2\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] Far_wl_monitor2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
  /* User FINALLY code for component 'Far_psd_monitor2'. */
  SIG_MESSAGE("Far_psd_monitor2 (Finally)");
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
{   /* Declarations of Far_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccFar_psd_monitor2_xmin;
MCNUM xmax = mccFar_psd_monitor2_xmax;
MCNUM ymin = mccFar_psd_monitor2_ymin;
MCNUM ymax = mccFar_psd_monitor2_ymax;
MCNUM xwidth = mccFar_psd_monitor2_xwidth;
MCNUM yheight = mccFar_psd_monitor2_yheight;
MCNUM radius = mccFar_psd_monitor2_radius;
#line 165 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
   free(PSD_N[0]);
   free(PSD_N);
   free(PSD_p[0]);
   free(PSD_p);
   free(PSD_p2[0]);
   free(PSD_p2);
}
#line 13929 "Monochromator_testing.c"
}   /* End of Far_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[20]) fprintf(stderr, "Warning: No xray could reach Component[20] Far_psd_monitor2\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] Far_psd_monitor2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'source_flat'. */
  SIG_MESSAGE("source_flat (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source_flat");
#define mccompcurname  source_flat
#define mccompcurtype  Source_flat
#define mccompcurindex 2
#define spectrum_file mccsource_flat_spectrum_file
#define prms mccsource_flat_prms
#define square mccsource_flat_square
{   /* Declarations of source_flat=Source_flat() SETTING parameters. */
MCNUM radius = mccsource_flat_radius;
MCNUM yheight = mccsource_flat_yheight;
MCNUM xwidth = mccsource_flat_xwidth;
MCNUM xmin = mccsource_flat_xmin;
MCNUM xmax = mccsource_flat_xmax;
MCNUM dist = mccsource_flat_dist;
MCNUM focus_xw = mccsource_flat_focus_xw;
MCNUM focus_yh = mccsource_flat_focus_yh;
MCNUM E0 = mccsource_flat_E0;
MCNUM dE = mccsource_flat_dE;
MCNUM lambda0 = mccsource_flat_lambda0;
MCNUM dlambda = mccsource_flat_dlambda;
MCNUM flux = mccsource_flat_flux;
MCNUM gauss = mccsource_flat_gauss;
MCNUM randomphase = mccsource_flat_randomphase;
MCNUM phase = mccsource_flat_phase;
#line 213 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Source_flat.comp"
{
  if (square == 1) {
    magnify("xy");
    rectangle("xy",0,0,0,xwidth,yheight);
  } else {
    magnify("xy");
    circle("xy",0,0,0,radius);
  }
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
}
#line 14000 "Monochromator_testing.c"
}   /* End of source_flat=Source_flat() SETTING parameter declarations. */
#undef square
#undef prms
#undef spectrum_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'dmu'. */
  SIG_MESSAGE("dmu (McDisplay)");
  printf("MCDISPLAY: component %s\n", "dmu");
#define mccompcurname  dmu
#define mccompcurtype  Undulator
#define mccompcurindex 3
#define verbose mccdmu_verbose
#define prms mccdmu_prms
#define alpha mccdmu_alpha
#define MELE mccdmu_MELE
{   /* Declarations of dmu=Undulator() SETTING parameters. */
MCNUM E0 = mccdmu_E0;
MCNUM dE = mccdmu_dE;
MCNUM phase = mccdmu_phase;
MCNUM randomphase = mccdmu_randomphase;
MCNUM Ee = mccdmu_Ee;
MCNUM dEe = mccdmu_dEe;
MCNUM Ie = mccdmu_Ie;
MCNUM tbunch = mccdmu_tbunch;
MCNUM t0 = mccdmu_t0;
MCNUM B = mccdmu_B;
MCNUM K = mccdmu_K;
MCNUM gap = mccdmu_gap;
int Nper = mccdmu_Nper;
MCNUM lu = mccdmu_lu;
MCNUM sigey = mccdmu_sigey;
MCNUM sigex = mccdmu_sigex;
MCNUM sigepx = mccdmu_sigepx;
MCNUM sigepy = mccdmu_sigepy;
MCNUM focus_xw = mccdmu_focus_xw;
MCNUM focus_yh = mccdmu_focus_yh;
MCNUM dist = mccdmu_dist;
MCNUM gauss_t = mccdmu_gauss_t;
MCNUM quick_integ = mccdmu_quick_integ;
MCNUM E1st = mccdmu_E1st;
#line 314 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../sources/Undulator.comp"
{
  magnify("xz");
  double zz,dz;
  const double xwidth=1e-2;
  const double D=dist;
  double x0,z0,x1,z1;

  zz=-(prms.length+lu)/2.0;
  dz=lu/2.0;

  while (zz<=(prms.length-lu)/2.0){
    box(0.0,gap/2.0+5e-4,zz,xwidth,1e-3,lu/2.0);
    box(0.0,-gap/2.0-5e-4,zz,xwidth,1e-3,lu/2.0);
    zz+=dz;
  }

  line(0.0,0.0,0.0, K*D*sin(prms.igamma), 0.0, D);
  line(0.0,0.0,0.0,-K*D*sin(prms.igamma), 0.0, D);
  line(0.0,0.0,0.0, 0.0, D*sin(prms.igamma), D);
  line(0.0,0.0,0.0, 0.0,-D*sin(prms.igamma), D);
  
  double phi,dphi;  
  phi =-prms.igamma;
  dphi= 2.0*prms.igamma/32;
  while(phi<prms.igamma){
    x0=D*sin(phi);
    x1=D*sin(phi+dphi);
    z0=D*cos(phi);
    z1=D*cos(phi+dphi);
    line(K*x0,0.0,z0,K*x1,0.0,z1);
    line(0.0,x0,z0,0.0,x1,z1);
    phi+=dphi;
  }
}
#line 14079 "Monochromator_testing.c"
}   /* End of dmu=Undulator() SETTING parameter declarations. */
#undef MELE
#undef alpha
#undef prms
#undef verbose
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mbefore_E_monitor'. */
  SIG_MESSAGE("Mbefore_E_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mbefore_E_monitor");
#define mccompcurname  Mbefore_E_monitor
#define mccompcurtype  E_monitor
#define mccompcurindex 4
#define nE mccMbefore_E_monitor_nE
#define filename mccMbefore_E_monitor_filename
#define E_N mccMbefore_E_monitor_E_N
#define E_p mccMbefore_E_monitor_E_p
#define E_p2 mccMbefore_E_monitor_E_p2
{   /* Declarations of Mbefore_E_monitor=E_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_E_monitor_xmin;
MCNUM xmax = mccMbefore_E_monitor_xmax;
MCNUM ymin = mccMbefore_E_monitor_ymin;
MCNUM ymax = mccMbefore_E_monitor_ymax;
MCNUM xwidth = mccMbefore_E_monitor_xwidth;
MCNUM yheight = mccMbefore_E_monitor_yheight;
MCNUM Emin = mccMbefore_E_monitor_Emin;
MCNUM Emax = mccMbefore_E_monitor_Emax;
MCNUM restore_xray = mccMbefore_E_monitor_restore_xray;
#line 119 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14119 "Monochromator_testing.c"
}   /* End of Mbefore_E_monitor=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mbefore_wl_monitor'. */
  SIG_MESSAGE("Mbefore_wl_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mbefore_wl_monitor");
#define mccompcurname  Mbefore_wl_monitor
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccMbefore_wl_monitor_nL
#define filename mccMbefore_wl_monitor_filename
#define L_N mccMbefore_wl_monitor_L_N
#define L_p mccMbefore_wl_monitor_L_p
#define L_p2 mccMbefore_wl_monitor_L_p2
{   /* Declarations of Mbefore_wl_monitor=L_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_wl_monitor_xmin;
MCNUM xmax = mccMbefore_wl_monitor_xmax;
MCNUM ymin = mccMbefore_wl_monitor_ymin;
MCNUM ymax = mccMbefore_wl_monitor_ymax;
MCNUM xwidth = mccMbefore_wl_monitor_xwidth;
MCNUM yheight = mccMbefore_wl_monitor_yheight;
MCNUM Lmin = mccMbefore_wl_monitor_Lmin;
MCNUM Lmax = mccMbefore_wl_monitor_Lmax;
MCNUM restore_xray = mccMbefore_wl_monitor_restore_xray;
#line 117 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14160 "Monochromator_testing.c"
}   /* End of Mbefore_wl_monitor=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mbefore_psd_monitor'. */
  SIG_MESSAGE("Mbefore_psd_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mbefore_psd_monitor");
#define mccompcurname  Mbefore_psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define nx mccMbefore_psd_monitor_nx
#define ny mccMbefore_psd_monitor_ny
#define nr mccMbefore_psd_monitor_nr
#define filename mccMbefore_psd_monitor_filename
#define restore_xray mccMbefore_psd_monitor_restore_xray
#define PSD_N mccMbefore_psd_monitor_PSD_N
#define PSD_p mccMbefore_psd_monitor_PSD_p
#define PSD_p2 mccMbefore_psd_monitor_PSD_p2
{   /* Declarations of Mbefore_psd_monitor=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccMbefore_psd_monitor_xmin;
MCNUM xmax = mccMbefore_psd_monitor_xmax;
MCNUM ymin = mccMbefore_psd_monitor_ymin;
MCNUM ymax = mccMbefore_psd_monitor_ymax;
MCNUM xwidth = mccMbefore_psd_monitor_xwidth;
MCNUM yheight = mccMbefore_psd_monitor_yheight;
MCNUM radius = mccMbefore_psd_monitor_radius;
#line 174 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14202 "Monochromator_testing.c"
}   /* End of Mbefore_psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PG1_arm'. */
  SIG_MESSAGE("PG1_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PG1_arm");
#define mccompcurname  PG1_arm
#define mccompcurtype  Arm
#define mccompcurindex 7
#line 45 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14230 "Monochromator_testing.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'M2_rotation_arm1'. */
  SIG_MESSAGE("M2_rotation_arm1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "M2_rotation_arm1");
#define mccompcurname  M2_rotation_arm1
#define mccompcurtype  Arm
#define mccompcurindex 8
#line 45 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14249 "Monochromator_testing.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'M2_rotation_arm2'. */
  SIG_MESSAGE("M2_rotation_arm2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "M2_rotation_arm2");
#define mccompcurname  M2_rotation_arm2
#define mccompcurtype  Arm
#define mccompcurindex 9
#line 45 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14268 "Monochromator_testing.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'M2_rotation_arm3'. */
  SIG_MESSAGE("M2_rotation_arm3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "M2_rotation_arm3");
#define mccompcurname  M2_rotation_arm3
#define mccompcurtype  Arm
#define mccompcurindex 10
#line 45 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14287 "Monochromator_testing.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'mirror2'. */
  SIG_MESSAGE("mirror2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "mirror2");
#define mccompcurname  mirror2
#define mccompcurtype  Mirror
#define mccompcurindex 11
#define reflec mccmirror2_reflec
#define reflec_table mccmirror2_reflec_table
#define prms mccmirror2_prms
{   /* Declarations of mirror2=Mirror() SETTING parameters. */
MCNUM zdepth = mccmirror2_zdepth;
MCNUM xwidth = mccmirror2_xwidth;
MCNUM R0 = mccmirror2_R0;
#line 130 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Mirror.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(-xwidth/2.0,0,-zdepth/2.0, xwidth/2.0,0,-zdepth/2.0);
  line(-xwidth/2.0,0, zdepth/2.0, xwidth/2.0,0, zdepth/2.0);
  line(-xwidth/2.0,0,-zdepth/2.0,-xwidth/2.0,0, zdepth/2.0);
  line( xwidth/2.0,0,-zdepth/2.0, xwidth/2.0,0, zdepth/2.0);
}
#line 14314 "Monochromator_testing.c"
}   /* End of mirror2=Mirror() SETTING parameter declarations. */
#undef prms
#undef reflec_table
#undef reflec
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Plane_grating'. */
  SIG_MESSAGE("Plane_grating (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Plane_grating");
#define mccompcurname  Plane_grating
#define mccompcurtype  MultiPurpose_grating
#define mccompcurindex 12
#define blazed mccPlane_grating_blazed
#define display mccPlane_grating_display
#define zdepth mccPlane_grating_zdepth
#define xwidth mccPlane_grating_xwidth
{   /* Declarations of Plane_grating=MultiPurpose_grating() SETTING parameters. */
MCNUM MCangleVariation = mccPlane_grating_MCangleVariation;
MCNUM R0 = mccPlane_grating_R0;
MCNUM r_rho = mccPlane_grating_r_rho;
MCNUM b = mccPlane_grating_b;
MCNUM N_slits = mccPlane_grating_N_slits;
MCNUM d = mccPlane_grating_d;
MCNUM blazed_angle = mccPlane_grating_blazed_angle;
#line 240 "MultiPurpose_grating.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(-xwidth/2.0,0,-zdepth/2.0, xwidth/2.0,0,-zdepth/2.0);
  line(-xwidth/2.0,0, zdepth/2.0, xwidth/2.0,0, zdepth/2.0);
  line(-xwidth/2.0,0,-zdepth/2.0,-xwidth/2.0,0, zdepth/2.0);
  line( xwidth/2.0,0,-zdepth/2.0, xwidth/2.0,0, zdepth/2.0);
}
#line 14350 "Monochromator_testing.c"
}   /* End of Plane_grating=MultiPurpose_grating() SETTING parameter declarations. */
#undef xwidth
#undef zdepth
#undef display
#undef blazed
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Monochromator_Monitor_arm'. */
  SIG_MESSAGE("Monochromator_Monitor_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Monochromator_Monitor_arm");
#define mccompcurname  Monochromator_Monitor_arm
#define mccompcurtype  Arm
#define mccompcurindex 13
#line 45 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 14374 "Monochromator_testing.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Close_wl_monitor1'. */
  SIG_MESSAGE("Close_wl_monitor1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Close_wl_monitor1");
#define mccompcurname  Close_wl_monitor1
#define mccompcurtype  L_monitor
#define mccompcurindex 14
#define nL mccClose_wl_monitor1_nL
#define filename mccClose_wl_monitor1_filename
#define L_N mccClose_wl_monitor1_L_N
#define L_p mccClose_wl_monitor1_L_p
#define L_p2 mccClose_wl_monitor1_L_p2
{   /* Declarations of Close_wl_monitor1=L_monitor() SETTING parameters. */
MCNUM xmin = mccClose_wl_monitor1_xmin;
MCNUM xmax = mccClose_wl_monitor1_xmax;
MCNUM ymin = mccClose_wl_monitor1_ymin;
MCNUM ymax = mccClose_wl_monitor1_ymax;
MCNUM xwidth = mccClose_wl_monitor1_xwidth;
MCNUM yheight = mccClose_wl_monitor1_yheight;
MCNUM Lmin = mccClose_wl_monitor1_Lmin;
MCNUM Lmax = mccClose_wl_monitor1_Lmax;
MCNUM restore_xray = mccClose_wl_monitor1_restore_xray;
#line 117 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14409 "Monochromator_testing.c"
}   /* End of Close_wl_monitor1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'close_E_monitor1'. */
  SIG_MESSAGE("close_E_monitor1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "close_E_monitor1");
#define mccompcurname  close_E_monitor1
#define mccompcurtype  E_monitor
#define mccompcurindex 15
#define nE mccclose_E_monitor1_nE
#define filename mccclose_E_monitor1_filename
#define E_N mccclose_E_monitor1_E_N
#define E_p mccclose_E_monitor1_E_p
#define E_p2 mccclose_E_monitor1_E_p2
{   /* Declarations of close_E_monitor1=E_monitor() SETTING parameters. */
MCNUM xmin = mccclose_E_monitor1_xmin;
MCNUM xmax = mccclose_E_monitor1_xmax;
MCNUM ymin = mccclose_E_monitor1_ymin;
MCNUM ymax = mccclose_E_monitor1_ymax;
MCNUM xwidth = mccclose_E_monitor1_xwidth;
MCNUM yheight = mccclose_E_monitor1_yheight;
MCNUM Emin = mccclose_E_monitor1_Emin;
MCNUM Emax = mccclose_E_monitor1_Emax;
MCNUM restore_xray = mccclose_E_monitor1_restore_xray;
#line 119 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14450 "Monochromator_testing.c"
}   /* End of close_E_monitor1=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'close_psd_monitor2'. */
  SIG_MESSAGE("close_psd_monitor2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "close_psd_monitor2");
#define mccompcurname  close_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define nx mccclose_psd_monitor2_nx
#define ny mccclose_psd_monitor2_ny
#define nr mccclose_psd_monitor2_nr
#define filename mccclose_psd_monitor2_filename
#define restore_xray mccclose_psd_monitor2_restore_xray
#define PSD_N mccclose_psd_monitor2_PSD_N
#define PSD_p mccclose_psd_monitor2_PSD_p
#define PSD_p2 mccclose_psd_monitor2_PSD_p2
{   /* Declarations of close_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccclose_psd_monitor2_xmin;
MCNUM xmax = mccclose_psd_monitor2_xmax;
MCNUM ymin = mccclose_psd_monitor2_ymin;
MCNUM ymax = mccclose_psd_monitor2_ymax;
MCNUM xwidth = mccclose_psd_monitor2_xwidth;
MCNUM yheight = mccclose_psd_monitor2_yheight;
MCNUM radius = mccclose_psd_monitor2_radius;
#line 174 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14492 "Monochromator_testing.c"
}   /* End of close_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Far_E_monitor2'. */
  SIG_MESSAGE("Far_E_monitor2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Far_E_monitor2");
#define mccompcurname  Far_E_monitor2
#define mccompcurtype  E_monitor
#define mccompcurindex 17
#define nE mccFar_E_monitor2_nE
#define filename mccFar_E_monitor2_filename
#define E_N mccFar_E_monitor2_E_N
#define E_p mccFar_E_monitor2_E_p
#define E_p2 mccFar_E_monitor2_E_p2
{   /* Declarations of Far_E_monitor2=E_monitor() SETTING parameters. */
MCNUM xmin = mccFar_E_monitor2_xmin;
MCNUM xmax = mccFar_E_monitor2_xmax;
MCNUM ymin = mccFar_E_monitor2_ymin;
MCNUM ymax = mccFar_E_monitor2_ymax;
MCNUM xwidth = mccFar_E_monitor2_xwidth;
MCNUM yheight = mccFar_E_monitor2_yheight;
MCNUM Emin = mccFar_E_monitor2_Emin;
MCNUM Emax = mccFar_E_monitor2_Emax;
MCNUM restore_xray = mccFar_E_monitor2_restore_xray;
#line 119 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/E_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14536 "Monochromator_testing.c"
}   /* End of Far_E_monitor2=E_monitor() SETTING parameter declarations. */
#undef E_p2
#undef E_p
#undef E_N
#undef filename
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Fardivergence_monitor'. */
  SIG_MESSAGE("Fardivergence_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Fardivergence_monitor");
#define mccompcurname  Fardivergence_monitor
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 18
#define nh mccFardivergence_monitor_nh
#define nv mccFardivergence_monitor_nv
#define rad mccFardivergence_monitor_rad
#define Div_N mccFardivergence_monitor_Div_N
#define Div_p mccFardivergence_monitor_Div_p
#define Div_p2 mccFardivergence_monitor_Div_p2
{   /* Declarations of Fardivergence_monitor=Divergence_monitor() SETTING parameters. */
char* filename = mccFardivergence_monitor_filename;
MCNUM xmin = mccFardivergence_monitor_xmin;
MCNUM xmax = mccFardivergence_monitor_xmax;
MCNUM ymin = mccFardivergence_monitor_ymin;
MCNUM ymax = mccFardivergence_monitor_ymax;
MCNUM xwidth = mccFardivergence_monitor_xwidth;
MCNUM yheight = mccFardivergence_monitor_yheight;
MCNUM maxdiv_h = mccFardivergence_monitor_maxdiv_h;
MCNUM maxdiv_v = mccFardivergence_monitor_maxdiv_v;
MCNUM restore_xray = mccFardivergence_monitor_restore_xray;
MCNUM nx = mccFardivergence_monitor_nx;
MCNUM ny = mccFardivergence_monitor_ny;
MCNUM nz = mccFardivergence_monitor_nz;
#line 148 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
    magnify("xy");
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 14582 "Monochromator_testing.c"
}   /* End of Fardivergence_monitor=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef rad
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Far_wl_monitor2'. */
  SIG_MESSAGE("Far_wl_monitor2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Far_wl_monitor2");
#define mccompcurname  Far_wl_monitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mccFar_wl_monitor2_nL
#define filename mccFar_wl_monitor2_filename
#define L_N mccFar_wl_monitor2_L_N
#define L_p mccFar_wl_monitor2_L_p
#define L_p2 mccFar_wl_monitor2_L_p2
{   /* Declarations of Far_wl_monitor2=L_monitor() SETTING parameters. */
MCNUM xmin = mccFar_wl_monitor2_xmin;
MCNUM xmax = mccFar_wl_monitor2_xmax;
MCNUM ymin = mccFar_wl_monitor2_ymin;
MCNUM ymax = mccFar_wl_monitor2_ymax;
MCNUM xwidth = mccFar_wl_monitor2_xwidth;
MCNUM yheight = mccFar_wl_monitor2_yheight;
MCNUM Lmin = mccFar_wl_monitor2_Lmin;
MCNUM Lmax = mccFar_wl_monitor2_Lmax;
MCNUM restore_xray = mccFar_wl_monitor2_restore_xray;
#line 117 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/L_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14624 "Monochromator_testing.c"
}   /* End of Far_wl_monitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef filename
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Far_psd_monitor2'. */
  SIG_MESSAGE("Far_psd_monitor2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Far_psd_monitor2");
#define mccompcurname  Far_psd_monitor2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define nx mccFar_psd_monitor2_nx
#define ny mccFar_psd_monitor2_ny
#define nr mccFar_psd_monitor2_nr
#define filename mccFar_psd_monitor2_filename
#define restore_xray mccFar_psd_monitor2_restore_xray
#define PSD_N mccFar_psd_monitor2_PSD_N
#define PSD_p mccFar_psd_monitor2_PSD_p
#define PSD_p2 mccFar_psd_monitor2_PSD_p2
{   /* Declarations of Far_psd_monitor2=PSD_monitor() SETTING parameters. */
MCNUM xmin = mccFar_psd_monitor2_xmin;
MCNUM xmax = mccFar_psd_monitor2_xmax;
MCNUM ymin = mccFar_psd_monitor2_ymin;
MCNUM ymax = mccFar_psd_monitor2_ymax;
MCNUM xwidth = mccFar_psd_monitor2_xwidth;
MCNUM yheight = mccFar_psd_monitor2_yheight;
MCNUM radius = mccFar_psd_monitor2_radius;
#line 174 "/zhome/7c/9/7041/mcxtrace/1.4/tools/Python/mxgui/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 14666 "Monochromator_testing.c"
}   /* End of Far_psd_monitor2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_xray
#undef filename
#undef nr
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
/* end of generated C code Monochromator_testing.c */
