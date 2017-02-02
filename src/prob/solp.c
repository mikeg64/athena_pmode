#include "copyright.h"

/*Generating configurations for solar atmosphere using VALIIIc model with SAC
 * becareful the read routine has the ordering wrong normally bb1,bb2,bb3 are the last field read
 * here the last field read is the background density and energy!

 * check and correct the python conversion routine used to generate the input file at
 * https://github.com/mikeg64/smaug/blob/master/smaug/python/slice3dto2d.py

/*see turb.c and par_collision.c for examples of setting up a problem*/
/*============================================================================*/
/*! \file rt.c
 *  \brief Problem generator for RT instabilty.
 *
 * PURPOSE: Problem generator for RT instabilty.  Gravitational pot. is
 *   hardwired to be 0.1z. Density difference is hardwired to be 2.0 in 2D, and
 *   is set by the input parameter <problem>/rhoh in 3D (default value is 3.0).
 *   This reproduces 2D results of Liska & Wendroff, 3D results of
 *   Dimonte et al.
 * 
 * FOR 2D HYDRO:
 * Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to
 * match Liska & Wendroff. Interface is at y=0; perturbation added to Vy
 * Gravity acts in the y-direction.  Special reflecting boundary conditions
 *   added in x2 to improve hydrostatic eqm (prevents launching of weak waves)
 * Atwood number A = (d2-d1)/(d2+d1) = 1/3
 *
 * FOR 3D:
 * Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1
 * Use gamma=5/3 to match Dimonte et al.
 * Interface is at z=0; perturbation added to Vz
 * Gravity acts in the z-direction.  Special reflecting boundary conditions
 *   added in x3 to improve hydrostatic eqm (prevents launching of weak waves)
 * Atwood number A = (d2-d1)/(d2+d1) = 1/2
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - ran2() - random number generator from NR
 * - reflect_ix2() - sets BCs on L-x2 (left edge) of grid used in 2D
 * - reflect_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 * - reflect_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * - reflect_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * - grav_pot2() - gravitational potential for 2D problem (accn in Y)
 * - grav_pot3() - gravitational potential for 3D problem (accn in Z)
 *
 * REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)    */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2() - random number generator from NR
 * reflect_ix2() - sets BCs on L-x2 (left edge) of grid used in 2D
 * reflect_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 * reflect_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * reflect_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * grav_pot2() - gravitational potential for 2D problem (accn in Y)
 * grav_pot3() - gravitational potential for 3D problem (accn in Z)
 *============================================================================*/

static double ran2(long int *idum);
/*static void reflect_ix2(GridS *pGrid);
static void reflect_ox2(GridS *pGrid);
static void reflect_ix3(GridS *pGrid);
static void reflect_ox3(GridS *pGrid);*/
static Real grav_pot2(const Real x1, const Real x2, const Real x3);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);
static int readasciivacconfig(DomainS *pD, char *sacfilename);
static void freadl(FILE *stream, char **string);

char name[50];








/* ========================================================================== */





/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  long int iseed = -1;
  Real amp,x1,x2,x3,lx,ly,lz,rhoh,L_rot,fact;

  double val3c[132][4],rho,pres,val;

  char st1[100],st2[100],st3[100],st4[100],sacfilename[100];
  int ntt;
#ifdef MHD
  Real b0,angle;
#endif
  int ixs, jxs, kxs;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Read perturbation amplitude, problem switch, background density */
  amp = par_getd("problem","amp");
  iprob = par_geti("problem","iprob");
  rhoh  = par_getd_def("problem","rhoh",3.0);
/* Distance over which field is rotated */
  L_rot  = par_getd_def("problem","L_rot",0.0);

/* Read magnetic field strength, angle [should be in degrees, 0 is along +ve
 * X-axis (no rotation)] */
#ifdef MHD
  b0 = par_getd("problem","b0");
  angle = par_getd("problem","angle");
  angle = (angle/180.)*PI;
#endif


  strcpy(sacfilename,par_gets("problem","sacfilename"));


/*w[encode3_uin(p,i,j,ii[2],energyb)]=((pres[i][j]-((bsq[i][j])/2))/(gamma-1))+(bsq[i][j]/2);*/
/* 2D PROBLEM --------------------------------------------------------------- */
/* Initialize two fluids with interface at y=0.0.  Pressure scaled to give a
 * sound speed of 1 at the interface in the light (lower, d=1) fluid 
 * Perturb V2 using single (iprob=1) or multiple (iprob=2) mode 
 */

if (pGrid->Nx[2] == 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        rho=0.0001;
        pres=100;
	pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].E = (pres)/Gamma_1;
	pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        
	pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
#ifdef MHD
	pGrid->B1i[k][j][i] = b0;
	pGrid->U[k][j][i].B1c = b0;
        pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif
      }
#ifdef MHD
    pGrid->B1i[k][j][ie+1] = b0;
#endif
    }
  





/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/

  StaticGravPot = grav_pot2;
  //if (pDomain->Disp[1] == 0) bvals_mhd_fun(pDomain, left_x2,  reflect_ix2);
  //if (pDomain->MaxX[1] == pDomain->RootMaxX[1])
  //  bvals_mhd_fun(pDomain, right_x2, reflect_ox2);

  } /* end of 2D initialization  */



/* Enroll gravitational potential to give accn in z-direction for 3D
 * Use special boundary condition routines.  In 3D, gravity is in the
 * z-direction, so special boundary conditions needed for x3
 */

  StaticGravPot = grav_pot3;

 // if (pDomain->Disp[2] == 0) bvals_mhd_fun(pDomain, left_x3,  reflect_ix3);
 // if (pDomain->MaxX[2] == pDomain->RootMaxX[2])
 //   bvals_mhd_fun(pDomain, right_x3, reflect_ox3);

  } /* end of 3D initialization */


if(readasciivacconfig(pDomain, sacfilename)!=0)
     ath_error("[main]: Bad Restart sac filename: %s\n",sacfilename);





  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll special boundary value functions,
 *    and initialize gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  if (pM->Nx[2] == 1) {
    StaticGravPot = grav_pot2;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      ;//  bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  reflect_ix2);
       ;// bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, reflect_ox2);
      }
    }
  }
 
  if (pM->Nx[2] > 1) {
    StaticGravPot = grav_pot3;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
       ;// bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  reflect_ix3);
       ;// bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, reflect_ox3);
      }
    }
  }

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{



  DomainS *pDomain = (DomainS*)&(pM->Domain[0][0]);
  GridS *pGrid = pM->Domain[0][0].Grid;


int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  Real newtime;

  Real qt,tdep,s_period,AA;
  Real delta_x, delta_z, xxmax, yymax, xxmin, yymin;
  Real exp_x,exp_y,exp_z,exp_xyz;
  Real r1,r2,xp, yp,zp;
  Real vvz;
  Real x1,x2,x3;

  Real xcz,xcx;

  int n1,n2;

  n1=0;
  n2=0;


  s_period=30.0; //Driver period
  AA=350.0;       //Driver amplitude
  //AA=1;
  xcz=0.5e6;
  xcx=2.0e6;
  //delta_z=0.004e6;
  delta_z=0.016e6;
  delta_x=0.016e6;




  if (isnan(pGrid->dt)) ath_error("Time step is NaN!");


	qt=pGrid->time;

	tdep=sin(qt*2.0*3.1415927/s_period);



	if (pM->Nx[2] == 1)
	{
		cc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x2;
		yymax=x1;
		cc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x2;
		yymax=yymax-x1;
		xxmin=x2;
		yymin=x1;
	}

        /*printf("positional info\n");
        printf("%d %d %d \n",is,js,ks);
        printf("%d %d %d \n",ie,je,ke);
        printf("%d %d %d \n", pGrid->Nx[0],pGrid->Nx[1], pGrid->Nx[2]);
        printf("%g %g %g %g\n",xxmin,xxmax,yymin,yymax);*/
	if (pGrid->Nx[2] == 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		zp=x1-xxmin;
		xp=x2-yymin;
		yp=x3;

		r2=(zp-xcz)*(zp-xcz);
                r1=(xp-xcx)*(xp-xcx);
		
                exp_y=exp(-r1/(delta_x*delta_x));
		exp_z=exp(-r2/(delta_z*delta_z));
		//exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
		exp_xyz=exp_y*exp_z;

		vvz=100*AA*exp_xyz*tdep;
                //vvz=0;
                //if(j==12)
                //    printf("%d %d %d %f %f %f %f %f %f\n",i,j,k,xp,yp,zp,xcz,exp_y,exp_z);

//if(i>60 && i<68)
//if(i>is && i<ie)
//{

               // if(i==23  && qt>0)
               //     printf("%d %d %d %g %g %g %g %g \n",i,j,k,AA,exp_z,tdep,exp_y,(pGrid->U[k][j][i].d));


		pGrid->U[k][j][i].M2 += (pGrid->dt)*vvz*(pGrid->U[k][j][i].d);
		//pGrid->U[k][j][i].M2 += (pGrid->dt)*vvz*1.0e-6;
		pGrid->U[k][j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[k][j][i].d)/2.0;
//}
	      }
              //printf("\n");

	    }
	  }
        }

	//for 3D model
	if (pM->Nx[2] > 1)
	{
		cc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x1;
		yymax=x2;
		cc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x1;
		yymax=yymax-x2;
		xxmin=x1;
		yymin=x2;
	}



	if (pGrid->Nx[2] > 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		yp=x2-yymin;
		zp=x3;

		r2=(x3-xcz)*(x3-xcz);
		
		exp_z=exp(-r2/(delta_z*delta_z));
		exp_xyz=sin(PI*xp*(n1+1)/xxmax)*sin(PI*yp*(n2+1)/yymax)*exp_z;

		vvz=AA*exp_xyz*tdep;
                vvz=0;

		pGrid->U[k][j][i].M3 += (pGrid->dt)*vvz*(pGrid->U[k][j][i].d);
		pGrid->U[k][j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[k][j][i].d)/2.0;
	      }

	    }
	  }
      }

	//newtime = pGrid->time + pGrid->dt;



  return;


}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief  Extracted from the Numerical Recipes in C (version 2) code.  
 *   Modified to use doubles instead of floats. - T. A. Gardiner - Aug. 12, 2003
 *   
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  //return 287*x2;
  return 0;
}
/*! \fn static Real grav_pot3(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */
static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
  //return 287*x3;
  return 0;
}




/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x2 for 2D sims
 */

static void reflect_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
        pGrid->U[k][js-j][i].E +=  
	  pGrid->U[k][js+(j-1)][i].d*0.1*(2*j-1)*pGrid->dx2/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js+(j-1)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x2 for 2D sims
 */

static void reflect_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
        pGrid->U[k][je+j][i].E -=
          pGrid->U[k][je-(j-1)][i].d*0.1*(2*j-1)*pGrid->dx2/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* j=je+1 is not set for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 2D sims
 */

static void reflect_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ks-k][j][i].E +=
          pGrid->U[ks+(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 3D sims
 */

static void reflect_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ke+k][j][i].E -=
          pGrid->U[ke-(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif

  return;
}


static void freadl(FILE *stream, char **string)
{    
    unsigned long counter = 0;
    char *line = NULL;
    int next = fgetc(stream);

    do {
        next = fgetc(stream);
        if (next == EOF) {
            free(line);
            break;
        }
        ++counter;
        line = (char*)realloc(line, counter + 1);
        if (line == NULL) {
            puts("line == NULL");
            exit(EXIT_FAILURE);
        }
        line[counter - 1] = (char)next;
    } while (next != '\n');
    line[counter - 1] = '\0';
    *string = (char *)malloc(strlen(line) + 1);
    if (*string == NULL) {
        puts("*string == NULL");
    } else {
        strcpy(*string, line);
    }
    free(line);

}



static int readasciivacconfig(DomainS *pD, char *sacfilename)
{
  int status=0;

  int i,j,k;
  int i1,j1,k1;
  int tn1,tn2,tn3;
  int ni,nj,nk;
  int iif,jf,kf;
  int is,js,ks;
  int ie,je,ke;
  int shift;
  Real x,y,z,val;
  Real t;
  Real x1,x2,x3;
  int it;
  int Nx3T, Nx2T, Nx1T;
  Real *w, *wd;

  char **hlines;

   int ii1,ii2,ii3;
   int pos1, pos2;

   pos1=0;
   pos2=1;

   //ni=p.n[0];
   //nj=p.n[1];
   is=0;
   js=0;
   ks=0;

   iif=ni;
   jf=nj;

  // #ifdef USE_SAC_3D
  // nk=p.n[2];
  // kf=nk;
  // #endif

   GridS *pGrid = pD->Grid;
   is = pGrid->is;  ie = pGrid->ie;
   js = pGrid->js;  je = pGrid->je;
   ks = pGrid->ks;  ke = pGrid->ke;

is=is-nghost;
js=js-nghost;

printf("is,ie %d %d %d %d %d %d \n",is,js,ks,ie,je,ke);
printf("Nx %d %d %d\n",pGrid->Nx[0],pGrid->Nx[1],pGrid->Nx[2]);


 // cc_pos(pGrid,is,i+js,ks,&x1,&x2,&x3);


    
  
printf("reading %s %d %d\n",sacfilename,ni,nj);
   FILE *fdt=fopen(sacfilename,"r+");
//FILE *fdt=fopen("zero1_np0201_001.ini","r+");
   //char **hlines;
   char *line;
   hlines=(char **)calloc(5, sizeof(char*));
   for(i=0; i<5; i++)
       hlines[i]=(char *)calloc(200,sizeof(char));
freadl(fdt, &hlines[0]);
     printf("%s\n", hlines[0]);

fscanf(fdt,"%d %lG %d %d %d\n",&(it),&(t),&ii1,&ii2,&ii3);

  /* Calculate physical size of grid */


  if (pGrid->Nx[1] > 1)
    Nx2T = pGrid->Nx[1] + 2*nghost;
  else
    Nx2T = 1;

  if (pGrid->Nx[0] > 1)
    Nx1T = pGrid->Nx[0] + 2*nghost;
  else
    Nx1T = 1;

ni=Nx1T;
nj=Nx2T;
//nk=;

//   #define NVAR 10
   //#define NDERV 19
//   #define NDERV 19

 w=(Real *)calloc(10*Nx1T*Nx2T,sizeof(Real));
 wd=(Real *)calloc(2*Nx1T*Nx2T,sizeof(Real));


  // typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
  // typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;

fscanf(fdt,"%d %d %d\n",&tn1,&tn2,&tn3);

//printf("here %d %d %d %d %d\n",tn1,tn2,tn3, Nx1T, Nx2T);


   //read 5 header lines
   for(i=3;i<5;i++)
   {
     freadl(fdt, &hlines[i]);
     printf("%s\n", hlines[i]);
   }
//printf("read ascii header %d %d %d %d\n" , p.ipe, is,iif, js,jf);
printf("read ascii header %d %d %d %d \n" , is,iif, js,jf);
  //fscanf(fdt,"%f",&val);
 //printf("%f",val);

//for( j1=js;j1<(je);j1++)
//for( i1=is;i1<(ie);i1++)
for( j1=js;j1<(je);j1++)
for( i1=is;i1<(ie);i1++)
             {

			shift=((j1)*ni+(i1));                         
			//shift=((j1-js)*ni+(i1-is));
                         fscanf(fdt,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&wd[shift],&wd[shift+(ni*nj)],&w[shift],&w[shift+(ni*nj)],&w[shift+(ni*nj*2)],&w[shift+(ni*nj*3)],&w[shift+(ni*nj*4)],&w[shift+(ni*nj*5)],&w[shift+(ni*nj*6)],&w[shift+(ni*nj*7)],&w[shift+(ni*nj*8)],&w[shift+(ni*nj*9)]);

for(k=ks; k<=ke; k++)
{
#ifdef MHD
//pGrid->B1i[k][j][i] =w[shift+(ni*nj*4)]+w[shift+(ni*nj*8)];
//pGrid->B2i[k][j][i] =w[shift+(ni*nj*5)]+w[shift+(ni*nj*9)];
i=i1;
j=j1;
pGrid->B1i[k][j][i] =w[shift+(ni*nj*4)]+w[shift+(ni*nj*6)];
pGrid->B2i[k][j][i] =w[shift+(ni*nj*5)]+w[shift+(ni*nj*7)];
pGrid->B3i[k][j][i] =0.0;
#endif

        //add background contributions from sac input file
	//pGrid->U[k][j][i].d = w[shift]+w[shift+(ni*nj*7)];
        //pGrid->U[k][j][i].E = w[shift+(ni*nj*3)]+w[shift+(ni*nj*6)];
	pGrid->U[k][j][i].d = w[shift]+w[shift+(ni*nj*9)];
        pGrid->U[k][j][i].E = w[shift+(ni*nj*3)]+w[shift+(ni*nj*8)];

	pGrid->U[k][j][i].M1 = w[shift+(ni*nj)];
        pGrid->U[k][j][i].M2 = w[shift+(ni*nj*2)];
        pGrid->U[k][j][i].M3 =0 ;


              //  if(j==63)
              //      printf("%d %d %d %g %g %g %g  \n",i,j,k,(pGrid->dt), w[shift]+w[shift+(ni*nj*9)],w[shift+(ni*nj*3)]+w[shift+(ni*nj*8)] ,(pGrid->U[k][j][i].d));



}

              }


	      fclose(fdt);


        free(w);
        free(wd);

//printf("here end:%d %d %d %d\n",is,js,ie,je);

	return status;
}

