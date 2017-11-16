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
        rho=1.0e-12;
        pres=0;
        b0=0;
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

 // StaticGravPot = grav_pot3;

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
  int jea,jsa;
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


  s_period=180.0; //Driver period
  AA=350.0;       //Driver amplitude
  //AA=0.0;
  //AA=1;
  xcz=0.5e6;
  xcx=2.0e6;
  //delta_z=0.004e6;
  //delta_z=0.016e6;
  //delta_x=0.016e6;
delta_z=0.08e6;
  delta_x=0.08e6;



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

        //printf("positional info\n");
        //printf("%d %d %d \n",is,js,ks);
        //printf("%d %d %d \n",ie,je,ke);
        /*printf("%d %d %d \n", pGrid->Nx[0],pGrid->Nx[1], pGrid->Nx[2]);
        printf("%g %g %g %g\n",xxmin,xxmax,yymin,yymax);*/
	if (pGrid->Nx[2] == 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		zp=x2-yymin;
		yp=x3;

		r2=(zp-xcz)*(zp-xcz);
                r1=(xp-xcx)*(xp-xcx);
		
                exp_y=exp(-r1/(delta_x*delta_x));
		exp_z=exp(-r2/(delta_z*delta_z));
		//exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
		exp_xyz=exp_y*exp_z;

                vvz=0.0;

                //if(qt>20)
		   vvz=AA*exp_xyz*tdep;
                //vvz=0;
                //if(j==14)
                //    printf("%d %d %d %f %f %f %f %f %f %f %f\n",i,j,k,xp-xcx,zp-xcz,tdep,xcz,xcx,exp_y,exp_z,delta_x, delta_z);

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


	  for (k=0; k<pGrid->Nx[2]; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=0; i<pGrid->Nx[0]; i++) {

                jea=pGrid->Nx[1];
                jsa=0;
                //pGrid->U[k][jsa+3][i].M2
                pGrid->U[k][jsa+6][i].M2 = pGrid->U[k][jsa+7][i].M2;
                pGrid->U[k][jsa+5][i].M2 = pGrid->U[k][jsa+6][i].M2;
		pGrid->U[k][jsa+4][i].M2 = pGrid->U[k][jsa+5][i].M2;
                pGrid->U[k][jsa+3][i].M2=pGrid->U[k][jsa+4][i].M2;
		pGrid->U[k][jsa+2][i].M2 = pGrid->U[k][jsa+3][i].M2;
		pGrid->U[k][jsa+1][i].M2 = pGrid->U[k][jsa+2][i].M2;
		pGrid->U[k][jsa][i].M2 = pGrid->U[k][jsa+1][i].M2;

                //pGrid->U[k][jea-5][i].M2 = pGrid->U[k][jea-6][i].M2;

                pGrid->U[k][jea-6][i].M2 = pGrid->U[k][jea-7][i].M2;

                pGrid->U[k][jea-5][i].M2 = pGrid->U[k][jea-6][i].M2;
		pGrid->U[k][jea-4][i].M2 = pGrid->U[k][jea-5][i].M2;
                
		pGrid->U[k][jea-3][i].M2 = pGrid->U[k][jea-4][i].M2;
		pGrid->U[k][jea-2][i].M2 = pGrid->U[k][jea-3][i].M2;
		pGrid->U[k][jea-1][i].M2 = pGrid->U[k][jea-2][i].M2;
		pGrid->U[k][jea][i].M2 = pGrid->U[k][jea-1][i].M2;
}
}
}



	  //for (k=ks; k<=ke; k++) {
	    for (j=0; j<pM->Nx[1]; j++) {
	      //for (i=is; i<=ie; i++) {
                 
                 k=ke;
                 printf("%d %d %g %g %g\n", j,k, pGrid->U[k][j][6].M2/pGrid->U[k][j][6].d, pGrid->U[k][j][60].M2/pGrid->U[k][j][60].d,pGrid->U[k][j][120].M2/pGrid->U[k][j][120].d);
//printf("%d %d %g %g %g\n", i,k, pGrid->U[k][j][6].d, pGrid->U[k][j][60].d,pGrid->U[k][j][110].d);

		//	}
		//}
	}

	//newtime = pGrid->time + pGrid->dt;



  return;


}

void Userwork_after_loop(MeshS *pM)
{
}


/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot1(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot1(const Real x1, const Real x2, const Real x3)
{
  //return 287*x2;
  return 0;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  //return 28.7*x2;
  // 1  return 259*x2;  -0.1
  // 2 return 256*x2;  0.054
  // 3 return 257*x2;   0.01
  // 4 return 257.5*x2; -0.03
  // 5 return 257.2*x2;
  //6  return 257.15*x2; -0.18 after10
  //7 return 257.05*x2;  -76t -107b 0.02
  //8 return 257.1*x2;  76t -107b
   return 257.045*x2;
  //return 0;
}
/*! \fn static Real grav_pot3(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */
static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
  //return 287*x3;
  return 0;
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

   //iif=ni;
   //jf=nj;

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
//w=(Real *)calloc(10*ie*je,sizeof(Real));
// wd=(Real *)calloc(2*ie*je,sizeof(Real));

  // typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
  // typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;

fscanf(fdt,"%d %d %d\n",&tn1,&tn2,&tn3);
k=0;

//printf("here %d %d %d %d %d\n",tn1,tn2,tn3, Nx1T, Nx2T);
printf("here in config %d %d %d %d %d %d %d\n", k,is,js,ks,ie,je,ke);


   //read 5 header lines
   for(i=3;i<5;i++)
   {
     freadl(fdt, &hlines[i]);
     printf("%s\n", hlines[i]);
   }
//printf("read ascii header %d %d %d %d\n" , p.ipe, is,iif, js,jf);
printf("read ascii header %d %d %d %d %d\n" , nghost, ni,nj,is, js);
  //fscanf(fdt,"%f",&val);
 //printf("%f",val);

//for( j1=js;j1<(je);j1++)
//for( i1=is;i1<(ie);i1++)

//for(k=ks; k<=ke; k++)
for( i1=is;i1<(ni);i1++)
for( j1=js;j1<(nj);j1++)
//for( i1=is;i1<(ie);i1++)
             {

			//shift=((j1)*ni+(i1));                         
			//shift=((j1-js)*ni+(i1-is));
                        // fscanf(fdt,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&wd[shift],&wd[shift+(ni*nj)],&w[shift],&w[shift+(ni*nj)],&w[shift+(ni*nj*2)],&w[shift+(ni*nj*3)],&w[shift+(ni*nj*4)],&w[shift+(ni*nj*5)],&w[shift+(ni*nj*6)],&w[shift+(ni*nj*7)],&w[shift+(ni*nj*8)],&w[shift+(ni*nj*9)]);
//for( i1=is;i1<(ie);i1++)
//for(k=ks; k<=ke; k++)
{

//shift=((j1)*ni+(i1));
shift=((i1)*nj+(j1));

 fscanf(fdt,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&wd[shift],&wd[shift+(ni*nj)],&w[shift],&w[shift+(ni*nj)],&w[shift+(ni*nj*2)],&w[shift+(ni*nj*3)],&w[shift+(ni*nj*4)],&w[shift+(ni*nj*5)],&w[shift+(ni*nj*6)],&w[shift+(ni*nj*7)],&w[shift+(ni*nj*8)],&w[shift+(ni*nj*9)]);


//if(i1==0)
//   printf("here in config1 %d %d  %d %g %g\n", i1,j1,shift,w[shift+(ni*nj*8)],w[shift+(ni*nj*9)]);


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

//top boundary
/*for( i1=is;i1<(ie);i1++)
{
	pGrid->B1i[k][j][i] =;
	pGrid->B2i[k][j][i] =;
	pGrid->B3i[k][j][i] =;

	pGrid->U[k][j][i].d = ;
        pGrid->U[k][j][i].E = ;

	pGrid->U[k][j][i].M1 = ;
        pGrid->U[k][j][i].M2 = ;
        pGrid->U[k][j][i].M3 = ;

}*/


//bottom boundary
/*for( i1=is;i1<(ie);i1++)
{
	pGrid->B1i[k][j][i] =;
	pGrid->B2i[k][j][i] =;
	pGrid->B3i[k][j][i] =;

	pGrid->U[k][j][i].d = ;
        pGrid->U[k][j][i].E = ;

	pGrid->U[k][j][i].M1 = ;
        pGrid->U[k][j][i].M2 = ;
        pGrid->U[k][j][i].M3 = ;
}*/


//left boundary
/*for( j=js;j<(je);j++)
{
	pGrid->B1i[k][j][i] =;
	pGrid->B2i[k][j][i] =;
	pGrid->B3i[k][j][i] =;

	pGrid->U[k][j][i].d = ;
        pGrid->U[k][j][i].E = ;

	pGrid->U[k][j][i].M1 = ;
        pGrid->U[k][j][i].M2 = ;
        pGrid->U[k][j][i].M3 = ;
}*/



//right boundary



	      fclose(fdt);


        free(w);
        free(wd);

//printf("here end:%d %d %d %d\n",is,js,ie,je);

	return status;
}




