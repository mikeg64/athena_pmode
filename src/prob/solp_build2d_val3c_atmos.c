#include "copyright.h"
/*============================================================================*/
/*! \file rotor.c
 *  \brief Sets up 2D rotor test problem.
 *
 * PURPOSE: Sets up 2D rotor test problem.  The center of the grid is assumed to
 *   have coordinates (x1,x2) = [0,0]; the grid initialization must be
 *   consistent with this
 *
 * REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes",
 *   JCP, 161, 605 (2000)						      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"



#ifndef MHD
#error : The rotor problem can only be run with MHD.
#endif
#ifdef ISOTHERMAL 
#error : The rotor problem can only be run with an ADIABATIC eos.
#endif


static void ry_bc(GridS *pGrid);
static Real grav_pot2(const Real x1, const Real x2, const Real x3);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);
static int readasciivacconfig(DomainS *pD, char *sacfilename);
static void freadl(FILE *stream, char **string);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,is,ie,js,je,ks,ke;
  Real v0,p0,bx0,x1,x2,x3,rad,frac,r0,r1;

  Real amp,lx,ly,lz;
  Real a,v1,fac,w;
  double val3c[132][4],rho,pres,val;
  
  char st1[100],st2[100],st3[100],st4[100],sacfilename[100];
  int ntt;

#ifdef MHD
  Real b0,angle;
#endif


/* Read initial conditions from 'athinput' */

  v0 = par_getd("problem","v0");
  p0 = par_getd("problem","p0");
  bx0 = par_getd("problem","bx0");
  r0 = par_getd("problem","r0");
  r1 = par_getd("problem","r1");

  
  strcpy(sacfilename,par_gets("problem","sacfilename"));
  

/* Initialize the grid.  Note the center is always assumed to have coordinates
 * x1=0, x2=0; the grid range in the input file must be consistent with this */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];


if(readasciivacconfig(pDomain, sacfilename)!=0)
     ath_error("[main]: Bad Restart sac filename: %s\n",sacfilename);


             for(i=0; i<132; i++)
                 if((i+js)<=je)
		{
		  cc_pos(pGrid,is,i+js,ks,&x1,&x2,&x3);
                 // printf("%f %f %f %f %f\n",x2, val3c[i][0], val3c[i][1], val3c[i][2], val3c[i][3]);
                //if(p->ipe==1)
                 }

/*lagrange interpolation*/
/*   % t1=(xval-x(i+1))/(x(i)-x(i-1));
     % t2=(xval-x(i))/(x(i+1)-x(i));
     % y =t1*f(i)+t2*f(i+1); */  








//bvals_mhd_fun(pDomain, right_x2, ry_bc);


/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/


StaticGravPot = grav_pot2;

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

void problem_read_restart(MeshS *pM, FILE *fp)
{
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
  Real delta_x, delta_y, delta_z, xxmax, yymax, xxmin, yymin;
  Real exp_x,exp_y,exp_z,exp_xyz;
  Real r1,r2,xp, yp,zp;
  Real vvz;
  Real x1,x2,x3;

  Real xcz,xcx;

  int n1,n2;

  n1=2;
  n2=2;


  s_period=180.0; //Driver period
  AA=350.0;       //Driver amplitude
  //AA=1;
  xcz=0.5e6;
  xcx=2.0e6;
  delta_z=0.016e6;
  delta_x=0.016e6;
  delta_y=0.016e6;

printf("source term1 \n");


  if (isnan(pGrid->dt)) ath_error("Time step is NaN!");


	qt=pGrid->time;

	tdep=sin(qt*2.0*PI/s_period);
        //tdep=1.0;


	if (pM->Nx[2] == 1)
	{
		fc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x1;
		yymax=x3;
		fc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x1;
		yymax=yymax-x3;
		xxmin=x1;
		yymin=x3;
	}

        /*printf("%d %d %d \n",is,js,ks);
        printf("%d %d %d \n",ie,je,ke);
        printf("%d %d %d \n", pGrid->Nx[0],pGrid->Nx[1], pGrid->Nx[2]);*/
	if (pGrid->Nx[2] == 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		fc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		yp=x3-yymin;
		zp=x2;

		r2=(zp-xcz)*(zp-xcz);
                r1=(xp-xcx)*(xp-xcx);
		
                exp_y=exp(-r1/(delta_x*delta_x));
		exp_z=exp(-r2/(delta_z*delta_z));
                exp_x=exp(-r1/(delta_y*delta_y));

		exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
		//exp_xyz=exp_y*exp_z;
                //exp_xyz=exp_x*exp_z;

		vvz=AA*exp_xyz*tdep;
                //vvz=0;
                //if(j==12)
                //    printf("%d %d %d %f %f %f %f %f %f %f\n",i,j,k,xp,yp,zp,xcz,exp_x,exp_z,vvz);

//if(i>60 && i<68)
//if(i>is && i<ie)
//{

//                if(j>8 && j<16 && qt<2)
//                    printf("%d %d %d %g %g %g %g  \n",i,j,k,vvz,exp_x,exp_z,(pGrid->dt)*vvz*(pGrid->U[k][j][i].d));


		pGrid->U[k][j][i].M2 += (pGrid->dt)*vvz*(pGrid->U[k][j][i].d);
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
		fc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x1;
		yymax=x2;
		fc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x1;
		yymax=yymax-x2;
		xxmin=x1;
		yymin=x2;
	}



	if (pGrid->Nx[2] > 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		fc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		yp=x2-yymin;
		zp=x3;

		r2=(x3-xcz)*(x3-xcz);
		
		exp_z=exp(-r2/(delta_z*delta_z));
		exp_xyz=sin(PI*xp*(n1+1)/xxmax)*sin(PI*yp*(n2+1)/yymax)*exp_z;

		vvz=AA*exp_xyz*tdep;
                //vvz=0;

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

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return 287*x2;
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



/*----------------------------------------------------------------------------*/
/*! \fn static void ry_bc(GridS *pG)
 *  \brief  Apply boundary condition in right-y direction
 */

static void ry_bc(GridS *pGrid)
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
        pGrid->U[k][je+j][i].M2 = 0; /* reflect 2-mom. */
        //pGrid->U[k][je+j][i].E -=
        //  pGrid->U[k][je-(j-1)][i].d*0.1*(2*j-1)*pGrid->dx2/Gamma_1;
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

for( j1=js;j1<(je);j1++)
for( i1=is;i1<(ie);i1++)
             {

                         shift=((j1-js)*ni+(i1-is));
                         fscanf(fdt,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&wd[shift],&wd[shift+(ni*nj)],&w[shift],&w[shift+(ni*nj)],&w[shift+(ni*nj*2)],&w[shift+(ni*nj*3)],&w[shift+(ni*nj*4)],&w[shift+(ni*nj*5)],&w[shift+(ni*nj*6)],&w[shift+(ni*nj*7)],&w[shift+(ni*nj*8)],&w[shift+(ni*nj*9)]);

for(k=ks; k<=ke; k++)
{
#ifdef MHD
pGrid->B1i[k][j][i] =w[shift+(ni*nj*4)]+w[shift+(ni*nj*8)];
pGrid->B2i[k][j][i] =w[shift+(ni*nj*5)]+w[shift+(ni*nj*9)];
pGrid->B3i[k][j][i] =0.0;
#endif

        //add background contributions from sac input file
	pGrid->U[k][j][i].d = w[shift]+w[shift+(ni*nj*3)];
        pGrid->U[k][j][i].E = w[shift+(ni*nj*3)]+w[shift+(ni*nj*6)];
	pGrid->U[k][j][i].M1 = w[shift+(ni*nj)];
        pGrid->U[k][j][i].M2 = w[shift+(ni*nj*2)];
        pGrid->U[k][j][i].M3 =0 ;
}

              }


	      fclose(fdt);


        free(w);
        free(wd);

//printf("here end:%d %d %d %d\n",is,js,ie,je);

	return status;
}
