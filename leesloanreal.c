#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesloanreal.h"

void readsloanreal()
{

  FILE *pf;
  char filename[200];
  int i,j;
  int ngalssloan;
  float disred;
  double dismin = 1.E26, dismax = -1.E26;
  double alfamin = 1.E26, alfamax = -1.E26;
  double deltamin = 1.E26, deltamax = -1.E26;

  cp.omegam = 0.3                       ;  /* OMEGA MATERIA                              */
  cp.omegal = 0.7                       ;  /* OMEGA LAMBDA                               */
  cp.omegak = 1.0-cp.omegam-cp.omegal   ;  /* OMEGA CURVATURA                            */
  cp.h0     = 100.                      ;  /* ESTO DEJA TODO EN UNIDADES DE H^-1         */

  RED("Read Sloan...\n");

  sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");

  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  ngalssloan = cp.npart;

  /* ALOCATACION Y LECTURA */
  gal = (struct galsloanreal *) malloc(ngalssloan*sizeof(struct galsloanreal));
  P   = (struct particle_data *) malloc(ngalssloan*sizeof(struct particle_data));

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  i=0;
  for(j=0;j<ngalssloan;j++)
  {
		fscanf(pf,"%d %d %d %f %f %f %f %f %f %f %f %f %f %d %f %f %f", \
    &gal[i].id,&gal[i].id_parent,&gal[i].id_NYU_VAGC,&gal[i].alfa,&gal[i].delta,&gal[i].red,&gal[i].mr,&gal[i].mr_lim,\
    &gal[i].completenness,&gal[i].KE_petro[0],&gal[i].KE_petro[1],&gal[i].KE_model[0],&gal[i].KE_model[1], \
    &gal[i].id_red_source_type,&gal[i].red_FOG,&gal[i].red_kaiser,&gal[i].red_real);

    //if(gal[i].red<=redmin || gal[i].red >=redmax || gal[i].mr > rmaplim || gal[i].mr < rmapmin) continue;

    gal[i].alfa  *= PI180;
    gal[i].delta *= PI180;

    disred = red2dis(gal[i].red_real); /*EN MPC*/

    P[i].Pos[0] =  disred*cos(gal[i].delta)*cos(gal[i].alfa) ;
    P[i].Pos[1] =  disred*cos(gal[i].delta)*sin(gal[i].alfa) ;
    P[i].Pos[2] =  disred*sin(gal[i].delta)                  ;
    P[i].Dis    = disred;

    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];

    if(P[i].Dis > dismax) dismax = P[i].Dis;
    if(P[i].Dis < dismin) dismin = P[i].Dis;

    if(gal[i].alfa > alfamax) alfamax = gal[i].alfa;
    if(gal[i].alfa < alfamin) alfamin = gal[i].alfa;
    if(gal[i].delta > deltamax) deltamax = gal[i].delta;
    if(gal[i].delta < deltamin) deltamin = gal[i].delta;

    i++;
  }

  cp.npart = i;
  cp.vol = (alfamax-alfamin)*(cos(deltamin)-cos(deltamax))*(pow(dismax,3.)-pow(dismin,3.))/3.0;
  cp.vol = fabs(cp.vol);

  fprintf(stdout,"DisMin %.8e DisMax %.8e\n",dismin,dismax);
  fprintf(stdout,"AlfaMin %.8e AlfaMax %.8e\n",alfamin,alfamax);
  fprintf(stdout,"DeltaMin %.8e DeltaMax %.8e\n",deltamin,deltamax);
  fprintf(stdout,"Volumen aprox %.8e\n",cp.vol);
  fprintf(stdout,"Num Total %d\n",cp.npart);
  fflush(stdout);
  RED("End Read Sloan\n");

  fclose(pf);

}

double funlum(double x, void *p)
{
  struct paramfl *par=(struct paramfl *)p ;
  double ma=(par->ma)                     ;
  double alfa=(par->alfa)                 ;
  double t                                ;

  t=pow(10.0,0.4*(ma-x)) ;
  t=pow(t,1.0+alfa)*exp(-t) ;

  return(t) ;
}

double intfl(double x1, double x2)
{
  double resultado, error;
  size_t neval;
  struct paramfl pfl;
  gsl_function F; 

  pfl.ma=flma;
  pfl.alfa=flalfa;
  pfl.fia=flfia;

  F.function = &funlum;
  F.params = &pfl;

  if(x1<MAGMENOSINF)x1=MAGMENOSINF;
  gsl_integration_qng(&F,x1,x2,1.0e-7,1.0e-7,&resultado,&error,&neval);
  return(resultado);
}

double f(double z, void *p) /* FUNCION A INTEGRAR PARA LA DISTANCIA EN FUNCION DE Z*/
{
  struct paramcos *par=(struct paramcos *)p ;
  double om=(par->omegam);
  double ol=(par->omegal);
  double ok=(par->omegak);
  double q;
  q=pow(1.+z,3.)*om+pow(1.+z,2.)*ok+ol;
  return(1.0/sqrt(q));
}

double red2dis(double z)
{
  double resultado, error;
  size_t neval;
  struct paramcos pcos;
  gsl_function F;

  pcos.omegam=cp.omegam;
  pcos.omegal=cp.omegal;
  pcos.omegak=cp.omegak;

  F.function=&f;
  F.params=&pcos;

  gsl_integration_qng(&F,0.0,z,1.0e-7,1.0e-7,&resultado,&error,&neval);
 
  return(CVEL/cp.h0*resultado);
}

void change_positions(int n)
{
  int ip, idim;
  #ifdef PRINT_XYZ
  FILE *pfout;
  char filename[200];
  #endif
  RED("Inicio Change Positions\n");

  #ifdef PRINT_XYZ
  sprintf(filename,"sloanreal_xyz.bin");
  pfout=fopen(filename,"w");
  fwrite(&n,sizeof(int),1,pfout);
  #endif

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip = 0; ip < n; ip++)
  {
    #ifdef PRINT_XYZ
    fwrite(&P[ip].Pos[0],sizeof(float),1,pfout);
    fwrite(&P[ip].Pos[1],sizeof(float),1,pfout);
    fwrite(&P[ip].Pos[2],sizeof(float),1,pfout);
    #endif

    for(idim = 0; idim < 3; idim++)
      P[ip].Pos[idim] -= pmin[idim];
  }

  #ifdef PRINT_XYZ
  fclose(pfout);
  #endif

  cp.lbox = pmax[0] - pmin[0];
  for(idim = 1; idim < 3; idim++)
    if(cp.lbox < (pmax[idim] - pmin[idim])) cp.lbox = (pmax[idim] - pmin[idim]);

  cp.lbox *= 1.001;
  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
}

