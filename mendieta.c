#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesloanreal.h"
#include "timer.h"
#include "iden.h"
#include "colores.h"
#include "grid.h"

void Write_Groups(double *fof);
void Write_Sloan(double *fof);

int main(int argc, char **argv)
{
  int    i;
  double start,end;
  #ifdef LEN_MANUEL
  double rmablim;  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  double rmabmin;  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  double rintlim;  // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
  #endif

  TIMER(start);
  
  init_variables(argc,argv);

  omp_set_nested(1);
  
  readsloanreal();

  // Cambia origen de coordenadas
  change_positions(cp.npart);

  for(i=0;i<cp.npart;i++) P[i].sub = 0;

  #ifdef LEN_MANUEL
  rmablim = rmaplim-25.0-5.0*log10(red2dis(zcut)*(1.+zcut));  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = rmapmin-25.0-5.0*log10(red2dis(zcut)*(1.+zcut));  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  rintlim = intfl(rmabmin,rmablim);                           // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
  fprintf(stdout, "INTEGRAL LIMITES\n");
  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmabmin, rmapmin);
  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmablim, rmaplim);
  fflush(stdout);
  #endif

  for(i=0;i<nfrac;i++)
  {
    fprintf(stdout, "\nBegins Identification : Step %d of %d \n",i+1,nfrac);
   
    #ifdef LEN_MANUEL
      iden.r0 = 3.0/(4.0*M_PI*(fof[i]+1)*(0.4*log(10.0)*flfia*rintlim));
      iden.r0 = cbrt(iden.r0);
    #else
      iden.r0 = cbrt(cp.vol/(cp.npart*(fof[i]+1))); // cbrt(vol/(npart*(fof+1)))
    #endif

    iden.step = i;
    iden.nobj = cp.npart;

    fprintf(stdout,"Linking length = %f \n",iden.r0);

    identification();

    Write_Groups(fof);
    
    free(Temp.head);
    free(Temp.npgrup);
    free(Temp.ll);
  }

  /************* TERMINO LA IDENTIFICACION ***************/
  Write_Sloan(fof);

  free(P);
  free(gal);
  grid_free();
  free(fof);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

void Write_Groups(double *fof)
{
  int i,j,k,id,npar,gn,save_sub;
  float dx,dy,dz,xc,yc,zc;
  char filename[200];
  FILE *pfout, *pfcentros, *pfcentros_ascii;

  //////////////////////////////////////////////////////
  if(iden.step!=0)
  {
    sprintf(filename,"%.2f_%.2f_fof.bin",fof[1],fof[0]);
    pfout=fopen(filename,"w");
    i = iden.ngrupos-1;
    fwrite(&i,sizeof(int),1,pfout);
    sprintf(filename,"%.2f_%.2f_centros.bin",fof[1],fof[0]);
    pfcentros=fopen(filename,"w");
    fwrite(&i,sizeof(int),1,pfcentros);
    sprintf(filename,"%.2f_%.2f_centros.dat",fof[1],fof[0]);
    pfcentros_ascii=fopen(filename,"w");
  }else{
    sprintf(filename,"%.2f_%.2f_fof.bin",fof[0],fof[1]);
    pfout=fopen(filename,"w");
    i = iden.ngrupos-1;
    fwrite(&i,sizeof(int),1,pfout);
  }
  //////////////////////////////////////////////////////

  npar = gn = 0;

  for(i=1;i<iden.ngrupos;i++)
  {

    j = 0;
    id = k = Temp.head[i];
    if(iden.step!=0)
    {
      xc = yc = zc = 0.0;
      save_sub = P[k].sub;
      fwrite(&save_sub,sizeof(int),1,pfout);
    }
    fwrite(&i,sizeof(int),1,pfout);
    fwrite(&Temp.npgrup[i],sizeof(int),1,pfout);

    while(k != -1)
    {
      if(iden.step!=0)
      {
        // cuidado con el orden {pos[i]-centro} en este caso
        dx = P[k].Pos[0] - P[id].Pos[0];
        dy = P[k].Pos[1] - P[id].Pos[1];
        dz = P[k].Pos[2] - P[id].Pos[2];

        #ifdef PERIODIC
        dx = dx > cp.lbox*0.5 ? dx-cp.lbox : dx;
        dy = dy > cp.lbox*0.5 ? dy-cp.lbox : dy;
        dz = dz > cp.lbox*0.5 ? dz-cp.lbox : dz;
  
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif

        xc += dx;
        yc += dy;
        zc += dz;
      }else{   
        P[k].sub = P[k].gr;
      }
      //fwrite(&P[k].id,sizeof(int),1,pfout);
      fwrite(&k,sizeof(int),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    #ifdef DEBUG
    assert(j == Temp.npgrup[i]);
    #endif

    if(iden.step!=0)
    {
      xc /= (float)Temp.npgrup[i];
      yc /= (float)Temp.npgrup[i];
      zc /= (float)Temp.npgrup[i];

      xc += P[id].Pos[0];
      yc += P[id].Pos[1];
      zc += P[id].Pos[2];
      
      xc += pmin[0];
      yc += pmin[1];
      zc += pmin[2];
           
      #ifdef PERIODIC
      xc = xc<0 ? cp.lbox+(float)fmod(xc,cp.lbox) : (float)fmod(xc,cp.lbox);
      yc = yc<0 ? cp.lbox+(float)fmod(yc,cp.lbox) : (float)fmod(yc,cp.lbox);
      zc = zc<0 ? cp.lbox+(float)fmod(zc,cp.lbox) : (float)fmod(zc,cp.lbox);
      #endif

      fwrite(&save_sub,sizeof(int),1,pfcentros);
      fwrite(&i,sizeof(int),1,pfcentros);
      fwrite(&xc,sizeof(float),1,pfcentros);
      fwrite(&yc,sizeof(float),1,pfcentros);
      fwrite(&zc,sizeof(float),1,pfcentros);
      fwrite(&Temp.npgrup[i],sizeof(int),1,pfcentros);
      fprintf(pfcentros_ascii,"%d %d %f %f %f %d\n",save_sub,i,xc,yc,zc,Temp.npgrup[i]);
    }  

    npar+=j;
    gn++;
  }

  assert(gn == (iden.ngrupos-1));
  fclose(pfout);

  if(iden.step!=0)
  {
    fclose(pfcentros);
    fclose(pfcentros_ascii);
  }

  fprintf(stdout,"num de grupos %d num de particulas en grupos %d\n",gn,npar);
  fflush(stdout);

  return;
}

void Write_Sloan(double *fof)
{
  int i;
  char filename[200];
  FILE *pfout;
  int **npgrupo;
  struct galsloang galg;

  npgrupo = (int **) calloc(cp.npart,sizeof(int *));
  for(i=0;i<cp.npart;i++)
    npgrupo[i] = (int *) calloc(nfrac,sizeof(int));

  for(i=0;i<cp.npart;i++)
  {
    npgrupo[P[i].sub][0]++;
    npgrupo[P[i].gr][1]++;
  }

  sprintf(filename,"catreal_%.2f_%.2f.dat",fof[0],fof[1]);
  pfout=fopen(filename,"w");
  fwrite(&cp.npart,sizeof(int),1,pfout);

  for(i=0;i<cp.npart;i++)
  {
    galg.gal       = gal[i];
    galg.gal.alfa  /= PI180;
    galg.gal.delta /= PI180;
    galg.grupo[0]   = P[i].sub;
    galg.grupo[1]   = P[i].gr;
    galg.npgrupo[0] = npgrupo[P[i].sub][0];
    galg.npgrupo[1] = npgrupo[P[i].gr][1];
    fwrite(&galg,sizeof(struct galsloang),1,pfout);
  } 

  for(i=0;i<cp.npart;i++)
    free(npgrupo[i]);
  free(npgrupo);

  return;
}
