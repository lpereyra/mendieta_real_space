#ifndef LEESLOAN_H
#define LEESLOAN_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define PI180 0.01745329251994329577
#define CVEL  299792.458
#define MAGMENOSINF -26.0

struct SnapST
{
  int nfiles;
  char root[200], name[200]; 
} snap;

struct paramfl
{
  double ma   ;
  double alfa ;
  double fia  ;
};

struct paramcos
{
  double omegam   ;
  double omegal   ;
  double omegak   ;
};

struct galsloanreal
{
  int id,id_parent,id_NYU_VAGC; //galaxy ID, galaxy ID in parent catalog, NYU-VAGC ID : object ID (row in NYU-VAGC files)
  float alfa,delta,red;         //ra (in degree), dec (in degree), z redshift
  float mr, mr_lim;             //apparent magnitude r_band (SDSS magnitude), magnitude limit r_band
  float completenness;          //completeness in this region
  float KE_petro[2];            //^{0.1}M_r-5\log h  (K+E corrected to z=0.1) -- Petro ; ^{0.1}(g-r)  color (K+E corrected to z=0.1) -- Petro
  float KE_model[2];            //^{0.1}M_r-5\log h  (K+E corrected to z=0.1) -- Model; ^{0.1}(g-r)  color (K+E corrected to z=0.1) -- Model
  int   id_red_source_type;     //redshift source type: 1 SDSS; 2 other; 3 KSG-VAGC; 4 nearest 5 from X-ray clusters (last 30) for those with negative values (469), photometries are update according to their parents. 
  float red_FOG;                //re-FOG space redshift (correct for the Kaiser effect only)
  float red_kaiser;             //re-Kaiser redshift (correct for the FOG effect only)
  float red_real;               //re-real redshift (correct for both the Kaiser and FOG effect)
} *gal;

struct galsloang
{
  struct galsloanreal gal;
  int grupo[2], npgrupo[2];
} *gr;

void readsloanreal();
double funlum(double x, void *p);
double intfl(double x1, double x2);
double f(double z, void *p);
double red2dis(double z);
void change_positions(int n);

#endif
