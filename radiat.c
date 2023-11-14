/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute right hand side for Tabulated cooling

  \authors A. Mignone (mignone@ph.unito.it)\n
           M. Sormani\n

 \b References

  \date   Apr 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
extern double UNIT_TIME, UNIT_MASS, UNIT_ENERGY;
/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  mu, T, Tmid, scrh, dT, prs;
  static double *L_tab, *T_tab, E_cost;
  double nH = UNIT_DENSITY/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
  double ne = nH*(1.0 + 0.5*CONST_AZ*FRAC_Z);
/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */
/*
  if (T_tab == NULL){

          [Mach = 2.779757]
    FILE *fcool;
    printLog (" > Reading table from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      printLog ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
0
    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    fclose(fcool);
  }
*/
/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */
/*
  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }

  mu  = MeanMolecularWeight(v);
  T   = prs/v[RHO]*KELVIN*mu;

  if (T != T){
    printf ("! Radiat(): Nan found: rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }
*/
/*
  if (T < g_minCoolingTemp) { 
    rhs[RHOE] = 0.0;
    return;
  }
*/
/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */
/*
  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi] || T < T_tab[klo]){
    rhs[RHOE] = 0.0;
    return;
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }
*/
/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */
/*
  dT        = T_tab[khi] - T_tab[klo];
  scrh      = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
  rhs[RHOE] = -ne*nH*scrh*E_cost;
  */
/* ----------------------------------------------
    Temperature cutoff
   ---------------------------------------------- */
  
  if( g_inputParam[EDOT])
  {
    double rhoep, init_time;
    init_time = 1./UNIT_TIME; //want to start at 1 second
    rhoep = 2.e10/UNIT_ENERGY/UNIT_MASS*UNIT_TIME*v[TRC]*v[RHO]* pow( (g_time+init_time)/(86400/UNIT_TIME), -1.3);
    rhs[RHOE] = rhoep;
  }
  else{
    rhs[RHOE] = 0;
  }

}
