/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
const double UNIT_MASS    = UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
const double UNIT_TIME    = UNIT_LENGTH/UNIT_VELOCITY;
const double UNIT_ENERGY  = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
const double UNIT_G       = (6.67e-8)/UNIT_ENERGY/UNIT_LENGTH*UNIT_MASS*UNIT_MASS;
/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  
  g_smallDensity  = 0.1*g_inputParam[RHOM]/UNIT_DENSITY * g_inputParam[RHOMINF];
  g_smallPressure = 1e-10;  
  v[RHO] = g_inputParam[RHOM]/UNIT_DENSITY*g_inputParam[RHOMINF]; // Sets units density
  v[PRS] = UNIT_G*v[RHO]*g_inputParam[MBH]/UNIT_MASS/x1; // set ambient medium pressure profile
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  double entA = UNIT_G*g_inputParam[MBH]/UNIT_MASS/2./(g_inputParam[NGAM]+1.)/pow(g_inputParam[RHOM]/UNIT_DENSITY, 1./g_inputParam[NGAM])*(1.-1./g_inputParam[DIST]); // sets entropy constant
  double rhoc = pow(UNIT_G*g_inputParam[MBH]/UNIT_MASS/(g_inputParam[NGAM]+1.)/entA/g_inputParam[RDISK]*UNIT_LENGTH*(g_inputParam[RDISK]/x1/UNIT_LENGTH-0.5*(g_inputParam[RDISK]*g_inputParam[RDISK]/UNIT_LENGTH/UNIT_LENGTH/(sin(x2)+1.e-16)/(sin(x2)+1.e-16))/x1/x1-0.5/g_inputParam[DIST]),g_inputParam[NGAM]);// sets density, need to clean up
  if (rhoc > v[RHO]){
    v[RHO] = rhoc;
    v[VX3] = sqrt(UNIT_G*g_inputParam[MBH] * g_inputParam[RDISK]/UNIT_LENGTH/UNIT_MASS)/x1/sin(x2); //sets velocity
    v[TRC] = 1.;
  }

  
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 * DCARREL Initialize to be in equilibrium
 *********************************************************************** */
{
  // Initial the pressure so that the domain is initially in
  // hydrostatic equilibrium
  int i,j,k, I,J,K, id;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  if (g_inputParam[RESTART]){
    printf("Resetting domain from /restart...\n");
    //density
    id = InputDataOpen("restart/rho.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[RHO][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    //pressure
    id = InputDataOpen("restart/prs.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[PRS][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    //velocity 1
    id = InputDataOpen("restart/vx1.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[VX1][k][j][i] = g_inputParam[RESET_VEL]*InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    //velocity 2
    id = InputDataOpen("restart/vx2.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[VX2][k][j][i] = g_inputParam[RESET_VEL]*InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    //velocity 3
    id = InputDataOpen("restart/vx3.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[VX3][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    // tracer
    id = InputDataOpen("restart/tr1.dbl", "restart/grid.out", " ",0, CENTER);
    TOT_LOOP(k,j,i){d->Vc[TRC][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);}
    InputDataClose(id);
    
    return;
  }

  double vsq, rhocurr, rcurr, drcurr;
  DOM_LOOP(k,j,i){
    if (i == IBEG){continue;} //outer is already set
    K=k;J=j;
    I = IEND - (i-IBEG);
    vsq = d->Vc[VX3][K][J][I]*d->Vc[VX3][K][J][I];
    rcurr = grid->xgc[IDIR][I];
    rhocurr = d->Vc[RHO][K][J][I];
    drcurr  = grid->xgc[IDIR][I+1] - grid->xgc[IDIR][I];
    // sets pressure
    d->Vc[PRS][K][J][I] = d->Vc[PRS][K][J][I+1] + drcurr*rhocurr/rcurr*(UNIT_G*g_inputParam[MBH]/UNIT_MASS/rcurr-vsq);
  }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
       BOX_LOOP(box,k,j,i){
        double vx1 = d->Vc[VX1][k][j][IBEG];
        if (vx1 > 0.){
          d->Vc[VX1][k][j][i] = 0.0;
        }
        else{
          d->Vc[VX1][k][j][i] = vx1;
        }
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IBEG];
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][IBEG];
	d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
      }
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
     BOX_LOOP(box,k,j,i){
        double vx1 = d->Vc[VX1][k][j][IEND];
        if (vx1 < 0.){
          d->Vc[VX1][k][j][i] = 0.0;
        }
        else{
          d->Vc[VX1][k][j][i] = vx1;
        }
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND];
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][IEND];
	d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND];
      }
  

    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return -UNIT_G*g_inputParam[MBH]/UNIT_MASS/x1;
}
#endif
