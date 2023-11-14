#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            11

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  NGAM                           0
#define  MBH                            1
#define  RHOM                           2
#define  RDISK                          3
#define  DIST                           4
#define  ALPHA                          5
#define  EDOT                           6
#define  RHOMINF                        7
#define  PRSMINF                        8
#define  RESTART                        9
#define  RESET_VEL                      10
/* [Beg] user-defined constants (do not change this line) */

#define WARNING_MESSAGES NO
#define SHOCK_FLATTENING MULTID
#define LIMITER          VANLEER_LIM

#define UNIT_DENSITY     1.e9
#define UNIT_LENGTH      1.e8
#define UNIT_VELOCITY    3.e10

/* [End] user-defined constants (do not change this line) */
