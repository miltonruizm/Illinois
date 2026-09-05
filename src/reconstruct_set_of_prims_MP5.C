/*****************************************
 * MP5 Reconstruction Interface.
 *
 * Fifth-order monotonicity-preserving
 * reconstruction following
 * Suresh & Huynh (1997).
 *****************************************/

#define MINUS2 0
#define MINUS1 1
#define PLUS0  2
#define PLUS1  3
#define PLUS2  4
#define MAXNUMINDICES 5
//      ^^^^^^^^^^^^^ Be _sure_ to define MAXNUMINDICES appropriately!

// You'll find the #define's for LOOP_DEFINE and SET_INDEX_ARRAYS inside:
#include "loop_defines_reconstruction.h"

static inline CCTK_REAL minmod(const CCTK_REAL a,
                               const CCTK_REAL b)
{
    if(a*b <= 0.0) return 0.0;
    return (fabs(a) < fabs(b)) ? a : b;
}

static inline CCTK_REAL minmod4(
        const CCTK_REAL a,
        const CCTK_REAL b,
        const CCTK_REAL c,
        const CCTK_REAL d)
{
    if(a*b<=0.0) return 0.0;
    if(a*c<=0.0) return 0.0;
    if(a*d<=0.0) return 0.0;

    const CCTK_REAL s = (a>0.0 ? 1.0 : -1.0);

    return s*std::min(std::min(fabs(a),fabs(b)),
                      std::min(fabs(c),fabs(d)));
}

static inline CCTK_REAL median(const CCTK_REAL a,
                               const CCTK_REAL b,
                               const CCTK_REAL c)
{
    return a + minmod(b-a,c-a);
}

//===============================================================
// Unlimited fifth-order MP5 interface reconstruction
//===============================================================

static inline CCTK_REAL MP5_right(
    const CCTK_REAL Um2,
    const CCTK_REAL Um1,
    const CCTK_REAL U0,
    const CCTK_REAL Up1,
    const CCTK_REAL Up2) {
    const CCTK_REAL U5 = (2.0*Um2 - 13.0*Um1 + 47.0*U0 + 27.0*Up1 - 3.0*Up2)/60.0;

    // Monotone predictor
    const CCTK_REAL alpha = 4.0;
    const CCTK_REAL UMP   = U0 + minmod(Up1 - U0,alpha*(U0 - Um1));

    constexpr CCTK_REAL epsm = 1.0e-10;
    if ((U5 - U0)*(U5 - UMP) <= epsm)
      return U5;
    
    const CCTK_REAL Dm2  = Um2 - 2.0*Um1 + U0;
    const CCTK_REAL Dm1  = Um1 - 2.0*U0  + Up1;
    const CCTK_REAL D0   = U0  - 2.0*Up1 + Up2;

    const CCTK_REAL Dm1M = minmod4(4.0*Dm1 - D0, 4.0*D0  - Dm1,Dm1,D0);
    const CCTK_REAL Dm2M = minmod4(4.0*Dm2 - Dm1,4.0*Dm1 - Dm2,Dm2,Dm1);


    const CCTK_REAL UUL  = U0 + alpha*(U0 - Um1);
    const CCTK_REAL UAV  = 0.5*(U0 + Up1);
    const CCTK_REAL UMD  = UAV - 0.5*Dm1M;
    const CCTK_REAL ULC  = U0 + 0.5*(U0 - Um1) + (4.0/3.0)*Dm2M;

    const CCTK_REAL UMIN = std::max(std::min(U0,std::min(Up1,UMD)),std::min(U0,std::min(UUL,ULC)));
    const CCTK_REAL UMAX = std::min(std::max(U0,std::max(Up1,UMD)),std::max(U0,std::max(UUL,ULC)));

    return std::max(UMIN,std::min(U5,UMAX));
}

static inline CCTK_REAL MP5_left(
    const CCTK_REAL Um2,
    const CCTK_REAL Um1,
    const CCTK_REAL U0,
    const CCTK_REAL Up1,
    const CCTK_REAL Up2) {
    const CCTK_REAL U5 = (-3.0*Um2 + 27.0*Um1 + 47.0*U0  - 13.0*Up1  + 2.0*Up2)/60.0;

    // Monotone predictor
    const CCTK_REAL alpha = 4.0;
    const CCTK_REAL UMP   = U0 + minmod(Um1 - U0,alpha*(U0 - Up1));

    constexpr CCTK_REAL epsm = 1.0e-10;
    if ((U5 - U0)*(U5 - UMP) <= epsm)
      return U5;

    const CCTK_REAL Dm2  = Um2 - 2.0*Um1 + U0;
    const CCTK_REAL Dm1  = Um1 - 2.0*U0  + Up1;
    const CCTK_REAL D0   = U0  - 2.0*Up1 + Up2;

    const CCTK_REAL Dm1M = minmod4(4.0*Dm1 - D0, 4.0*D0  - Dm1,Dm1,D0);
    const CCTK_REAL Dm2M = minmod4(4.0*Dm2 - Dm1,4.0*Dm1 - Dm2,Dm2,Dm1);
    
    const CCTK_REAL UUL =  U0 + alpha*(U0 - Up1);
    const CCTK_REAL UAV =  0.5*(U0 + Um1);
    const CCTK_REAL UMD = UAV - 0.5*Dm1M;
    const CCTK_REAL ULC = U0 + 0.5*(U0 - Up1) + (4.0/3.0)*Dm2M;

    const CCTK_REAL UMIN  = std::max(std::min(Um1,std::min(U0,UMD)),std::min(U0,std::min(UUL,ULC))); 
    const CCTK_REAL UMAX  = std::min(std::max(Um1,std::max(U0,UMD)),std::max(U0,std::max(UUL,ULC)));

    return std::max(UMIN,std::min(U5,UMAX));
}


static void reconstruct_set_of_prims_MP5( const igm_eos_parameters eos,
                                          const cGH *restrict cctkGH,
                                          const int *restrict cctk_lsh,
                                          const int flux_dirn,
                                          const int num_prims_to_reconstruct,
                                          const int *restrict which_prims_to_reconstruct,
                                          gf_and_gz_struct *restrict in_prims,
                                          gf_and_gz_struct *restrict out_prims_r,
                                          gf_and_gz_struct *restrict out_prims_l,
                                          CCTK_REAL *restrict temporary ) {

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES],dU[MAXNUMVARS][MAXNUMINDICES],slope_lim_dU[MAXNUMVARS][MAXNUMINDICES],
    Ur[MAXNUMVARS][MAXNUMINDICES],Ul[MAXNUMVARS][MAXNUMINDICES];
  int ijkgz_lo_hi[4][2];

  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];

    if(in_prims[whichvar].gz_lo[flux_dirn]!=0 || in_prims[whichvar].gz_hi[flux_dirn]!=0) {
      CCTK_VError(VERR_DEF_PARAMS,"TOO MANY GZ'S! WHICHVAR=%d: %d %d %d : %d %d %d DIRECTION %d",whichvar,
		  in_prims[whichvar].gz_lo[1],in_prims[whichvar].gz_lo[2],in_prims[whichvar].gz_lo[3],
		  in_prims[whichvar].gz_hi[1],in_prims[whichvar].gz_hi[2],in_prims[whichvar].gz_hi[3],flux_dirn);
    }


    // *** LOOP 1: Interpolate to Ur and Ul, which are face values ***
    //  You will find that Ur depends on U at MINUS1,PLUS0, PLUS1,PLUS2, and
    //                     Ul depends on U at MINUS2,MINUS1,PLUS0,PLUS1.
    //  However, we define the below loop from MINUS2 to PLUS2. Why not split
    //     this up and get additional points? The reason is that later on,
    //     Ur and Ul depend on ftilde, which is defined from MINUS2 to PLUS2,
    //     so we would lose those points anyway.
    LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-2,2,flux_dirn);
      /* *** LOOP 1a: READ INPUT *** */
      // Read in a primitive at all gridpoints between m = MINUS2 & PLUS2, where m's direction is given by flux_dirn. Store to U.
      for(int ii=MINUS2;ii<=PLUS2;ii++) U[whichvar][ii] = in_prims[whichvar].gf[index_arr[flux_dirn][ii]];


      /* *** LOOP 1b: DO COMPUTATION *** */
      /* First, compute simple dU = U(i) - U(i-1), where direction of i
       *         is given by flux_dirn, and U is a primitive variable:
       *         {rho_b,P,vx,vy,vz,Bx,By,Bz}. */
      // Note that for Ur and Ul at i, we must compute dU(i-1),dU(i),dU(i+1),
      //         and dU(i+2)

      // Compute MP5 reconstruction (Temporary placeholder)


      const CCTK_REAL Um2 = U[whichvar][MINUS2];
      const CCTK_REAL Um1 = U[whichvar][MINUS1];
      const CCTK_REAL U0  = U[whichvar][PLUS0];
      const CCTK_REAL Up1 = U[whichvar][PLUS1];
      const CCTK_REAL Up2 = U[whichvar][PLUS2];

      Ur[whichvar][PLUS0] = MP5_right(Um2,Um1,U0,Up1,Up2);
      Ul[whichvar][PLUS0] = MP5_left(Um2,Um1,U0,Up1,Up2);
      
      // Prevent new extrema (temporary monotonicity limiter)
      /*
      CCTK_REAL Umin = std::min(U0,Up1);
      CCTK_REAL Umax = std::max(U0,Up1);

      Ur[whichvar][PLUS0] = std::min(std::max(Ur[whichvar][PLUS0],Umin),Umax);
      Umin = std::min(Um1,U0);
      Umax = std::max(Um1,U0);

      Ul[whichvar][PLUS0] = std::min(std::max(Ul[whichvar][PLUS0],Umin),Umax);
      */

      //Ur = fifth-order polynomial;
      //Ul = fifth-order polynomial;
 
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ur[whichvar][PLUS0];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = Ul[whichvar][PLUS0];
    }


    // *** LOOP 2: STEEPEN RHOB ***
    // Note that this loop applies ONLY to RHOB.
    if(whichvar==RHOB) {
      LOOP_DEFINE(2,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
	SET_INDEX_ARRAYS(-2,2,flux_dirn);
	// Set rho and P separately, since within this loop,

	// Read in all primitives between MINUS2 & PLUS2. Store to U.
	for(int ii=MINUS2;ii<=PLUS2;ii++) U[RHOB][ii]     = in_prims[RHOB    ].gf[index_arr[flux_dirn][ii]];
	for(int ii=MINUS1;ii<=PLUS1;ii++) U[PRESSURE][ii] = in_prims[PRESSURE].gf[index_arr[flux_dirn][ii]];
	Ur[RHOB][PLUS0] = out_prims_r[RHOB].gf[index_arr[flux_dirn][PLUS0]];
	Ul[RHOB][PLUS0] = out_prims_l[RHOB].gf[index_arr[flux_dirn][PLUS0]];

	dU[whichvar][MINUS1] = U[whichvar][MINUS1]- U[whichvar][MINUS2];
	dU[whichvar][PLUS0]  = U[whichvar][PLUS0] - U[whichvar][MINUS1];
	dU[whichvar][PLUS1]  = U[whichvar][PLUS1] - U[whichvar][PLUS0];
	dU[whichvar][PLUS2]  = U[whichvar][PLUS2] - U[whichvar][PLUS1];


	// Steepen rho
	// DEPENDENCIES: RHOB face values, RHOB(MINUS2,MINUS1,PLUS0,PLUS1,PLUS2), P(MINUS1,PLUS0,PLUS1), and slope_lim_dU[RHOB](MINUS1,PLUS1)

	// Output rho
	out_prims_r[RHOB].gf[index_arr[flux_dirn][PLUS0]] = Ur[RHOB][PLUS0];
	out_prims_l[RHOB].gf[index_arr[flux_dirn][PLUS0]] = Ul[RHOB][PLUS0];
      }
    }
  }


  // *** LOOP 4: SHIFT Ur AND Ul ***
  /* Currently face values are set so that
   *      a) Ur(i) represents U(i+1/2), and
   *      b) Ul(i) represents U(i-1/2)
   *    Here, we shift so that the indices are consistent:
   *      a) U(i-1/2+epsilon) = oldUl(i)   = newUr(i)
   *      b) U(i-1/2-epsilon) = oldUr(i-1) = newUl(i)
   *    Note that this step is not strictly necessary if you keep
   *      track of indices when computing the flux. */
  for(int ww=0;ww<num_prims_to_reconstruct;ww++) {
    int whichvar=which_prims_to_reconstruct[ww];
    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(-1,0,flux_dirn);
      temporary[index_arr[flux_dirn][PLUS0]] = out_prims_r[whichvar].gf[index_arr[flux_dirn][MINUS1]];
    }

    LOOP_DEFINE(3,2,  cctk_lsh,flux_dirn,  ijkgz_lo_hi,in_prims[whichvar].gz_lo,in_prims[whichvar].gz_hi) {
      SET_INDEX_ARRAYS(0,0,flux_dirn);
      // Then shift so that Ur represents the gridpoint at i-1/2+epsilon,
      //                and Ul represents the gridpoint at i-1/2-epsilon.
      // Ur(i-1/2) = Ul(i-1/2)     = U(i-1/2+epsilon)
      // Ul(i-1/2) = Ur(i+1/2 - 1) = U(i-1/2-epsilon)
      out_prims_r[whichvar].gf[index_arr[flux_dirn][PLUS0]] = out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]];
      out_prims_l[whichvar].gf[index_arr[flux_dirn][PLUS0]] = temporary[index_arr[flux_dirn][PLUS0]];
    }
    // Ul was just shifted, so we lost another ghostzone.
    out_prims_l[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_l[whichvar].gz_hi[flux_dirn]+=0;
    // As for Ur, we didn't need to get rid of another ghostzone,
    //    but we did ... seems wasteful!
    out_prims_r[whichvar].gz_lo[flux_dirn]+=1;
    out_prims_r[whichvar].gz_hi[flux_dirn]+=0;

  }
}


// standard Colella-Woodward parameters:
//    K0 = 0.1d0, eta1 = 20.0, eta2 = 0.05, epsilon = 0.01d0
#define K0      0.1
#define ETA1   20.0
#define ETA2    0.05
#define PPM_EPSILON 0.01


  CCTK_REAL contact_discontinuity_check = Gamma*K0*fabs(U[RHOB][PLUS1]-U[RHOB][MINUS1])*
    MIN(U[PRESSURE][PLUS1],U[PRESSURE][MINUS1])
    -fabs(U[PRESSURE][PLUS1]-U[PRESSURE][MINUS1])*MIN(U[RHOB][PLUS1],U[RHOB][MINUS1]);
  CCTK_REAL second_deriv_check = -d2rho_b_p1*d2rho_b_m1;
  CCTK_REAL relative_change_check = fabs(2.0*d1rho_b) - PPM_EPSILON*MIN(U[RHOB][PLUS1],U[RHOB][MINUS1]);

  if(contact_discontinuity_check >= 0.0 && second_deriv_check >= 0.0
     && relative_change_check >= 0.0) {

    CCTK_REAL eta_tilde=0.0;
    if (fabs(d1rho_b) > 0.0) {
      eta_tilde = -(1.0/6.0)*(d2rho_b_p1-d2rho_b_m1)/(2.0*d1rho_b);
    }
    CCTK_REAL eta = MAX(0.0,MIN(ETA1*(eta_tilde - ETA2),1.0));

    // Next compute Urp1 and Ul for RHOB, using the MC prescription:
    // Ur_p1 = U_p1   - 0.5*slope_lim_dU_p1
    CCTK_REAL rho_br_mc_p1 = U[RHOB][PLUS1] - 0.5*slope_lim_dU[RHOB][PLUS1];
    // Ul = U_m1 + 0.5*slope_lim_dU_m1
    // Based on this line of code, Ur[index] = a_j - \delta_m a_j / 2. (cf. Eq. 65 in Marti & Muller's "PPM Method for 1D Relativistic Hydro." paper)
    //    So: Ur[indexp1] = a_{j+1} - \delta_m a_{j+1} / 2. This is why we have rho_br_mc[indexp1]
    CCTK_REAL rho_bl_mc    = U[RHOB][MINUS1] + 0.5*slope_lim_dU[RHOB][MINUS1];

    rho_bl_ppm[PLUS0] = rho_bl_ppm[PLUS0]*(1.0-eta) + rho_bl_mc*eta;
    rho_br_ppm[PLUS0] = rho_br_ppm[PLUS0]*(1.0-eta) + rho_br_mc_p1*eta;

  }
}

