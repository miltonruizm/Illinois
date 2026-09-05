#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "con2prim_helpers.h"
#include "inlined_functions.h"

/*****************************************************************************/
/*********** IMPORTED FROM A STRIPPED VERSION OF THE HEADER FILES ************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>


static bool use_entropy_final(
    const igm_eos_parameters eos,
    const double rho,
    const double temp_energy,
    const double depsdT,
    const double rel_T,
    const double rel_P,
    const double rel_eps);

static void dual_energy_compare(
    const igm_eos_parameters eos,
    const double x,
    const double *restrict con,
    const double *restrict param,
    const double temp_guess,
    double *restrict prim, 
    output_stats& stats);


void palenzuela( const igm_eos_parameters eos,
                 const double S_squared,
                 const double BdotS,
                 const double B_squared,
                 const double *restrict con,
                 double *restrict prim,
                 const double *restrict SU,
                 const double tol_x,
                 output_stats& stats );

double zbrent( double (*func)(const igm_eos_parameters, double, double *restrict, double *restrict, output_stats&),
               const igm_eos_parameters eos,
               double *restrict param,
               double *restrict temp_guess,
               double x1,
               double x2,
               double tol_x,
               output_stats& stats );

static void recover_from_energy(
    const igm_eos_parameters eos,
    double rho,
    double Ye,
    double eps,
    double *press,
    double *ent,
    double *temp);

static void recover_from_entropy(
    const igm_eos_parameters eos,
    double rho,
    double Ye,
    double ent,
    double *press,
    double *eps,
    double *temp); 

void calc_prim( const igm_eos_parameters eos,
                const double x,
                const double *restrict con,
                const double * param,
                const double temp_guess,
                double *restrict prim,
                const double S_squared,
                const double BdotS,
                const double B_squared,
                const double *restrict SU,
                output_stats& stats );

double func_root( const igm_eos_parameters eos,
                  double x,
                  double *restrict param,
                  double *restrict temp_guess,
                  output_stats& stats );

inline int check_depsdT_condition( const igm_eos_parameters eos,
                                   const CCTK_REAL *restrict param,
                                   const CCTK_REAL T_guess,
                                   const CCTK_REAL x ) {
  // Now check whether or not to use the entropy equation
  double q       = param[par_q];
  double r       = param[par_r];
  double s       = param[par_s];
  double t       = param[par_t];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2        = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  double W       = pow(Wminus2, -0.5);

  double rho     = param[conDD]/W;
  double ye      = param[conYE]/param[conDD];
  double temp    = T_guess;
  double eps     = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );
  double press   = 0.0;
  double ent     = 0.0;
  double depsdT  = 0.0;

  // Now perform an EOS call. The goal here is to determine
  // deps/dT, which should allow us to select the appropriate
  // variable that we will use to recover the temperature.
  enforce_table_bounds_rho_Ye_eps( eos,&rho,&ye,&eps );
    WVU_EOS_P_S_T_and_depsdT_from_rho_Ye_eps( rho,ye,eps, &press,&ent,&depsdT,&temp);

    //-- debug

    //    static int checkprint=0; if(checkprint<20) {printf("CHECK: rho=%e Ye=%e eps=%e depsdT=%e threshold=%e\n",rho,ye,eps,depsdT,eos.depsdT_threshold);checkprint++;}
    //if (depsdT >= eos.depsdT_threshold) {printf("ENERGY SELECTED: " "rho=%e Ye=%e eps=%e " "depsdT=%e threshold=%e\n",rho, ye, eps,depsdT,eos.depsdT_threshold);}

    int con2prim_key;
    if( eos.evolve_entropy && (depsdT < eos.depsdT_threshold) ) {
      con2prim_key = Palenzuela1D_entropy;
    }
    else {
      con2prim_key = Palenzuela1D;
    }


    /*/-
    
    printf("C2P_CHECK: "
           "x=%e "
           "rho=%e "
           "Ye=%e "
           "eps=%e "
           "T=%e "
           "depsdT=%e "
           "threshold=%e "
           "key=%d\n",
           x,
           rho,
           ye,
           eps,
           temp,
           depsdT,
           eos.depsdT_threshold,
           con2prim_key);
    //--
    */


    static int c2p_check_count = 0;
    if (c2p_check_count < 500) {
      printf("C2P_CHECK: x=%e rho=%e Ye=%e eps=%e T=%e depsdT=%e threshold=%e branch=%s\n",
	     x, rho, ye, eps, temp, depsdT, eos.depsdT_threshold,
	     (con2prim_key == Palenzuela1D_entropy ? "entropy" : "energy"));
      c2p_check_count++;
    }
    
    return con2prim_key;

}

/*****************************************************************************/
/************** ILLINOISGRMHD INTERFACE FOR PALENZUELA CON2PRIM **************/
/*****************************************************************************/
//
// This routine expects 8 conservs as input:
//
// -> D     = rho_star      / sqrt(gamma) = W * rho (W is the Lorentz factor)
// -> DYe   = Ye_Star       / sqrt(gamma)
// -> B^{i} = \tilde{B}^{i} / sqrt(gamma) / sqrt(4pi)
// -> S_{i} = \tilde{S}_{i} / sqrt(gamma)
//
// From the input quantities, we compute B_{i} and S^{i}
int con2prim_Palenzuela1D( const igm_eos_parameters eos,
                           const CCTK_REAL *restrict adm_quantities,
                           const CCTK_REAL *restrict con,
                           CCTK_REAL *restrict prim,
                           output_stats& stats ) {

  //--

  // TEMP DIAGNOSTIC  -- raw specific entropy from conservatives, pre-con2prim

  if (con[DD] > 1e-7 && con[DD] < 1e-4) {
      const double ent_specific_raw = con[WS] / con[DD];  // ~ S/D, W cancels approx at low v
        fprintf(stderr,
		"RAW_WS: D=%e WS=%e ent_specific=%e\n",
		con[DD], con[WS], ent_specific_raw);}
  
  //--
  
  // Set gamma_{ij} and gamma^{ij}
  CCTK_REAL gammaDD[3][3],gammaUU[3][3];
  set_gammaDD_and_gammaUU_from_ADM_quantities(adm_quantities,gammaDD,gammaUU);

  // Read in B^{i}
  CCTK_REAL BU[3];
  for(int i=0;i<3;i++) BU[i] = con[B1_con+i];

  // Read in S_{i}
  CCTK_REAL SD[3];
  for(int i=0;i<3;i++) SD[i] = con[S1_cov+i];

  // Compute B_{i} = gamma_{ij}B^{j}
  CCTK_REAL BD[3];
  raise_or_lower_indices_3d( BU,gammaDD, BD );

  // Compute S^{i} = gamma^{ij}S_{j}
  CCTK_REAL SU[3];
  raise_or_lower_indices_3d( SD,gammaUU, SU );

  // S^2 = S^i * S_i
  CCTK_REAL S_squared = 0.0;
  for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];

  // Enforce ceiling on S^{2} (A5 of Palenzuela et al. https://arxiv.org/pdf/1505.01607.pdf)
  CCTK_REAL S_squared_max = SQR( con[DD] + con[TAU] );
  if( S_squared > 0.9999 * S_squared_max ) {
    // Compute rescaling factor
    CCTK_REAL rescale_factor_must_be_less_than_one = sqrt(0.9999*S_squared_max/S_squared);
    // Rescale S_{i}
    for(int i=0;i<3;i++) SD[i] *= rescale_factor_must_be_less_than_one;
    // S_{i} has been rescaled. Recompute S^{i}.
    raise_or_lower_indices_3d( SD,gammaUU, SU );
    // Now recompute S^{2} := gamma^{ij}S^{i}S_{j}.
    S_squared = 0.0;
    for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];
    // Check if the fix was successful
    if( simple_rel_err(S_squared,0.9999*S_squared_max) > 1e-12 ) CCTK_VError(VERR_DEF_PARAMS,"Incompatible values of S_squared after rescaling: %.15e %.15e\n",S_squared,0.9999*S_squared_max);
  }

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];

  // B^2 = B^i * B_i
  CCTK_REAL B_squared = 0.0;
  for(int i=0;i<3;i++) B_squared += BU[i] * BD[i];

  // Now proceed as usual. In this first function call,
  // we allow the entropy to be used (if needed).
  const CCTK_REAL tolerance = 1e-10;
  stats.which_routine       = None;

  stats.c2p_failed          = true;


  palenzuela( eos, S_squared,BdotS,B_squared, con, prim, SU, tolerance, stats );

  
  if( (stats.c2p_failed == true) && (stats.which_routine == Palenzuela1D_entropy) ) {

    printf("Entropy routine failed. Retrying with energy.\n");

    // If the entropy con2prim failed, then try again
    // using the standard Palenzuela1D con2prim
    //    stats.which_routine = Palenzuela1D;

    stats.which_routine = Palenzuela1D;
    
    palenzuela( eos, S_squared,BdotS,B_squared, con, prim, SU, tolerance, stats );


  }

  return stats.c2p_failed;

}

/*****************************************************************************/
/*********************** PALENZUELA CON2PRIM FUNCTIONS ***********************/
/*****************************************************************************/
void palenzuela( const igm_eos_parameters eos,
                 const double S_squared,
                 const double BdotS,
                 const double B_squared,
                 const double * con,
                 double *restrict prim,
                 const CCTK_REAL *restrict SU,
                 const double tol_x,
                 output_stats& stats ) {

  // main root finding routine
  // using scheme of Palenzuela et al. 2015, PRD, 92, 044045

  // some quantities computed from the conservatives,
  // defined as order unity quantities
  double param[7];
  param[par_q] = con[TAU]/con[DD];
  param[par_r] = S_squared/(con[DD]*con[DD]);
  param[par_s] = B_squared/con[DD];
  param[par_t] = BdotS/(pow(con[DD],1.5));
  param[conDD] = con[DD];
  param[conYE] = con[YE];
  param[conWS] = con[WS];

  // q = tau/D
  // s = B^{2}/D
  //
  // x := hW^{2} + B^{2} - P - 0.5*( (B.v)^{2} + (B/W)^{2} )
  //
  //

  //fprintf(stderr,
  //	   ">>>>>>>>>> ENTERED con2prim_Palenzuela1D <<<<<<<<<<\n");
  //fflush(stderr);
  
  // bracket for x
  double xlow = 1.0+param[par_q]-param[par_s];
  double xup  = 2.0+2.0*param[par_q]-param[par_s];
  // initial guess for temperature
  double temp_guess=prim[TEMP];


  if( stats.which_routine == None ) {
    const int xlow_entropy_key =
        check_depsdT_condition(eos,param,temp_guess,xlow);

    const int xup_entropy_key =
        check_depsdT_condition(eos,param,temp_guess,xup);

    if( (xlow_entropy_key == Palenzuela1D_entropy) ||
        (xup_entropy_key  == Palenzuela1D_entropy) ) {

      /*static int entprint = 0;
        if(entprint < 100) {
            printf("USING ENTROPY RECOVERY\n");
            entprint++;
        }*/

        stats.which_routine = Palenzuela1D_entropy;
    }
    else {

        stats.which_routine = Palenzuela1D;
    }

    //--
    static int c2p_selection_count = 0;
    if (c2p_selection_count < 500) {
      printf("C2P_SELECTION: xlow=%e xup=%e xlow_key=%d xup_key=%d selected=%d\n",
	     xlow, xup,
	     xlow_entropy_key, xup_entropy_key,
	     stats.which_routine);
      c2p_selection_count++;
    }
    //--
    
  }

  
  // find x, this is the recovery process
  double x = zbrent(*func_root, eos, param, &temp_guess, xlow, xup, tol_x, stats);

  // calculate final set of primitives
  calc_prim(eos, x, con, param, temp_guess, prim, S_squared, BdotS, B_squared, SU, stats);

  // Execute dual-energy diagnostics (evaluates BOTH branches independently)
  dual_energy_compare(eos, x, con, param, temp_guess, prim, stats);

}


static void recover_from_energy(
    const igm_eos_parameters eos,
    double rho,
    double Ye,
    double eps,
    double *press,
    double *ent,
    double *temp)
{
    // Save the values before EOS table-bound enforcement
    const double rho_in = rho;
    const double Ye_in  = Ye;
    const double eps_in = eps;

    enforce_table_bounds_rho_Ye_eps(
        eos,
        &rho,
        &Ye,
        &eps);

    // Diagnostic: did the EOS table bounds modify the state?
    /*if (rho != rho_in || Ye != Ye_in || eps != eps_in) {
        fprintf(stderr,
                "ENERGY_BOUNDS: "
                "rho=%e -> %e "
                "Ye=%e -> %e "
                "eps=%e -> %e\n",
                rho_in, rho,
                Ye_in,  Ye,
                eps_in, eps);
		}*/

    WVU_EOS_P_S_and_T_from_rho_Ye_eps(
        rho,
        Ye,
        eps,
        press,
        ent,
        temp);
}


static bool use_entropy_final(
    const igm_eos_parameters eos,
    const double rho,
    const double temp_energy,
    const double depsdT,
    const double rel_T,
    const double rel_P,
    const double rel_eps)
{
    // Never switch if entropy evolution is not active
    if (!eos.evolve_entropy) return false;

    // --- Criterion 1: energy branch is pinned at the temperature floor ---
    const bool energy_at_floor      = (temp_energy <= 1.01 * eos.T_min);

    // --- Criterion 2: EOS is weakly temperature-sensitive ---
    const bool weak_depsdT =  (depsdT < eos.depsdT_threshold);

    // --- Criterion 3: the two branches strongly disagree in T or eps ---
    const bool large_T_discrepancy   = (rel_T   > 0.1);
    const bool large_eps_discrepancy = (rel_eps  > 0.1);

    // Density gate removed: HOT2 diagnostic log showed severe energy/entropy
    // disagreement (T_E up to ~2.3, T_S at/near floor) at rho as high as
    // ~4e-4, four orders of magnitude above the old rho_switch = 1e-6.
    // The disagreement magnitude now drives the decision directly,
    // with no assumed density band.

    const bool energy_unreliable =
      energy_at_floor ||
      large_T_discrepancy ||
      large_eps_discrepancy ||
      weak_depsdT;

    const bool switch_to_entropy = energy_unreliable;   // no density gate here

    return switch_to_entropy;    

}



static void dual_energy_compare(
    const igm_eos_parameters eos,
    const double x,
    const double *restrict con,
    const double *restrict param,
    const double temp_guess,
    double *restrict prim,
    output_stats& stats)
{
    // If c2p failed during zbrent, do not run on garbage state
    if (stats.c2p_failed || !eos.evolve_entropy) return;


    static int build_stamp = 0;
    if (build_stamp == 0) { fprintf(stderr, "BUILD_MARKER: v_B_no_density_gate\n"); build_stamp = 1; }
    
    const double rho = prim[RHO];
    const double Ye  = prim[YE];
    const double W   = prim[WLORENTZ];

    //--------------------------------------------------
    // Energy branch: reconstruct eps from mechanical x
    //--------------------------------------------------
    const double q = param[par_q];
    const double s = param[par_s];
    const double t = param[par_t];

    double eps_energy =
        W - 1.0
      + (1.0 - W * W) * x / W
      + W * (q - s + t * t / (2.0 * x * x) + s / (2.0 * W * W));

    double press_energy = 0.0;
    double ent_energy   = 0.0;
    double temp_energy  = prim[TEMP];

    recover_from_energy(eos, rho, Ye, eps_energy,
                        &press_energy, &ent_energy, &temp_energy);


    // Recompute depsdT at the converged state for use in the selector
    double depsdT_local = 0.0;
    {
      double press_tmp = 0.0;
      double ent_tmp   = 0.0;
      double temp_tmp  = temp_energy;
      double rho_tmp   = rho;
      double Ye_tmp    = Ye;
      double eps_tmp   = eps_energy;
      enforce_table_bounds_rho_Ye_eps(eos, &rho_tmp, &Ye_tmp, &eps_tmp);
      WVU_EOS_P_S_T_and_depsdT_from_rho_Ye_eps(
					       rho_tmp, Ye_tmp, eps_tmp,
					       &press_tmp, &ent_tmp, &depsdT_local, &temp_tmp);
    }
    
    
    //--------------------------------------------------
    // Entropy branch: reconstruct from conserved S_star
    //--------------------------------------------------
    double ent_entropy   = con[WS] / W;
    double press_entropy = 0.0;
    double eps_entropy   = 0.0;
    double temp_entropy  = prim[TEMP];

    recover_from_entropy(eos, rho, Ye, ent_entropy,
                         &press_entropy, &eps_entropy, &temp_entropy);

    //--------------------------------------------------
    // Compute relative differences between branches
    //--------------------------------------------------
    const double dT   = temp_energy  - temp_entropy;
    const double deps = eps_energy   - eps_entropy;
    const double dS   = ent_energy   - ent_entropy;
    const double dP   = press_energy - press_entropy;

    const double rel_eps = fabs(deps) / (fabs(eps_energy) + 1e-15);
    const double rel_T   = fabs(dT)   / (fabs(temp_energy) + 1e-15);
    const double rel_P   = fabs(dP)   / (fabs(press_energy) + 1e-15);
    const double rel_S   = fabs(dS)   / (fabs(ent_energy) + 1e-15);

    //--------------------------------------------------
    // Decide which branch wins
    //--------------------------------------------------

    const bool entropy_reliability = use_entropy_final(
       eos,
       rho,
       temp_energy,
       depsdT_local,
       rel_T,
       rel_P,
       rel_eps);


    // Surface-only hot-cell test.
    // Keep this restricted to the outer layers so it does not override
    // the bulk where the entropy branch can be too aggressive.

    const double rho_switch = 1e-4;   // start here based on your plots
    const bool surface_band  = (rho < rho_switch);
    
    const double T_hot_cut  = 1e-2;
    const double T_cold_cut = 5.0 * eos.T_min;
    
    const bool energy_hot   = (temp_energy > T_hot_cut);
    const bool entropy_cold = (temp_entropy <= T_cold_cut);
    const bool hot_suspicious = surface_band && energy_hot && entropy_cold;
    
    // If both are hot and close, trust the energy branch
    const bool both_hot = (temp_energy > T_hot_cut) && (temp_entropy > T_hot_cut);
    const bool both_hot_and_close = both_hot && (rel_T < 0.1) && (rel_eps < 0.1);
    
    // Final decision
    const bool entropy_wins =
      !both_hot_and_close &&
      (entropy_reliability || hot_suspicious);
    
    //stats.which_routine = entropy_wins ? Palenzuela1D_entropy : Palenzuela1D;
    
    static int hot_diag_count = 0;
    if (hot_diag_count < 2000 && rho > 1e-10) {
      fprintf(stderr,
	      "HOTCLASS: rho=%e T_E=%e T_S=%e sel=%s "
	      "energy_hot=%d entropy_cold=%d rel_T=%e rel_eps=%e depsdT=%e\n",
	      rho,
	      temp_energy,
	      temp_entropy,
	      entropy_wins ? "S" : "E",
	      (int)energy_hot,
	      (int)entropy_cold,
	      rel_T,
	      rel_eps,
	      depsdT_local);
      hot_diag_count++;
    }
    
    // Surface-band diagnostic: isolate the density range where the hot rim
    // sits, and record BOTH branches plus the trigger flags, so we can tell
    // apart (a) entropy branch also contaminated, (b) switch not firing,
    // vs (c) something downstream overwriting prim[TEMP] later.
    static int surface_diag_count = 0;

    if (surface_diag_count < 20000 && rho > 1e-7 && rho < 1e-2) {
      fprintf(stderr,
	      "SURF: rho=%e T_E=%e T_S=%e sel=%s rel_T=%e rel_eps=%e "
	      "weakdepsdT=%d floor=%d\n",
	      rho, temp_energy, temp_entropy,
	      entropy_wins ? "S" : "E",
	      rel_T, rel_eps,
	      (int)(depsdT_local < eos.depsdT_threshold),
	      (int)(temp_energy <= 1.01 * eos.T_min));
      surface_diag_count++;
    }
    
    
    //--------------------------------------------------
    // Apply the selected branch to prim[]
    // Only thermal fields are overwritten.
    // Mechanical fields (RHO, WLORENTZ, velocities, B) are untouched.
    //--------------------------------------------------
    if (entropy_wins) {
        prim[EPS]   = eps_entropy;
        prim[TEMP]  = temp_entropy;
        prim[PRESS] = press_entropy;
        prim[ENT]   = ent_entropy;
    } else {
        prim[EPS]   = eps_energy;
        prim[TEMP]  = temp_energy;
        prim[PRESS] = press_energy;
        prim[ENT]   = ent_energy;
    }

    //--------------------------------------------------
    // Diagnostics
    //--------------------------------------------------
    /*if (rho > 1e-8 && (rel_eps > 1e-6 || rel_S > 1e-6 || rel_T > 1e-3)) {
        fprintf(stderr,
                "THERMO_DIFF: rho=%e W=%e dT=%e deps=%e dS=%e dP=%e\n",
                rho, W, dT, deps, dS, dP);
        fprintf(stderr,
                "REL_DIFF: eps=%e T=%e P=%e S=%e\n",
                rel_eps, rel_T, rel_P, rel_S);
    }

     fprintf(stderr,
            "DUAL_STATE: rho=%e W=%e Ye=%e rout=%d selected=%s\n"
            "  ENERGY : eps=%e T=%e P=%e S=%e\n"
            "  ENTROPY: eps=%e T=%e P=%e S=%e\n",
            rho, W, Ye, (int)stats.which_routine,
            entropy_wins ? "entropy" : "energy",
            eps_energy, temp_energy, press_energy, ent_energy,
            eps_entropy, temp_entropy, press_entropy, ent_entropy);*/

}


static void recover_from_entropy(
    const igm_eos_parameters eos,
    double rho,
    double Ye,
    double ent,
    double *press,
    double *eps,
    double *temp)
{
    // Save the values before EOS table-bound enforcement
    const double rho_in = rho;
    const double Ye_in  = Ye;
    const double ent_in = ent;

    enforce_table_bounds_rho_Ye_S(
        eos,
        &rho,
        &Ye,
        &ent);

    // Diagnostic: did the EOS table bounds modify the state?
    /*    if (rho != rho_in || Ye != Ye_in || ent != ent_in) {
        fprintf(stderr,
                "ENTROPY_BOUNDS: "
                "rho=%e -> %e "
                "Ye=%e -> %e "
                "S=%e -> %e\n",
                rho_in, rho,
                Ye_in,  Ye,
                ent_in, ent);
		}*/

    WVU_EOS_P_eps_and_T_from_rho_Ye_S(
        rho,
        Ye,
        ent,
        press,
        eps,
        temp);
}


void calc_prim( const igm_eos_parameters eos,
                const double x,
                const double *restrict con,
                const double *restrict param,
                const double temp_guess,
                double *restrict prim,
                const double S_squared,
                const double BdotS,
                const double B_squared,
                const double *restrict SU,
                output_stats& stats ) {

  const double q = param[par_q];
  const double r = param[par_r];
  const double s = param[par_s];
  const double t = param[par_t];

  double Wminus2 =
      1.0 - ( x*x*r + (2*x+s)*t*t )
      / ( x*x*(x+s)*(x+s) );

  Wminus2 = fmin(fmax(Wminus2,
                      eos.inv_W_max_squared),
                 1.0);

  const double W = pow(Wminus2, -0.5);

  double rho = con[DD]/W;
  double ye  = con[YE]/con[DD];

  double temp = temp_guess;
  double press = 0.0;
  double eps = 0.0;
  double ent = 0.0;

  if( eos.evolve_T ) {

    if (stats.which_routine == Palenzuela1D) {

      eps =
          W - 1.0
          + (1.0-W*W)*x/W
          + W*(q-s+t*t/(2*x*x)+s/(2*W*W));

      enforce_table_bounds_rho_Ye_eps(
          eos,&rho,&ye,&eps);

      WVU_EOS_P_S_and_T_from_rho_Ye_eps(
          rho,ye,eps,
          &press,&ent,&temp);

    }
    else {

      ent = con[WS]/W;

      enforce_table_bounds_rho_Ye_S(
          eos,&rho,&ye,&ent);

      WVU_EOS_P_eps_and_T_from_rho_Ye_S(
          rho,ye,ent,
          &press,&eps,&temp);
    }

  }
  else {

    enforce_table_bounds_rho_Ye_T(
        eos,&rho,&ye,&temp);

    WVU_EOS_P_eps_and_S_from_rho_Ye_T(
        rho,ye,temp,
        &press,&eps,&ent);
  }

  const double Z = x*rho*W;

  prim[UTCON1] =
      W*(SU[0] + (BdotS)*con[B1_con]/Z)
      /(Z+B_squared);

  prim[UTCON2] =
      W*(SU[1] + (BdotS)*con[B2_con]/Z)
      /(Z+B_squared);

  prim[UTCON3] =
      W*(SU[2] + (BdotS)*con[B3_con]/Z)
      /(Z+B_squared);

  prim[RHO]      = rho;
  prim[EPS]      = eps;
  prim[B1_con]   = con[B1_con];
  prim[B2_con]   = con[B2_con];
  prim[B3_con]   = con[B3_con];
  prim[TEMP]     = temp;
  prim[YE]       = ye;
  prim[PRESS]     = press;
  prim[WLORENTZ] = W;
  prim[ENT]      = ent;
}

  
double func_root( const igm_eos_parameters eos,
                  double x,
                  double *restrict param,
                  double *restrict temp_guess,
                  output_stats& stats ) {

  // computes f(x) from x and q,r,s,t

  const double q = param[par_q];
  const double r = param[par_r];
  const double s = param[par_s];
  const double t = param[par_t];

  // (i)
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));
  Wminus2        = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  const double W = pow(Wminus2, -0.5);

  // (ii)
  double rho     = param[conDD]/W;
  double ye      = param[conYE]/param[conDD];
  double temp    = *temp_guess;
  double P       = 0.0;
  double eps     = 0.0;
  double ent     = 0.0;

  if( eos.evolve_T ) {

    //--debug
    //static int rtprint=0; if(rtprint<20) {printf("FINAL: stats.which_routine = %d\n", stats.which_routine);rtprint++;}
    //
    
    if( stats.which_routine == Palenzuela1D ) {

      //-- debug
      //static int energy1=0; if(energy1<20) { printf("FINAL EOS: ENERGY branch\n");energy1++;}
      //-
      
      // Default Palenzuela con2prim: get Hydro quantities from (rho,Ye,eps)
      eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );

      //--
      //static int finaleps=0;if(finaleps<20) {printf("FINAL ENERGY INPUT: rho=%e Ye=%e eps=%e\n",rho,ye,eps);finaleps++;}
      //--
      
      enforce_table_bounds_rho_Ye_eps( eos,&rho,&ye,&eps );

      //--
      //static int finaleps2=0;if(finaleps2<20) {printf("FINAL ENERGY AFTER FLOOR: rho=%e Ye=%e eps=%e\n",rho,ye,eps);finaleps2++;}
      //--
      
      WVU_EOS_P_S_and_T_from_rho_Ye_eps( rho,ye,eps, &P,&ent,&temp );
      //--
      //static int finalout=0;if(finalout<20) {printf("FINAL ENERGY OUTPUT: P=%e S=%e T=%e\n",P,ent,temp);finalout++;}
      //--
      
    }
    else {
      //--
      //static int first_entropy_exec = 0; if(first_entropy_exec < 50) {printf("FINAL EOS: EXECUTING ENTROPY rho=%e Ye=%e\n", rho, ye);first_entropy_exec++;}
      //--
      // Modified Palenzuela con2prim: get Hydro quantities from (rho,Ye,S)
      ent = param[conWS]/W;
      enforce_table_bounds_rho_Ye_S( eos,&rho,&ye,&ent );

      //--
      //static int finalS2=0;if(finalS2<20) {printf("FINAL ENTROPY AFTER FLOOR: rho=%e Ye=%e S=%e\n",rho,ye,ent);finalS2++;}
      //--
      
      //--
      //static int finalS=0;if(finalS<20) {printf("FINAL ENTROPY INPUT: rho=%e Ye=%e S=%e conWS=%e W=%e\n",rho,ye,ent,param[conWS],W);finalS++;}
      //--
      
      WVU_EOS_P_eps_and_T_from_rho_Ye_S( rho,ye,ent, &P,&eps,&temp );
      //--
      //static int finalSout=0;if(finalSout<20) {printf("FINAL ENTROPY OUTPUT: P=%e eps=%e T=%e\n",P,eps,temp);finalSout++;}
      //--
    }
  }
  else {
    enforce_table_bounds_rho_Ye_T( eos,&rho,&ye,&temp );
    WVU_EOS_P_eps_and_S_from_rho_Ye_T( rho,ye,temp, &P,&eps,&ent );
  }

  // (iv)
  double ans = x - (1.0 + eps + P/rho)*W;

  return ans;
}

/****************************************************************************/
/************************* BRENT's METHOD FUNCTIONS *************************/
/****************************************************************************/

#define ITMAX 300   //Maximum allowed number of iterations.
#define PREC 3.0e-15
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double zbrent( double (*func)(const igm_eos_parameters, double, double *, double *restrict, output_stats&),
               const igm_eos_parameters eos,
               double *restrict param,
               double * temp_guess,
               double x1,
               double x2,
               double tol_x,
               output_stats& stats ) {

  // Implementation of Brent’s method to find the root of a function func to accuracy tol_x. The root
  // must be bracketed by x1 and x2.
  // Based on W. H. Press et al. 2007, Numerical Recipes in C: The Art of Scientific Computing, 3rd edition,
  // Cambridge University Press, ISBN 978-0-521-88068-8

  int iter = 0;

  //termination critera definitions
  // int MAX_NEWT_ITER   = c2p.max_iterations;
  // int EXTRA_NEWT_ITER = c2p.extra_iterations;
  double ans = 0.0;

  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1=0.0,min2=0.0;
  double temp_guess_a = *temp_guess;
  double temp_guess_b = *temp_guess;
  double temp_guess_c = *temp_guess;
  double fa=(*func)(eos, a, param, &temp_guess_a, stats);
  double fb=(*func)(eos, b, param, &temp_guess_b, stats);
  double fc,p,q,r,s,tol1;

  // root must be bracketed
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    stats.c2p_failed = true;
    return b;
  }
  fc=fb;
  temp_guess_c = temp_guess_b;

  int i_extra = 0;
  int doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating = 1;
  double maxerror = 0.0;

  while (keep_iterating) {

    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      // reset c, adjust bracket interval d
      c=a;
      fc=fa;
      e=d=b-a;
      temp_guess_c = temp_guess_a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
      temp_guess_a = temp_guess_b;
      temp_guess_b = temp_guess_c;
      temp_guess_c = temp_guess_a;
    }

    tol1=2.0*PREC*fabs(b)+0.5*tol_x;
    maxerror=0.5*(c-b);
    // check convergence
    if (fabs(maxerror) <= tol1 || fb == 0.0){
      // we are done
      stats.c2p_failed = false;
      keep_iterating = 0;
      ans = b;
    }
    else if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      // inverse quadratic interpolation
      s=fb/fa;

      if (a == c) {
        // fake quad interpolation
        p=2.0*maxerror*s;
        q=1.0-s;
      } else {
        // actual quad interpolation
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*maxerror*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }

      if (p > 0.0) q = -q;  // check whether in bounds
      p=fabs(p);
      min1=3.0*maxerror*q-fabs(tol1*q);
      min2=fabs(e*q);

      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        // interpolation successful
        e=d;
        d=p/q;
      } else {
        // interpolation failed, use bisection
        d=maxerror;
        e=d;
      }

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1,maxerror);
      fb=(*func)(eos, b, param, &temp_guess_b, stats);


    } else {

      // Bounds decreasing too slowly, use bisection

      d=maxerror;
      e=d;

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1, maxerror);
      fb=(*func)(eos, b, param, &temp_guess_b, stats);
    }

    ++iter;

    // termination criterion
    if( (fabs(maxerror) <= tol1 ) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(maxerror) <= tol1)&&(doing_extra == 0))
        || (i_extra > EXTRA_NEWT_ITER) || (iter >= (ITMAX)) ) {
      keep_iterating = 0;
    }

  } //end for loop

  if( (!robust_isfinite(fb)) ) {
    stats.c2p_failed = true;
    return b;
  } if( (fabs(maxerror) <= tol1 || fb == 0.0)){
    stats.c2p_failed = false;
    return ans;
  } else if( (fabs(maxerror) <= tol1) && (fabs(maxerror) > tol1) ){
    stats.c2p_failed = false;
    return ans;
  } else {
    stats.c2p_failed = true;
    return b;
  }

}
