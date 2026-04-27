/** @file magnetic.h Documented includes for trg module */

#include "primordial.h"
#include "trigonometric_integrals.h"

#ifndef __MAGNETIC__
#define __MAGNETIC__

enum Bk_outputs {Bk_linear,Bk_nonlinear};

/**
 * Structure containing all information TODO document
 *
 */

struct magnetic {

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int index_md_vectors; /**< set equal to phr->index_md_vectors
                           (useful since this module only deals with
                           vectors) */
  int ic_size;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int ic_ic_size;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short * is_non_zero; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */

  //@}

  //@{
  short has_Bk; /**< do we need magnetic Fourier spectrum? */

  int k_size;      /**< k_size = total number of k values */
  int k_size_pk;   /**< k_size = number of k values for P(k,z) and T(k,z) output) */
  double * k;      /**< k[index_k] = list of k values */
  double * ln_k;   /**< ln_k[index_k] = list of log(k) values */

  double * ln_tau;     /**< log(tau) array, only needed if user wants
                          some output at z>0, instead of only z=0.  This
                          array only covers late times, used for the
                          output of P(k) or T(k), and matching the
                          condition z(tau) < z_max_pk */

  int ln_tau_size;     /**< total number of values in this array */

  double * ln_Bk_ic_l;   /**< Magnetic power spectrum (linear).
                             Depends on indices index_ic1_ic2, index_k, index_tau as:
                             ln_Bk_ic_l[(index_tau * pfo->k_size + index_k)* pfo->ic_ic_size + index_ic1_ic2]
			     where index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since
                             the primordial spectrum is symmetric in (index_ic1, index_ic2)).
                             - for diagonal elements (index_ic1 = index_ic2) this arrays contains
                             ln[B(k)] where B(k) is positive by construction.
                             - for non-diagonal elements this arrays contains the k-dependent
                             cosine of the correlation angle, namely
                             B(k)_(index_ic1, index_ic2)/sqrt[B(k)_index_ic1 B(k)_index_ic2]
                             This choice is convenient since the sign of the non-diagonal cross-correlation
                             could be negative. For fully correlated or anti-correlated initial conditions,
                             this non-diagonal element is independent on k, and equal to +1 or -1.
                          */

  double * ddln_Bk_ic_l; /**< second derivative of above array with respect to log(tau), for spline interpolation. So:
                             - for index_ic1 = index_ic, we spline ln[B(k)] vs. ln(k), which is
                             good since this function is usually smooth.
                             - for non-diagonal coefficients, we spline
                             B(k)_(index_ic1, index_ic2)/sqrt[B(k)_index_ic1 B(k)_index_ic2]
                             vs. ln(k), which is fine since this quantity is often assumed to be
                             constant (e.g for fully correlated/anticorrelated initial conditions)
                             or nearly constant, and with arbitrary sign.
                          */

  double * ln_Bk_l;   /**< Total magnetic power spectrum summed over initial conditions (linear).
                          Only depends on indices index_k, index_tau as:
                          ln_Bk[index_tau * pfo->k_size + index_k]
                       */

  double * ddln_Bk_l; /**< second derivative of above array with respect to log(tau), for spline interpolation. */

  double sigma1;   /**< sigma1 */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int k_size_extra;/** total number of k values of extrapolated k array (high k)*/

  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values, covering
                      all the values of the perturbation module */

  //@}

  /** @name - parameters for the pk_eq method */
  /** @name - technical parameters */

  //@{

  short magnetic_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  short is_allocated; /**< flag is set to true if allocated */
  
  //@}
};

/**
 * Structure containing variables used only internally in fourier module by various functions.
 *
 */

/*struct magnetic_workspace {

  //@{


  //@}

};*/

/********************************************************************************/

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* external functions (meant to be called from other modules) */

  int magnetic_Bk_at_z(
                      struct background * pba,
                      struct magnetic *pma,
                      enum linear_or_logarithmic mode,
                      double z,
                      double * out_Bk,
                      double * out_Bk_ic
                      );
  
  int magnetic_sigmas_at_z(
                          struct precision * ppr,
                          struct background * pba,
                          struct magnetic * pma,
                          double R,
                          double z,
			  //                          enum out_sigmas sigma_output,
                          double * result
                          );

  
  /* internal functions */

  int magnetic_init(
                   struct precision *ppr,
                   struct background *pba,
                   struct thermodynamics *pth,
                   struct perturbations *ppt,
                   struct primordial *ppm,
                   struct magnetic *pma
                   );

  int magnetic_free(
                   struct magnetic *pma
                   );

  int magnetic_indices(
                      struct precision *ppr,
                      struct background *pba,
                      struct perturbations * ppt,
                      struct primordial * ppm,
                      struct magnetic * pma
                      );

  int magnetic_get_k_list(
                         struct precision *ppr,
                         struct perturbations * ppt,
                         struct magnetic * pma
                         );

  int magnetic_get_tau_list(
                           struct perturbations * ppt,
                           struct magnetic * pma
                           );


  int magnetic_Bk_linear(
                        struct background *pba,
                        struct perturbations *ppt,
                        struct primordial *ppm,
                        struct magnetic *pma,
			int index_tau,
                        int k_size,
                        double * lnBk,
                        double * lnBk_ic
                        );

  int magnetic_sigmas(
                     struct magnetic * pma,
                     double R,
                     double *lnpk_l,
                     double *ddlnpk_l,
                     int k_size,
                     double k_per_decade,
		     //                     enum out_sigmas sigma_output,
                     double * result
                     );

  
  int magnetic_get_source(
			  struct background * pba,
			  struct perturbations * ppt,
			  struct magnetic * pma,
			  int index_k,
			  int index_ic,
			  int index_tp,
			  int index_tau,
			  double ** sources,
			  double * source);
  

  
#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
