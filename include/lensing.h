/** @file lensing.h Documented includes for harmonic module */

#ifndef __LENSING__
#define __LENSING__

#include "harmonic.h"

/**
 * various types of lensing methods
 * The digit behind q is the order of C_2 (quadrupole of lensing correlation)
 * The order behind m is the order of C_0 (monopole of lensing correlation)
 * All sky is the method of 0502425. Simplified all sky is the method of 0601594
 */


/**
 * Structure containing everything about lensed spectra that other modules need to know.
 *
 * Once initialized by lensing_init(), contains a table of all lensed
 * \f$ C_l\f$'s for the all modes (scalar/tensor), all types (TT, TE...),
 * and all pairs of initial conditions (adiabatic, isocurvatures...).
 * FOR THE MOMENT, ASSUME ONLY SCALAR & ADIABATIC
 */

struct lensing {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  short has_lensed_cls; /**< do we need to compute lensed \f$ C_l\f$'s at all ? */

  //@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

  int has_tt; /**< do we want lensed \f$ C_l^{TT}\f$? (T = temperature) */
  int has_ee; /**< do we want lensed \f$ C_l^{EE}\f$? (E = E-polarization) */
  int has_te; /**< do we want lensed \f$ C_l^{TE}\f$? */
  int has_bb; /**< do we want \f$ C_l^{BB}\f$? (B = B-polarization) */
  int has_pp; /**< do we want \f$ C_l^{\phi\phi}\f$? (\f$ \phi \f$ = CMB lensing potential) */
  int has_tp; /**< do we want \f$ C_l^{T\phi}\f$? */
  int has_dd; /**< do we want \f$ C_l^{dd}\f$? (d = matter density) */
  int has_td; /**< do we want \f$ C_l^{Td}\f$? */
  int has_ll; /**< do we want \f$ C_l^{ll}\f$? (l = lensing potential) */
  int has_tl; /**< do we want \f$ C_l^{Tl}\f$? */

  int index_lt_tt; /**< index for type \f$ C_l^{TT} \f$*/
  int index_lt_ee; /**< index for type \f$ C_l^{EE} \f$*/
  int index_lt_te; /**< index for type \f$ C_l^{TE} \f$*/
  int index_lt_bb; /**< index for type \f$ C_l^{BB} \f$*/
  int index_lt_pp; /**< index for type \f$ C_l^{\phi\phi} \f$*/
  int index_lt_tp; /**< index for type \f$ C_l^{T\phi} \f$*/
  int index_lt_dd; /**< index for type \f$ C_l^{dd} \f$*/
  int index_lt_td; /**< index for type \f$ C_l^{Td} \f$*/
  int index_lt_ll; /**< index for type \f$ C_l^{dd} \f$*/
  int index_lt_tl; /**< index for type \f$ C_l^{Td} \f$*/

  int lt_size; /**< number of \f$ C_l\f$ types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{

  int l_unlensed_max;    /**< last multipole in all calculations (same as in harmonic module)*/

  int l_lensed_max;    /**< last multipole at which lensed spectra are computed */

  /* interpolable version: */

  int l_size;       /**< number of l values */

  int * l_max_lt;    /**< last multipole (given as an input) at which
                        we want to output \f$ C_l \f$'s for a given mode and type */

  double * l;       /**< table of multipole values l[index_l] */
  double * cl_lens; /**< table of anisotropy spectra for each
                       multipole and types,
                       cl[index_l * ple->lt_size + index_lt] */

  double * ddcl_lens; /**< second derivatives for interpolation */

  //@}

  /** @name - technical parameters */

  //@{

  short lensing_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  short is_allocated; /**< flag is set to true if allocated */

  int lensing_C2_order ; /**< lensing order of the quadrupolar contribution to the deflection correlation (called C^lg_2 in some papers) */

  int lensing_C0_order ; /**< lensing order of the quadrupolar contribution to the deflection correlation (called C^lg_0 in some papers) */

  //@}
};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int lensing_cl_at_l(
                      struct lensing * ple,
                      int l,
                      double * cl_lensed
                      );

  int lensing_init(
                   struct precision * ppr,
                   struct perturbations * ppt,
                   struct harmonic * phr,
                   struct fourier * pfo,
                   struct lensing * ple
                   );

  int lensing_free(
                   struct lensing * ple
                   );

  int lensing_indices(
                      struct precision * ppr,
                      struct harmonic * phr,
                      struct lensing * ple
                      );

  int lensing_lensed_cl_tt(
                           double *ksi,
                           double **d00,
                           double *w8,
                           int nmu,
                           struct lensing * ple
                           );

  int lensing_lensed_cl_te(
                           double *ksiX,
                           double **d20,
                           double *w8,
                           int nmu,
                           struct lensing * ple
                           );

  int lensing_lensed_cl_ee_bb(
                              double *ksip,
                              double *ksim,
                              double **d22,
                              double **d2m2,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              );
  int lensing_addback_cl_tt(
                            struct lensing *ple,
                            double *cl_tt
                            );

  int lensing_addback_cl_te(
                            struct lensing *ple,
                            double *cl_te
                            );

  int lensing_addback_cl_ee_bb(
                               struct lensing *ple,
                               double *cl_ee,
                               double *cl_bb
                               );


  /**int lensing_X000(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X000
                   );

  int lensing_Xp000(
                    double * mu,
                    int num_mu,
                    int lmax,
                    double * sigma2,
                    double ** Xp000
                    );

  int lensing_X220(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X220
                   );

  int lensing_X022(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X022
                   );

  int lensing_Xp022(
                    double * mu,
                    int num_mu,
                    int lmax,
                    double * sigma2,
                    double ** Xp022
                    );

  int lensing_X121(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X121
                   );

  int lensing_X132(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X132
                   );

  int lensing_X242(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X242
                   );*/

  int lensing_dm1m2(
                   double * mu,
                   int num_mu,
                   int lmax,
		   int m1,
		   int m2,
                   double ** dm1m2,
		   ErrorMsg erreur
                   );
  
#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
