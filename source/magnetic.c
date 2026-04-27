/** @file magnetic.c Documented magnetic power sectrum module
 *
 * Cyril Pitrou  june 2023
 *
 */

#include "magnetic.h"

/**
 * Initialize the magnetic structure
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to therodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pma Input/Output: pointer to initialized magnetic structure
 * @return the error status
 */

int magnetic_init(
                 struct precision *ppr,
                 struct background *pba,
                 struct thermodynamics *pth,
                 struct perturbations *ppt,
                 struct primordial *ppm,
                 struct magnetic *pma
                 ) {

  int index_k;
  int index_tau;
  int index_tau_sources;
  int index_tau_late;

  //  double **pk_nl;
  //  double **lnpk_l;
  //  double **ddlnpk_l;

  /** - preliminary tests */
  pma->has_Bk = ppt->has_magnetic_transfer;
  
  /** --> This module only makes sense for dealing with scalar
      perturbations, so it should do nothing if there are no
      scalars */
  if (ppt->has_vectors == _FALSE_) {
    if (pma->magnetic_verbose > 0)
      printf("No vectors modes requested. Magnetic module skipped.\n");
    return _SUCCESS_;
  }

  /** --> Nothing to be done if we don't want the matter power spectrum */

  if (pma->has_Bk == _FALSE_) {
    if (pma->magnetic_verbose > 0)
      printf("No magnetic transfer requested. Magnetic module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pma->magnetic_verbose > 0)
      printf("Computing linear magnetic spectra.\n");
  }

  /** - define indices in fourier structure (and allocate some arrays in the structure) */

  class_call(magnetic_indices(
                             ppr,
                             pba,
                             ppt,
                             ppm,
                             pma),
             pma->error_message,
             pma->error_message);

  for (index_tau=0; index_tau<pma->ln_tau_size;index_tau++) {

    index_tau_sources = ppt->tau_size-ppt->ln_tau_size+index_tau;

    class_call(magnetic_Bk_linear(
				  pba,
				  ppt,
				  ppm,
				  pma,
				  index_tau_sources,
				  pma->k_size,
				  &(pma->ln_Bk_l[index_tau * pma->k_size]),
				  &(pma->ln_Bk_ic_l[index_tau * pma->k_size * pma->ic_ic_size])
				  ),
	       pma->error_message,
	       pma->error_message);
    
	       }

    class_call(magnetic_sigmas_at_z(ppr,
                                   pba,
                                   pma,
                                   1./pba->h,
                                   0.,
				    //                                   out_sigma,
				    &(pma->sigma1)),
               pma->error_message,
               pma->error_message);
    

  if (pma->magnetic_verbose>0) {

    fprintf(stdout," -> sigma_1(B)=%g Gauss (computed till k = %g h/Mpc)\n",
              pma->sigma1,
              pma->k[pma->k_size-1]/pba->h);
    
  }
  

  //    free(pk_nl);
  //    free(lnpk_l);
  //    free(ddlnpk_l);


  pma->is_allocated = _TRUE_;
  return _SUCCESS_;
}

/**
 * Free all memory space allocated by magnetic_init().
 *
 *
 * @param pfo Input: pointer to fourier structure (to be freed)
 * @return the error status
 */

int magnetic_free(
                 struct magnetic *pma
                 ) {
  if (pma->has_Bk == _TRUE_) {
    free(pma->k);
    free(pma->ln_k);
    
    free(pma->ln_Bk_ic_l);
    free(pma->ln_Bk_l);

    if (pma->ln_tau_size>1) {
      free(pma->ddln_Bk_ic_l);
      free(pma->ddln_Bk_l);
      free(pma->ln_tau);
    }

    free(pma->is_non_zero);
  }

  pma->is_allocated = _FALSE_;
  return _SUCCESS_;
}

/**
 * Define indices in the fourier structure, and when possible, allocate
 * arrays in this structure given the index sizes found here
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pma Input/Output: pointer to magnetic structure
 * @return the error status
 */

int magnetic_indices(
                    struct precision *ppr,
                    struct background *pba,
                    struct perturbations * ppt,
                    struct primordial * ppm,
                    struct magnetic * pma
                    ) {

  int index_ic1_ic2;

  /** - define indices for initial conditions (and allocate related arrays) */
  pma->index_md_vectors = ppt->index_md_vectors;
  pma->ic_size = ppm->ic_size[pma->index_md_vectors];
  pma->ic_ic_size = ppm->ic_ic_size[pma->index_md_vectors];
  class_alloc(pma->is_non_zero,sizeof(short)*pma->ic_ic_size,pma->error_message);
  for (index_ic1_ic2=0; index_ic1_ic2 < pma->ic_ic_size; index_ic1_ic2++)
    pma->is_non_zero[index_ic1_ic2] = ppm->is_non_zero[pma->index_md_vectors][index_ic1_ic2];


  /** - get list of k values */

  class_call(magnetic_get_k_list(ppr,ppt,pma),
             pma->error_message,
             pma->error_message);

  /** - get list of tau values */

  class_call(magnetic_get_tau_list(ppt,pma),
             pma->error_message,
             pma->error_message);


  class_alloc(pma->ln_Bk_ic_l,pma->ln_tau_size*pma->k_size*pma->ic_ic_size*sizeof(double*),pma->error_message);
  class_alloc(pma->ln_Bk_l,pma->ln_tau_size*pma->k_size*sizeof(double*),pma->error_message);

  /** - if interpolation of \f$B(k,\tau)\f$ will be needed (as a function of tau),                                                                                                                                                                                                           
      compute also the array of second derivatives in view of spline interpolation */

  if (pma->ln_tau_size > 1) {

    class_alloc(pma->ddln_Bk_ic_l,pma->ln_tau_size*pma->k_size*pma->ic_ic_size*sizeof(double*),pma->error_message);
    class_alloc(pma->ddln_Bk_l,pma->ln_tau_size*pma->k_size*sizeof(double*),pma->error_message);
  }

  
  return _SUCCESS_;
}


/** Copy list of k.                                         
 *                                                                                                                                                                                                                                                                                           
 * @param ppr Input: pointer to precision structure                                                                                                                                                                                                                                          
 * @param ppt Input: pointer to perturbation structure                                                                                                                                                                                                                                       
 * @param pfo Input/Output: pointer to fourier structure                                                                                                                                                                                                                                     
 * @return the error status                                                                                                                                                                                                                                                                  
 */

int magnetic_get_k_list(
                       struct precision *ppr,
                       struct perturbations * ppt,
                       struct magnetic * pma
                       ) {

  double k=0;
  double k_max,exponent;
  int index_k;

  pma->k_size = ppt->k_size[pma->index_md_vectors];
  pma->k_size_pk = ppt->k_size_pk;
  k_max = ppt->k[pma->index_md_vectors][pma->k_size-1];

  /** - allocate array of k */
  class_alloc(pma->k,   pma->k_size*sizeof(double),pma->error_message);
  class_alloc(pma->ln_k,pma->k_size*sizeof(double),pma->error_message);

  /** - fill array of k (not extrapolated) */
  for (index_k=0; index_k<pma->k_size; index_k++) {
    k = ppt->k[pma->index_md_vectors][index_k];
    pma->k[index_k] = k;
    pma->ln_k[index_k] = log(k);
  }

  return _SUCCESS_;
}



/**                                                                                                                                                                                                                                                                                         
 * Copy list of tau from perturbation module                                                                                                                                                                                                                                                 
 *
 * @param ppt Input: pointer to perturbation structure
 * @param pma Input/Output: pointer to magnetic structure
 * @return the error status
 */

int magnetic_get_tau_list(
                         struct perturbations * ppt,
                         struct magnetic * pma
                         ) {

  int index_tau;

  /** -> for linear calculations: only late times are considered, given the value z_max_pk inferred from the input */
  pma->ln_tau_size = ppt->ln_tau_size;

  if (ppt->ln_tau_size > 1) {

    class_alloc(pma->ln_tau,pma->ln_tau_size*sizeof(double),pma->error_message);

    for (index_tau=0; index_tau<pma->ln_tau_size;index_tau++) {
      pma->ln_tau[index_tau] = ppt->ln_tau[index_tau];
    }
  }

  return _SUCCESS_;
}

/** TODO put documentation */

int magnetic_Bk_linear(
                      struct background *pba,
                      struct perturbations *ppt,
                      struct primordial *ppm,
                      struct magnetic *pma,
		      int index_tau,
                      int k_size,
                      double * lnBk,    //lnBk[index_k]
                      double * lnBk_ic  //lnBk[index_k * pfo->ic_ic_size + index_ic1_ic2]
                      ) {

  int index_k;
  int index_tp;
  int index_ic1,index_ic2,index_ic1_ic1,index_ic1_ic2,index_ic2_ic2;
  double * primordial_pk;
  double Bk;
  double * Bk_ic;
  double source_ic1;
  double source_ic2;
  double cosine_correlation;

  /** - allocate temporary vector where the primordial spectrum will be stored */

  class_alloc(primordial_pk,pma->ic_ic_size*sizeof(double),pma->error_message);

  class_alloc(Bk_ic,pma->ic_ic_size*sizeof(double),pma->error_message);

  index_tp = ppt->index_tp_magnetic;

  /** - loop over k values */

  for (index_k=0; index_k<k_size; index_k++) {

    /** --> get primordial spectrum */
    class_call(primordial_spectrum_at_k(ppm,pma->index_md_vectors,logarithmic,pma->ln_k[index_k],primordial_pk),
               ppm->error_message,
               pma->error_message);

    /** --> initialize a local variable for B(k) to zero */
    Bk = 0.;

    /** --> here we recall the relations relevant for the nomalization fo the power spectrum:
        For adiabatic modes, the curvature primordial spectrum thnat we just read was:
        P_R(k) = 1/(2pi^2) k^3 < R R >
        Thus the primordial curvature correlator is given by:
        < R R > = (2pi^2) k^-3 P_R(k)
        So the delta_m correlator reads:
        P(k) = < delta_m delta_m > = (source_m)^2 < R R > = (2pi^2) k^-3 (source_m)^2 P_R(k)

        For isocurvature or cross adiabatic-isocurvature parts,
        one would just replace one or two 'R' by 'S_i's */

    /** Here we consider vector modes and we sum over +-1 hence there is globale factor 2 to add on top of these considerations which are related to scalar stuff. However this is already included in the vector_to_scalar ratio so we do not put it */
    /** Be careful and clear about tehse conventions TODO */

    /** --> get contributions to P(k) diagonal in the initial conditions */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {

      index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pma->ic_size);

      class_call(magnetic_get_source(pba,
                                    ppt,
                                    pma,
                                    index_k,
                                    index_ic1,
                                    index_tp,
                                    index_tau,
                                    ppt->sources[pma->index_md_vectors],
                                    &source_ic1),
                 pma->error_message,
                 pma->error_message);

      Bk_ic[index_ic1_ic1] = 2.*_PI_*_PI_/exp(3.*pma->ln_k[index_k])
        *source_ic1*source_ic1
        *exp(primordial_pk[index_ic1_ic1]);

      Bk += Bk_ic[index_ic1_ic1];

      //printf("DEBUG we have added %e to Bk\n",Bk_ic[index_ic1_ic1]);

      if (lnBk_ic != NULL) {
        lnBk_ic[index_k * pma->ic_ic_size + index_ic1_ic1] = log(Bk_ic[index_ic1_ic1]);
      }
    }

    /** --> get contributions to P(k) non-diagonal in the initial conditions */
    for (index_ic1 = 0; index_ic1 < pma->ic_size; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < pma->ic_size; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);
        index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pma->ic_size);
        index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,pma->ic_size);

        if (pma->is_non_zero[index_ic1_ic2] == _TRUE_) {

          class_call(magnetic_get_source(pba,
                                        ppt,
                                        pma,
                                        index_k,
                                        index_ic1,
                                        index_tp,
                                        index_tau,
                                        ppt->sources[pma->index_md_vectors],
                                        &source_ic1),
                     pma->error_message,
                     pma->error_message);

          class_call(magnetic_get_source(pba,
                                        ppt,
                                        pma,
                                        index_k,
                                        index_ic2,
                                        index_tp,
                                        index_tau,
                                        ppt->sources[pma->index_md_vectors],
                                        &source_ic2),
                     pma->error_message,
                     pma->error_message);

          cosine_correlation = primordial_pk[index_ic1_ic2]*SIGN(source_ic1)*SIGN(source_ic2);

          Bk_ic[index_ic1_ic2] = cosine_correlation * sqrt(Bk_ic[index_ic1_ic1]*Bk_ic[index_ic2_ic2]);

          Bk += 2.*Bk_ic[index_ic1_ic2];

          if (lnBk_ic != NULL) {
            lnBk_ic[index_k * pma->ic_ic_size + index_ic1_ic2] = cosine_correlation;
          }
        }
        else {
          if (lnBk_ic != NULL) {
            lnBk_ic[index_k * pma->ic_ic_size + index_ic1_ic2] = 0.;
          }
        }
      }
    }

    lnBk[index_k] = log(Bk);
  }

  free(primordial_pk);
  free(Bk_ic);

  return _SUCCESS_;

}



/**
 * Get sources for a given wavenumber (and for a given time, type, ic,
 * mode...) either directly from precomputed valkues (computed ain
 * perturbation module), or by analytic extrapolation
 *
 * @param pba             Input: pointer to background structure
 * @param ppt             Input: pointer to perturbation structure
 * @param pfo             Input: pointer to fourier structure
 * @param index_k         Input: index of required k value
 * @param index_ic        Input: index of required ic value
 * @param index_tp        Input: index of required tp value
 * @param index_tau       Input: index of required tau value
 * @param sources         Input: array containing the original sources
 * @param source          Output: desired value of source
 * @return the error status
 */

int magnetic_get_source(
                       struct background * pba,
                       struct perturbations * ppt,
                       struct magnetic * pma,
                       int index_k,
                       int index_ic,
                       int index_tp,
                       int index_tau,
                       double ** sources,
                       double * source
                       ) {

  //double k,k_max,k_previous;
  //  double source_max,source_previous;
  //double scaled_factor,log_scaled_factor;

  /** - use precomputed values */
  if (index_k < pma->k_size) {
    *source = sources[index_ic * ppt->tp_size[pma->index_md_vectors] + index_tp][index_tau * pma->k_size + index_k];
  }
  /** - extrapolate **/
  else {

    //k = pma->k[index_k];

    /**
     * --> Get last source and k, which are used in (almost) all methods
     */
    //k_max = pma->k[pma->k_size-1];
    //source_max = sources[index_ic * ppt->tp_size[pma->index_md_vectors] + index_tp][index_tau * pma->k_size + pma->k_size - 1];

    /**
     * --> Get previous source and k, which are used in best methods
     */
    //k_previous = pma->k[pma->k_size-2];
    //source_previous = sources[index_ic * ppt->tp_size[pma->index_md_vectors] + index_tp][index_tau * pma->k_size + pma->k_size - 2];

    *source=0.0;

  }
 
  return _SUCCESS_;
}



/**
 * This routine computes the variance of magnetic fluctuations in a
 * sphere of radius R at redshift z, sigma(R,z), or other similar derived
 * quantitites.
 *
 * The integral is performed until the maximum value of k_max defined
 * in the perturbation module. Here there is not automatic checking
 * that k_max is large enough for the result to be well
 * converged. E.g. to get an accurate sigma8 at R = 8 Mpc/h, the user
 * should pass at least about P_k_max_h/Mpc = 1.
 *
 * @param ppr          Input: pointer to precision structure
 * @param pba          Input: pointer to background structure
 * @param pma          Input: pointer to magnetic structure
 * @param R            Input: radius in Mpc
 * @param z            Input: redshift
 * @param sigma_output Input: quantity to be computed (sigma, sigma', ...)
 * @param result       Output: result
 * @return the error status
 */

int magnetic_sigmas_at_z(
                        struct precision * ppr,
                        struct background * pba,
                        struct magnetic * pma,
                        double R,
                        double z,
			//			enum out_sigmas sigma_output,
                        double * result
                        ) {

  double * out_pk;
  double * ddout_pk;

  /** - allocate temporary array for P(k,z) as a function of k */

  class_alloc(out_pk, pma->k_size*sizeof(double), pma->error_message);
  class_alloc(ddout_pk, pma->k_size*sizeof(double), pma->error_message);

  /** - get P(k,z) as a function of k, for the right z */

  class_call(magnetic_Bk_at_z(pba,
                             pma,
                             logarithmic,
                             z,
			     out_pk,
                             NULL),
             pma->error_message,
             pma->error_message);

  /** - spline it along k */

  class_call(array_spline_table_columns(pma->ln_k,
                                        pma->k_size,
                                        out_pk,
                                        1,
                                        ddout_pk,
                                        _SPLINE_EST_DERIV_,
                                        pma->error_message),
             pma->error_message,
             pma->error_message);

  /** - calll the function computing the sigmas */

  class_call(magnetic_sigmas(pma,
                            R,
                            out_pk,
                            ddout_pk,
                            pma->k_size,
                            ppr->sigma_k_per_decade,
			     //			    sigma_output,
			    result),
             pma->error_message,
             pma->error_message);

  /** - free allocated arrays */

  free(out_pk);
  free(ddout_pk);

  return _SUCCESS_;
}


int magnetic_Bk_at_z(
                    struct background * pba,
                    struct magnetic *pma,
                    enum linear_or_logarithmic mode,
		    double z,
                    double * out_pk, // array out_pk[index_k]
                    double * out_pk_ic // array out_pk_ic[index_k * pfo->ic_ic_size + index_ic1_ic2]
                    ) {
  double tau;
  double ln_tau;
  int index_k;
  int index_ic1;
  int index_ic2;
  int index_ic1_ic1;
  int index_ic2_ic2;
  int index_ic1_ic2;
  int last_index;
  short do_ic = _FALSE_;

  /** - check whether we need the decomposition into contributions from each initial condition */

  if ((pma->ic_size > 1) && (out_pk_ic != NULL))
    do_ic = _TRUE_;

  /** - case z=0 requiring no interpolation in z */
  if (z == 0) {

    for (index_k=0; index_k<pma->k_size; index_k++) {
      
      out_pk[index_k] = pma->ln_Bk_l[(pma->ln_tau_size-1)*pma->k_size+index_k];
      
      if (do_ic == _TRUE_) {
	for (index_ic1_ic2 = 0; index_ic1_ic2 < pma->ic_ic_size; index_ic1_ic2++) {
	  out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic2] =
	    pma->ln_Bk_ic_l[((pma->ln_tau_size-1)*pma->k_size+index_k)*pma->ic_ic_size+index_ic1_ic2];
	}
      }
    }
  }

  /** - interpolation in z */
  else {

    class_test(pma->ln_tau_size == 1,
               pma->error_message,
               "You are asking for the magnetic power spectrum at z=%e but the code was asked to store it only at z=0. You probably forgot to pass the input parameter z_max_pk (see explanatory.ini)",z);

    /** --> get value of contormal time tau */
    class_call(background_tau_of_z(pba,
                                   z,
                                   &tau),
               pba->error_message,
               pma->error_message);

    ln_tau = log(tau);
    last_index = pma->ln_tau_size-1;

    /** -> check that tau is in pre-computed table */

    if (ln_tau <= pma->ln_tau[0]) {

      /** --> if ln(tau) much too small, raise an error */
      class_test(ln_tau<pma->ln_tau[0]-100.*_EPSILON_,
                 pma->error_message,
                 "requested z was not inside of tau tabulation range (Requested ln(tau_=%.10e, Min %.10e). Solution might be to increase input parameter z_max_pk (see explanatory.ini)",ln_tau,pma->ln_tau[0]);

      /** --> if ln(tau) too small but within tolerance, round it and get right values without interpolating */
      ln_tau = pma->ln_tau[0];

      for (index_k = 0 ; index_k < pma->k_size; index_k++) {
	out_pk[index_k] = pma->ln_Bk_l[index_k];
	if (do_ic == _TRUE_) {
	  for (index_ic1_ic2 = 0; index_ic1_ic2 < pma->ic_ic_size; index_ic1_ic2++) {
	    out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic2] = pma->ln_Bk_ic_l[index_k * pma->ic_ic_size + index_ic1_ic2];
	  }
	}
      }
    }

    else if (ln_tau >= pma->ln_tau[pma->ln_tau_size-1]) {

      /** --> if ln(tau) much too large, raise an error */
      class_test(ln_tau>pma->ln_tau[pma->ln_tau_size-1]+_EPSILON_,
                 pma->error_message,
                 "requested z was not inside of tau tabulation range (Requested ln(tau_=%.10e, Max %.10e) ",ln_tau,pma->ln_tau[pma->ln_tau_size-1]);

      /** --> if ln(tau) too large but within tolerance, round it and get right values without interpolating */
      ln_tau = pma->ln_tau[pma->ln_tau_size-1];

      for (index_k = 0 ; index_k < pma->k_size; index_k++) {
	out_pk[index_k] = pma->ln_Bk_l[(pma->ln_tau_size-1) * pma->k_size + index_k];
	if (do_ic == _TRUE_) {
	  for (index_ic1_ic2 = 0; index_ic1_ic2 < pma->ic_ic_size; index_ic1_ic2++) {
	    out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic2] = pma->ln_Bk_ic_l[((pma->ln_tau_size-1) * pma->k_size + index_k) * pma->ic_ic_size + index_ic1_ic2];
	  }
	}
      }
    }

    /** -> tau is in pre-computed table: interpolate */
    else {

      /** --> interpolate P_l(k) at tau from pre-computed array */
      class_call(array_interpolate_spline(pma->ln_tau,
					  pma->ln_tau_size,
					  pma->ln_Bk_l,
					  pma->ddln_Bk_l,
					  pma->k_size,
					  ln_tau,
					  &last_index,
					  out_pk,
					  pma->k_size,
					  pma->error_message),
                   pma->error_message,
                   pma->error_message);

        /** --> interpolate P_ic_l(k) at tau from pre-computed array */
      if (do_ic == _TRUE_) {
	class_call(array_interpolate_spline(pma->ln_tau,
					    pma->ln_tau_size,
					    pma->ln_Bk_ic_l,
					    pma->ddln_Bk_ic_l,
					    pma->k_size*pma->ic_ic_size,
					    ln_tau,
					    &last_index,
					    out_pk_ic,
					    pma->k_size*pma->ic_ic_size,
					    pma->error_message),
		   pma->error_message,
		   pma->error_message);
      }
    }
  }
  
  /** - so far, all output stored in logarithmic format. Eventually, convert to linear one. */

  if (mode == linear) {

    /** --> loop over k */
    for (index_k=0; index_k<pma->k_size; index_k++) {

      /** --> convert total spectrum */
      out_pk[index_k] = exp(out_pk[index_k]);

      if (do_ic == _TRUE_) {
        /** --> convert contribution of each ic (diagonal elements) */
        for (index_ic1=0; index_ic1 < pma->ic_size; index_ic1++) {
          index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pma->ic_size);

          out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic1] = exp(out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic1]);
        }

        /** --> convert contribution of each ic (non-diagonal elements) */
        for (index_ic1=0; index_ic1 < pma->ic_size; index_ic1++) {
          for (index_ic2=index_ic1+1; index_ic2 < pma->ic_size; index_ic2++) {
            index_ic1_ic1 = index_symmetric_matrix(index_ic1,index_ic1,pma->ic_size);
            index_ic2_ic2 = index_symmetric_matrix(index_ic2,index_ic2,pma->ic_size);
            index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,pma->ic_size);

            /* P_ic1xic2 = cos(angle) * sqrt(P_ic1 * P_ic2) */
            out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic2]
              = out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic2]
              *sqrt(out_pk_ic[index_k * pma->ic_ic_size + index_ic1_ic1]
                    *out_pk_ic[index_k * pma->ic_ic_size + index_ic2_ic2]);
          }
        }
      }
    }
  }
  

  return _SUCCESS_;
}



/**
 *
 * @param pma          Input: pointer to magnetic structure
 * @param R            Input: scale at which to compute sigma
 * @param lnpk_l       Input: array of ln(P(k))
 * @param ddlnpk_l     Input: its spline along k
 * @param k_size       Input: dimension of array lnpk_l, normally pma->k_size, but inside hmcode it its increased by extrapolation to pma->k_extra_size
 * @param k_per_decade Input: logarithmic step for the integral (recommended: pass ppr->sigma_k_per_decade)
 * @param sigma_output Input: quantity to be computed (sigma, sigma', ...)
 * @param result       Output: result
 * @return the error status
 */

int magnetic_sigmas(
                   struct magnetic * pma,
                   double R,
                   double * lnpk_l,
                   double * ddlnpk_l,
                   int k_size,
                   double k_per_decade,
		   //                   enum out_sigmas sigma_output,
                   double * result
                   ) {
  double pk, lnpk;

  double * array_for_sigma;
  int index_num;
  int index_x;
  int index_y;
  int index_ddy;
  int i=0;
  int integrand_size;
  int last_index=0;

  double k,W,W_prime,x,t;

  /** - allocate temporary array for an integral over y(x) */

  class_define_index(index_x,  _TRUE_,i,1); // index for x
  class_define_index(index_y,  _TRUE_,i,1); // index for integrand
  class_define_index(index_ddy,_TRUE_,i,1); // index for its second derivative (spline method)
  index_num=i;                              // number of columns in the array

  integrand_size=(int)(log(pma->k[k_size-1]/pma->k[0])/log(10.)*k_per_decade)+1;
  class_alloc(array_for_sigma,
              integrand_size*index_num*sizeof(double),
              pma->error_message);

  /** - fill the array with values of k and of the integrand */

  for (i=0; i<integrand_size; i++) {

    k=pma->k[0]*pow(10.,i/k_per_decade);

    if (i==0) {
      pk = exp(lnpk_l[0]);
    }
    else {
      class_call(array_interpolate_spline(
                                          pma->ln_k,
                                          k_size,
                                          lnpk_l,
                                          ddlnpk_l,
                                          1,
                                          log(k),
                                          &last_index,
                                          &lnpk,
                                          1,
                                          pma->error_message),
                 pma->error_message,
                 pma->error_message);

      pk = exp(lnpk);
    }

    t = 1./(1.+k);
    if (i == (integrand_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;

    //Old window function (in real space hence this is the Fourier transform kernel)
    /*if (x<0.01)
      W = 1.-x*x/10.;
    else
    W = 3./x/x/x*(sin(x)-x*cos(x));*/
    //For magnetic field we use rather a Gaussian smoothing
    W = exp(-1.*x*x/2.);
    array_for_sigma[(integrand_size-1-i)*index_num+index_x] = t;
    array_for_sigma[(integrand_size-1-i)*index_num+index_y] = k*k*k*pk*W*W/(t*(1.-t));
    
  }

  /** - spline the integrand */

  class_call(array_spline(array_for_sigma,
                          index_num,
                          integrand_size,
                          index_x,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pma->error_message),
             pma->error_message,
             pma->error_message);

  /** - integrate */

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                                  index_num,
                                                  integrand_size,
                                                  0, //integrand_size-1,
                                                  index_x,
                                                  index_y,
                                                  index_ddy,
                                                  result,
                                                  pma->error_message),
             pma->error_message,
             pma->error_message);

  /** - properly normalize the final result */

  *result = sqrt(*result/(2.*_PI_*_PI_));
  
  /** - free allocated array */

  free(array_for_sigma);

  return _SUCCESS_;
}
