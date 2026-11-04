 #' This file provides helper functionsto quickly generate 
#' synthetic spatial data, meshes, and latent fields to test the MCMC samplers.

#' Generate Spatial Toy Dataset with Optional Covariance Structure
#'
#' Generates synthetic spatial data on a planar domain using the hierarchical model:
#' $Y_i \sim \text{Normal}(\mu_i, \sigma_e^2)$, where the latent field is 
#' $\mu = \alpha + \beta$, and $\beta$ is a spatial random effect drawn 
#' either from an independent normal or a user-defined covariance matrix.
#'
#' @param mesh An `fm_mesh_2d` object representing the spatial domain.
#' @param alpha Numeric scalar. The intercept / baseline effect. Default is 2.0.
#' @param sigma_e Numeric scalar. Standard deviation of the observation noise $\epsilon$. Default is 0.1.
#' @param Sigma Optional numeric matrix. Covariance matrix for the spatial random effect $\beta$. 
#'   If `NULL`, an identity-based independent normal is used.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing the synthetic observations and the ground truth:
#' \itemize{
#'   \item `y`: Numeric vector of observations.
#'   \item `mu_true`: Numeric vector of the *true* latent field ($\alpha + \beta$).
#'   \item `beta_true`: Numeric vector of the *true* pure spatial random effect.
#'   \item `alpha_true`: Numeric scalar of the *true* intercept.
#' }
#' @export
generate_toy_data <- function(mesh, 
                              alpha = 2.0, 
                              sigma_e = 0.1, 
                              Sigma = NULL, 
                              seed = NULL) {
  
  # TODO: Implement implementation logic
  
}