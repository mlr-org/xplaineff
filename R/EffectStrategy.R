#' EffectStrategy: Abstract base for effect-based tree strategies.
#'
#' Defines the common interface and shared fit logic for AleStrategy and PdStrategy.
#' Effect strategies implement ALE (Accumulated Local Effects) or PD (Partial Dependence)
#' based tree growing. Not exported.
#'
#' @field name (`character(1)`) \cr
#'   Strategy name (e.g. \code{"ale"}, \code{"pd"}).
#' @field tree_ref (`GadgetTree` or `NULL`) \cr
#'   Reference to the fitted tree; set after \code{$fit()}.
#' @field fit_timing (`numeric()` or `NULL`) \cr
#'   Fit timing (seconds) for global/regional fits.
#'
#' @keywords internal
EffectStrategy = R6::R6Class(
  "EffectStrategy",
  public = list(
    name = NULL,
    tree_ref = NULL,
    fit_timing = NULL,

    #' @description
    #' Create an EffectStrategy instance.
    #' @param name (`character(1)`) \cr
    #'   Strategy identifier (e.g. \code{"ale"}, \code{"pd"}).
    initialize = function(name) {
      self$name = name
    },

    #' @description
    #' Optional post-fit cleanup to release retained objects.
    #' Default is no-op; subclasses override when needed (see \code{AleStrategy$clean()}).
    clean = function() {
      invisible(NULL)
    }
  ),
  private = list(
    fit_tree_internal = function(tree, Z, Y, grid, objective_value_root_j, objective_value_root,
      verbose, vecb_remaining_features = NULL, early_stopping_goal = NULL) {
      # Create new tree. The root objective is the total risk over all features; when selective
      # early stopping is on, objective_value_remaining additionally holds the risk over just the
      # features still flagged as interacting at the root.
      root_objective_value_remaining = if (is.null(vecb_remaining_features)) NA_real_ else
        sum(objective_value_root_j[vecb_remaining_features], na.rm = TRUE)
      tree$root = Node$new(
        id = 1, depth = 1, subset_idx = seq_len(nrow(Z)), grid = grid,
        objective_value_parent = NA, int_imp_j = NULL,
        objective_value = objective_value_root,
        objective_value_remaining = root_objective_value_remaining,
        objective_value_j = objective_value_root_j, improvement_met = FALSE, int_imp = NULL,
        vecb_remaining_features = vecb_remaining_features, strategy = self
      )
      t_regional = system.time({
        tree$root$split_node(
          Z = Z, Y = Y,
          objective_value_root_j = objective_value_root_j,
          objective_value_root = objective_value_root,
          min_node_size = tree$min_node_size,
          n_quantiles = tree$n_quantiles,
          impr_par = tree$impr_par,
          depth = 1,
          max_depth = tree$n_split + 1,
          early_stopping_goal = early_stopping_goal,
          verbose = verbose
        )
      })[["elapsed"]]
      t_regional
    }
  )
)
