#' GadgetTree: Generalized Additive Decomposition of Global Effects Tree
#'
#' @description
#' Wrapper for effect-based trees: given a strategy (pd/ale), fits tree via
#' \code{$fit()}, plots effects via \code{$plot()}, extracts splits via
#' \code{$extract_split_info()}.
#' Delegates all effect logic to the strategy.
#'
#' @usage NULL
#' @format [R6::R6Class] object.
#'
#' @section Construction:
#' ```
#' t = GadgetTree$new(strategy, n_split = 2, impr_par = 0.1, min_node_size = 10, n_quantiles = NULL)
#' ```
#' * `strategy` :: ([AleStrategy] | [PdStrategy])\cr
#'   Effect strategy object.
#' * `n_split` :: `integer(1)`\cr
#'   Maximum split depth, counted as splits along any root-to-leaf path.
#' * `impr_par` :: `numeric(1)`\cr
#'   Improvement threshold.
#' * `min_node_size` :: `integer(1)`\cr
#'   Minimum samples per node.
#' * `n_quantiles` :: `integer(1)` or `NULL`\cr
#'   Number of quantile cut points for numeric split candidates; `NULL` = use all unique values.
#'
#' @field strategy (PdStrategy | AleStrategy) \cr
#'   Effect-specific logic.
#' @field root (`Node`) \cr
#'   Root node.
#' @field n_split (`integer(1)`) \cr
#'   Maximum split depth, counted as splits along any root-to-leaf path.
#' @field impr_par (`numeric(1)`) \cr
#'   Improvement threshold.
#' @field min_node_size (`integer(1)`) \cr
#'   Minimum samples per node.
#' @field n_quantiles (`integer(1)` or `NULL`) \cr
#'   Number of quantile cut points for numeric split candidates.
#' @field split_benchmark (`list()`) \cr
#'   Internal split timing records.
#' @field tree_list_cache (`list()` or `NULL`) \cr
#'   Cached depth-based tree list; invalidated on \code{$fit()}.
#'
#' @details
#' This class manages the overall tree structure and delegates effect-specific operations
#' (such as splitting and plotting) to the provided strategy object. It is the main entry
#' point for fitting and visualizing effect-based decision trees in the xplaineff package.
#'
#' @examples
#' \dontrun{
#' # PD: `data` and `target_feature_name` come first; pass `effect` and/or `model` via `...`.
#' pd_strat = PdStrategy$new()
#' tree = GadgetTree$new(strategy = pd_strat, n_split = 2)
#' tree$fit(data = dat, target_feature_name = "y", effect = effect_obj)
#' tree$plot(data = dat, target_feature_name = "y", effect = effect_obj)
#' tree$plot(data = dat, target_feature_name = "y")  # uses cached effect from fit()
#' tree$plot_tree_structure()
#' split_info = tree$extract_split_info()
#'
#' # ALE: pass `model` (and optional `n_intervals`, `ale_engine`, etc.).
#' tree_ale = GadgetTree$new(strategy = AleStrategy$new(), n_split = 2)
#' tree_ale$fit(data = dat, target_feature_name = "y", model = fitted_model)
#' tree_ale$plot(data = dat, target_feature_name = "y")
#' }
#'
#' @export
GadgetTree = R6::R6Class(
  "GadgetTree",
  public = list(
    strategy = NULL,
    root = NULL,
    n_split = NULL,
    impr_par = NULL,
    min_node_size = NULL,
    n_quantiles = NULL,
    split_benchmark = NULL,
    tree_list_cache = NULL,

    #' @description
    #' Initialize tree parameters.
    #' @param strategy (PdStrategy | AleStrategy) \cr
    #'   Strategy object.
    #' @param n_split (`integer(1)`) \cr
    #'   Maximum split depth, counted as splits along any root-to-leaf path.
    #' @param impr_par (`numeric(1)`) \cr
    #'   Improvement threshold.
    #' @param min_node_size (`integer(1)`) \cr
    #'   Minimum node size.
    #' @param n_quantiles (`integer(1)` or `NULL`) \cr
    #'   Number of quantile cut points for numeric split candidates.
    initialize = function(strategy, n_split = 2, impr_par = 0.1, min_node_size = 10, n_quantiles = NULL) {
      checkmate::assert_r6(strategy, .var.name = "strategy")
      checkmate::assert_integerish(n_split, len = 1, lower = 0, any.missing = FALSE, .var.name = "n_split")
      checkmate::assert_numeric(impr_par, lower = 0, len = 1, any.missing = FALSE, .var.name = "impr_par")
      checkmate::assert_integerish(min_node_size, len = 1, lower = 1, any.missing = FALSE, .var.name = "min_node_size")
      checkmate::assert_integerish(n_quantiles, len = 1, null.ok = TRUE, .var.name = "n_quantiles")
      self$strategy = strategy
      self$n_split = n_split
      self$impr_par = impr_par
      self$min_node_size = min_node_size
      self$n_quantiles = n_quantiles
      self$split_benchmark = list()
    },

    #' @description
    #' Fit tree via \code{strategy$fit()}.
    #' @param data (`data.frame()`) \cr
    #'   Data with features and target.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target name.
    #' @param feature_set (`character()` or `NULL`) \cr
    #'   Features for effect; \code{NULL} = all.
    #' @param split_feature (`character()` or `NULL`) \cr
    #'   Features for splitting; \code{NULL} = all.
    #' @param gadget_improvements (`character(1)` or `NULL`) \cr
    #'   Optional selective early stopping method, which dynamically drops non-interacting
    #'   features from the split-candidate set Z and the risk set S. One of
    #'   \code{"plain_risk"} (Method 1: absolute normalized risk),
    #'   \code{"risk_reduction"} (Method 2: relative risk reduction of the split; being
    #'   reduction-based it cannot drop features before the first split), or
    #'   \code{"interaction_fraction"} (Method 3: interaction share
    #'   \eqn{R_j / (R_j + B_j + \delta)} of the feature's own local-effect variance).
    #'   \code{NULL} disables it. PD strategy only.
    #' @param gadget_impr_args (`list()`) \cr
    #'   Named list of method-specific arguments for \code{gadget_improvements}:
    #'   \code{tau} (`numeric(1)`, default 0.05) for all methods, and \code{delta}
    #'   (`numeric(1)`, default 1e-12) regularizing the denominator of Method 3.
    #' @param ... (`list()`) \cr
    #'   Strategy-specific arguments passed to \code{$fit()}.
    #'   For [AleStrategy]: \code{model} or \code{effect}, plus optional
    #'   \code{n_intervals}, \code{predict_fun}, \code{order_method},
    #'   and \code{ale_engine}.
    #'   For [PdStrategy]: \code{effect}, or \code{model} with optional
    #'   \code{predict_fun}, \code{n_grid}, and \code{pd_engine}.
    #' @return (`GadgetTree`) \cr
    #'   The tree, invisibly.
    fit = function(data, target_feature_name, feature_set = NULL, split_feature = NULL,
      gadget_improvements = NULL, gadget_impr_args = NULL, verbose = 0, ...) {
      checkmate::assert_data_frame(data, .var.name = "data")
      checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
      checkmate::assert_subset(target_feature_name, colnames(data), .var.name = "target_feature_name")
      checkmate::assert_choice(gadget_improvements,
        c("plain_risk", "risk_reduction", "interaction_fraction"),
        null.ok = TRUE, .var.name = "gadget_improvements")
      checkmate::assert_list(gadget_impr_args, names = "unique", null.ok = TRUE, .var.name = "gadget_impr_args")
      if (!is.null(gadget_improvements) && !inherits(self$strategy, "PdStrategy")) {
        cli::cli_abort("{.arg gadget_improvements} is currently only implemented for {.cls PdStrategy}.")
      }
      self$split_benchmark = list()
      self$tree_list_cache = NULL
      if (verbose > 0) {
        print("Starting gadgetTree$fit; using print")
        flush.console()
        message("Starting gadgetTree$fit; using message")
        # browser()
      }

      # The strategy is responsible for validating and handling its own arguments (via ...)
      result = self$strategy$fit(
        tree = self,
        data = data,
        target_feature_name = target_feature_name,
        feature_set = feature_set,
        split_feature = split_feature,
        gadget_improvements = gadget_improvements,
        gadget_impr_args = gadget_impr_args,
        verbose = verbose,
        ...
      )

      self$strategy$clean()
      invisible(result)
    },

    #' @description
    #' Plot tree via \code{strategy$plot()}.
    #' @param data (`data.frame()`) \cr
    #'   Data.
    #' @param target_feature_name (`character(1)`) \cr
    #'   Target name.
    #' @param effect (R6 or `list()` or `NULL`) \cr
    #'   Optional effect object; omit or \code{NULL} uses strategy-cached effect from \code{$fit()}.
    #' @param depth (`integer()` or `NULL`) \cr
    #'   Depths to plot.
    #' @param node_id (`integer()` or `NULL`) \cr
    #'   Node IDs to plot.
    #' @param features (`character()` or `NULL`) \cr
    #'   Features to plot.
    #' @param ... Plot arguments.
    #' @return (`list()`) \cr
    #'   Nested list (depth -> node -> patchwork).
    plot = function(data, target_feature_name, effect = NULL, depth = NULL, node_id = NULL, features = NULL, ...) {
      # browser()
      tree_list = self$get_tree_list()
      self$strategy$plot(tree = tree_list, effect = effect, data = data,
        target_feature_name = target_feature_name,
        depth = depth, node_id = node_id, features = features, ...)
    },

    #' @description
    #' Converts root to depth-list and calls \code{plot_tree_structure()}. Prints graph.
    #' @param label_wrap_width (`integer(1)` or `NULL`) \cr
    #'   Wrap node labels to this many characters per line; \code{NULL} disables wrapping.
    #' @param node_spread_x,node_spread_y (`numeric(1)`) \cr
    #'   Layout stretch factors for the ggraph \code{"tree"} layout (larger values separate nodes).
    plot_tree_structure = function(label_wrap_width = 34L, node_spread_x = 1.55, node_spread_y = 1.12) {
      checkmate::assert_integerish(label_wrap_width, len = 1L, lower = 8L, null.ok = TRUE,
        upper = 100L, any.missing = FALSE, .var.name = "label_wrap_width")
      checkmate::assert_numeric(node_spread_x, len = 1L, lower = 0.5, finite = TRUE, .var.name = "node_spread_x")
      checkmate::assert_numeric(node_spread_y, len = 1L, lower = 0.5, finite = TRUE, .var.name = "node_spread_y")
      # browser()
      plot_tree_structure(self$get_tree_list(),
        label_wrap_width = label_wrap_width,
        node_spread_x = node_spread_x,
        node_spread_y = node_spread_y)
    },

    #' @description
    #' Extract split info from tree.
    #' @param include_timing (`logical(1)`) \cr
    #'   Whether to include internal split timings in the output.
    #' @return (`data.frame()`) \cr
    #'   Split info: depth, id, split_feature, split_value, int_imp, etc.
    extract_split_info = function(include_timing = FALSE) {
      checkmate::assert_flag(include_timing, .var.name = "include_timing")
      # browser()
      split_benchmark = if (isTRUE(include_timing)) self$split_benchmark else NULL
      extract_split_info(self$get_tree_list(), split_benchmark = split_benchmark)
    },

    #' @description
    #' Get depth-based tree list (cached). Invalidated on \code{$fit()}.
    #' @return (`list()`) \cr
    #'   Depth-based list of nodes.
    #' @keywords internal
    get_tree_list = function() {
      if (is.null(self$tree_list_cache)) {
        self$tree_list_cache = convert_tree_to_list(self$root, self$n_split + 1)
      }
      self$tree_list_cache
    }
  )
)
