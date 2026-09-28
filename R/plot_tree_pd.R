#' Plot PD/ICE tree by depth and node.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects.
#' @param effect (R6 or `list()`) \cr
#'   Effect object (e.g. FeatureEffect).
#' @param data (`data.frame()` or `data.table()`) \cr
#'   Data with features and target.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column name.
#' @param color_ice,color_pd (`character(1)`) \cr
#'   Colors for ICE and PDP curves.
#' @param show_plot,show_point,mean_center (`logical(1)`) \cr
#'   Plot options.
#' @param depth (`integer()` or `NULL`) \cr
#'   Depth levels to render.
#' @param node_id (`integer()` or `NULL`) \cr
#'   Node IDs to render.
#' @param features (`character()` or `NULL`) \cr
#'   Feature subset.
#'
#' @return (`list()`) \cr
#'   Nested list (depth -> node -> patchwork).
#' @keywords internal
plot_tree_pd = function(tree, effect, data, target_feature_name,
  color_ice = "lightblue", color_pd = "lightcoral",
  show_plot = TRUE, show_point = FALSE, mean_center = TRUE,
  depth = NULL, node_id = NULL, features = NULL) {

  wide_mean_center = mean_center_ice(effect, feature_set = features, mean_center = mean_center)
  Y = wide_mean_center$Y
  grid_total = wide_mean_center$grid

  all_node_data = preprocess_node_data(tree, Y, grid_total, mean_center)

  depths_to_render = select_depths_to_render(tree, depth, node_id)
  if (length(depths_to_render) == 0) {
    cli::cli_warn("No valid depths to render")
    return(invisible(list()))
  }

  plot_list = list()

  for (depth_idx in depths_to_render) {
    prepared_data = all_node_data[[depth_idx]]
    if (is.null(node_id)) {
      nodes_to_render = seq_along(tree[[depth_idx]])
    } else {
      node_ids = suppressWarnings(as.integer(node_id))
      nodes_to_render = which(mlr3misc::map_lgl(tree[[depth_idx]], function(n) !is.null(n) && n$id %in% node_ids))
    }
    plots_at_depth = create_plots_for_depth(tree, prepared_data, data, target_feature_name,
      depth_idx, nodes_to_render, color_ice, color_pd,
      show_plot, show_point, mean_center)
    plot_list[[paste0("Depth_", depth_idx)]] = plots_at_depth
  }

  invisible(plot_list)
}
