#' Convert Recursive Tree to Depth-based List
#'
#' Given root_node and optional max_depth: traverses children recursively, collects nodes by depth.
#' Returns list of lists (depth 1, 2, ...), each element a list of Node objects at that depth.
#'
#' @param root_node (`Node` or `NULL`) \cr
#'   Root node of the tree.
#' @param max_depth (`integer(1)` or `NULL`) \cr
#'   Maximum depth to traverse; \code{NULL} = no limit.
#'
#' @return (`list()`) \cr
#'   List where each element is a depth level containing Node objects.
#' @keywords internal
convert_tree_to_list = function(root_node, max_depth = NULL) {
  if (is.null(root_node)) {
    return(list())
  }
  traverse_and_collect = function(node, current_depth, depth_list) {
    if (is.null(node)) {
      return(depth_list)
    }
    if (!is.null(max_depth) && current_depth > max_depth) {
      return(depth_list)
    }
    while (length(depth_list) < current_depth) {
      depth_list[[length(depth_list) + 1]] = list()
    }
    depth_list[[current_depth]][[length(depth_list[[current_depth]]) + 1]] = node
    if (!is.null(node$children)) {
      if (!is.null(node$children$left_child)) {
        depth_list = traverse_and_collect(node$children$left_child,
          current_depth + 1, depth_list)
      }
      if (!is.null(node$children$right_child)) {
        depth_list = traverse_and_collect(node$children$right_child,
          current_depth + 1, depth_list)
      }
    }
    depth_list
  }
  tree_depths = traverse_and_collect(root_node, 1, list())
  tree_depths = Filter(function(depth) length(depth) > 0, tree_depths)
  tree_depths
}
