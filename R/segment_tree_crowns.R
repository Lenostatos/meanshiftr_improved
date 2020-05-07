#' Calculate crown IDs for trees in a point cloud
#'
#' @param point_cloud A data.frame or data.table. Its first three columns are
#'   expected to hold coordinates.
#' @param version Character. One of "classic" or "improved".
#'
#' @export
segment_tree_crowns <- function(point_cloud,
                                version = "classic",
                                crown_diameter_2_tree_height,
                                crown_height_2_tree_height,
                                max_num_centroids_per_mode = 200,
                                min_num_neighbors_per_core,
                                neighborhood_radius) {

  assertthat::assert_that(version %in% c("classic", "improved"))

  if (version == "classic") {
    modes <- data.table::as.data.table(
      meanShiftClassic(as.matrix(point_cloud[, 1:3]),
                       crown_diameter_2_tree_height,
                       crown_height_2_tree_height,
                       max_num_centroids_per_mode))
  } else if (version == "improved") {
    modes <- data.table::as.data.table(
      meanShiftClassicImproved(as.matrix(point_cloud[, 1:3]),
                               crown_diameter_2_tree_height,
                               crown_height_2_tree_height,
                               max_num_centroids_per_mode))
  }

  crown_ids <- dbscan::dbscan(modes[, .(modeX, modeY, modeZ)],
                              eps = neighborhood_radius,
                              minPts = min_num_neighbors_per_core + 1)$cluster

  data.table::data.table(modes, crown_id = crown_ids)
}
