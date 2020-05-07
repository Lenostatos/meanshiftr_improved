#' Parallel mean shift clustering for individual tree crown delineation
#'
#' The function provides the frame work to apply the adaptive mean shift 3D
#' (AMS3D) algorithm on several sub point clouds of a large investigation area
#' in parallel. It requires a list of sub point clouds as input and returns one
#' large clustered point cloud as output. The input should have buffer zones
#' around the focal areas. The buffer width should correspond to at least the
#' maximal possible tree crown radius.
#'
#' @param point_clouds List of point clouds in data.table format containing
#'   columns X, Y and Z (produced by the \code{split_point_cloud_buffered}
#'   function).
#' @param used_fraction_of_cores Fraction of available cores to use for
#'   parallelization.
#' @param version of the AMS3D algorithm. Can be set to "classic" (slow but
#'   precise also with small trees) or "voxel" (fast but based on rounded
#'   coordinates of 1-m precision) or "classic improved" (like classic but
#'   faster).
#' @param crown_diameter_2_tree_height Factor for the ratio of height to crown
#'   width. Determines kernel diameter based on its height above ground.
#' @param crown_height_2_tree_height Factor for the ratio of height to crown
#'   length. Determines kernel height based on its height above ground.
#' @param max_num_centroids_per_mode Maximum number of iterations, i.e. steps
#'   that the kernel can move for each point. If centroid is not found after all
#'   iteration, the last position is assigned as centroid and the processing
#'   jumps to the next point
#' @param min_num_neighbors_per_core Integer Scalar. The minimum number of
#'   neighbors that a point needs to have in order to be considered as a core
#'   point by the DBSCAN clustering algorithm.
#' @param neighborhood_radius Numeric Scalar. The radius of the space around a
#'   point that is treated as the point's neighborhood.
#' @param buffer_width Width of the buffer around the core area in meters.
#' @param min_height Minimum height above ground for a point to be considered in
#'   the analysis. Has to be > 0.
#'
#' @return data.table of point cloud with points labelled with tree IDs
#'
#' @export
segment_tree_crowns_parallel <- function(point_clouds,
                                         used_fraction_of_cores = 0.5,
                                         version = "classic",
                                         crown_diameter_2_tree_height,
                                         crown_height_2_tree_height,
                                         max_num_centroids_per_mode = 200,
                                         min_num_neighbors_per_core,
                                         neighborhood_radius,
                                         buffer_width = 10,
                                         min_height = 2) {

  # Calculate the number of cores
  num_cores <- parallel::detectCores()

  # Initiate cluster
  my_cluster <- parallel::makeCluster(num_cores * used_fraction_of_cores)

  # Prepare the environment on each child worker
  # Pass parameters to each worker
  parallel::clusterExport(
    cl = my_cluster,
    varlist = c(
      "version",
      "crown_diameter_2_tree_height", "crown_height_2_tree_height",
      "max_num_centroids_per_mode", "buffer_width", "min_height"
    ),
    envir = environment()
  )

  # Wrapper function that runs mean shift and deals with buffers
  mean_shift_buffered <- function(buffered_point_cloud) {

    # Remove points below a minimum height (ground and near ground returns)
    buffered_point_cloud <-
      subset(buffered_point_cloud, Z >= min_height)

    # Get margins of the point cloud
    min_x <- floor(min(buffered_point_cloud$X))
    max_x <- ceiling(max(buffered_point_cloud$X))
    min_y <- floor(min(buffered_point_cloud$Y))
    max_y <- ceiling(max(buffered_point_cloud$Y))

    # Get margins of the core area
    core_min_x <- floor(min(buffered_point_cloud[Buffer == 0, X]))
    core_max_x <- ceiling(max(buffered_point_cloud[Buffer == 0, X]))
    core_min_y <- floor(min(buffered_point_cloud[Buffer == 0, Y]))
    core_max_y <- ceiling(max(buffered_point_cloud[Buffer == 0, Y]))

    # Convert to 3-column matrix
    point_cloud_matrix <- as.matrix(buffered_point_cloud)
    point_cloud_matrix <- point_cloud_matrix[, 1:3]

    # Run the requested version of the mean shift algorithm
    if (version == "classic") {
      modes <- meanShiftClassic(
        pointCloud = point_cloud_matrix,
        crownDiameter2TreeHeight = crown_diameter_2_tree_height,
        crownHeight2TreeHeight = crown_height_2_tree_height,
        maxNumCentroidsPerMode = max_num_centroids_per_mode
      )
    } else if (version == "voxel") {
      range_x <- max_x - min_x
      range_y <- max_y - min_y
      modes <- MeanShift_Voxels(
        pc = point_cloud_matrix,
        H2CW_fac = crown_diameter_2_tree_height,
        H2CL_fac = crown_height_2_tree_height,
        UniformKernel = FALSE, MaxIter = max_num_centroids_per_mode,
        maxx = range_x, maxy = range_y, maxz = max_z
      )
    } else if (version == "improved") {
      modes <- meanShiftClassicImproved(
        pointCloud = point_cloud_matrix,
        crownDiameter2TreeHeight = crown_diameter_2_tree_height,
        crownHeight2TreeHeight = crown_height_2_tree_height,
        maxNumCentroidsPerMode = max_num_centroids_per_mode
      )
    }

    modes_data_table <-
      data.table::data.table(modes)

    # Identify mode clusters with the DBSCAN algorithm
    crown_ids <- dbscan::dbscan(
      modes_data_table[, .(modeX, modeY, modeZ)],
      eps = neighborhood_radius, minPts = min_num_neighbors_per_core + 1
    )$cluster

    modes_data_table <- data.table::data.table(
      modes_data_table, "crown_id" = crown_ids
    )

    # Extract all modes that were not part of a cluster
    unclustered_modes <- modes_data_table[crown_id == 0]
    modes_data_table <- modes_data_table[crown_id != 0]

    # Keep all unclustered modes that are within the core area
    unclustered_core_modes <- unclustered_modes[
        core_min_x <= modeX & modeX <= core_max_x
      & core_min_y <= modeY & modeY <= core_max_y
    ]

    # Get the mean position of every cluster
    cluster_means <- modes_data_table[
      , lapply(.SD, mean),
      by = crown_id,
      .SDcols = c("modeX", "modeY")
    ]
    data.table::setnames(
      cluster_means,
      old = c("modeX", "modeY"),
      new = c("mean_cluster_x", "mean_cluster_y")
    )
    modes_data_table <-
      merge(modes_data_table, cluster_means, by = "crown_id")

    # Only keep points whose cluster's mean position lies inside the core area
    core_cluster_data_table <- modes_data_table[
        core_min_x <= mean_cluster_x & mean_cluster_x <= core_max_x
      & core_min_y <= mean_cluster_y & mean_cluster_y <= core_max_y
    ]

    segmented_point_cloud <- data.table::rbindlist(
      list(
        core_cluster_data_table[, -c("mean_cluster_x", "mean_cluster_y")],
        unclustered_core_modes
      ),
      use.names = TRUE
    )

    # Collect the clustered point cloud in the results list
    return(segmented_point_cloud)
  }

  # Apply the mean shift wrapper function in parallel using pblapply to display
  # a progress bar
  res_list <- pbapply::pblapply(
    cl = my_cluster, X = point_clouds, FUN = mean_shift_buffered
  )

  parallel::stopCluster(my_cluster)

  # Treat unclustered modes separately
  unclustered_points <- data.table::rbindlist(lapply(
    res_list,
    FUN = function(segmented_point_cloud) {
      segmented_point_cloud[crown_id == 0]
    }
  ))

  # Remove unclustered points from the main data set
  for (i in seq_along(res_list)) {
    res_list[[i]] <- res_list[[i]][crown_id != 0]
  }

  # Ensure that crown IDs don't overlap between clusters from different threads
  crown_id_increment <- 0
  for (i in seq_along(res_list)) {
    res_list[[i]][, crown_id := crown_id + crown_id_increment]
    crown_id_increment <- 1 + max(res_list[[i]]$crown_id)
  }

  # Put all clustered points into one data table
  res_data_table <- data.table::rbindlist(res_list)

  # Add all unclustered points with crown ID 0
  res_data_table <- rbind(res_data_table, unclustered_points)

  return(res_data_table)
}


# Clusters for interactive testing
# set.seed(665544)
# n <- 1000
# modes_data_table <- data.table::data.table(
#   modeX = runif(10, 0, 10) + rnorm(n, sd = 0.2),
#   modeY = runif(10, 0, 10) + rnorm(n, sd = 0.2),
#   modeZ = runif(10, 0, 10) + rnorm(n, sd = 0.2)
# )
#
# res_list <-
#   list(data.table::copy(modes_data_table),
#        data.table::copy(modes_data_table),
#        data.table::copy(modes_data_table)
#   )
#
# hist(res_list[[3]]$crown_id)
