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
#' @param buffer_width Width of the buffer around the core area in meters.
#' @param min_height Minimum height above ground for a point to be considered in
#'   the analysis. Has to be > 0.
#' @param mode_rounding Numeric Scalar. Specifies the rounding accuracy for mode
#'   positions. After rounding all modes with the same coordinates are
#'   considered to belong to one tree crown.
#'
#' @return data.table of point cloud with points labelled with tree IDs
#'
#' @export
parallel_mean_shift <- function(point_clouds,
                                used_fraction_of_cores = 0.5,
                                version = "classic",
                                crown_diameter_2_tree_height,
                                crown_height_2_tree_height,
                                max_num_centroids_per_mode = 200,
                                buffer_width = 10,
                                min_height = 2, mode_rounding = 2) {

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
      "max_num_centroids_per_mode", "buffer_width",
      "min_height", "mode_rounding"
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
    range_x <- max_x - min_x
    range_y <- max_y - min_y

    # Get margins of the core area
    core_min_x <- floor(min(buffered_point_cloud[Buffer == 0, X]))
    core_max_x <- ceiling(max(buffered_point_cloud[Buffer == 0, X]))
    core_min_y <- floor(min(buffered_point_cloud[Buffer == 0, Y]))
    core_max_y <- ceiling(max(buffered_point_cloud[Buffer == 0, Y]))

    # Shift to coordinate origin
    buffered_point_cloud[, X := X - min_x]
    buffered_point_cloud[, Y := Y - min_y]

    # Convert to 3-column matrix
    point_cloud_matrix <- as.matrix(buffered_point_cloud)
    point_cloud_matrix <- point_cloud_matrix[, 1:3]

    # Run the requested version of the mean shift algorithm
    if (version == "classic") {
      clustered_point_cloud <- meanShiftClassic(
        pointCloud = point_cloud_matrix,
        crownDiameter2TreeHeight = crown_diameter_2_tree_height,
        crownHeight2TreeHeight = crown_height_2_tree_height,
        maxNumCentroidsPerMode = max_num_centroids_per_mode
      )
    } else if (version == "voxel") {
      clustered_point_cloud <- MeanShift_Voxels(
        pc = point_cloud_matrix,
        H2CW_fac = crown_diameter_2_tree_height,
        H2CL_fac = crown_height_2_tree_height,
        UniformKernel = FALSE, MaxIter = max_num_centroids_per_mode,
        maxx = range_x, maxy = range_y, maxz = my.maxz
      )
    } else if (version == "classic improved") {
      clustered_point_cloud <- meanShiftClassicImproved(
        pointCloud = point_cloud_matrix,
        crownDiameter2TreeHeight = crown_diameter_2_tree_height,
        crownHeight2TreeHeight = crown_height_2_tree_height,
        maxNumCentroidsPerMode = max_num_centroids_per_mode
      )
    }

    # Round the centroid coordinates
    clustered_point_cloud_data_table <-
      data.table::data.table(clustered_point_cloud)
    clustered_point_cloud_data_table[
      , rounded_mode_X := plyr::round_any(modeX, accuracy = mode_rounding)
    ]
    clustered_point_cloud_data_table[
      , rounded_mode_Y := plyr::round_any(modeY, accuracy = mode_rounding)
    ]
    clustered_point_cloud_data_table[
      , rounded_mode_Z := plyr::round_any(modeZ, accuracy = mode_rounding)
    ]

    # Shift back to original positions
    clustered_point_cloud_data_table[, X := X + min_x]
    clustered_point_cloud_data_table[, Y := Y + min_y]
    clustered_point_cloud_data_table[, modeX := modeX + min_x]
    clustered_point_cloud_data_table[, modeY := modeY + min_y]
    clustered_point_cloud_data_table[
      , rounded_mode_X := rounded_mode_X + min_x
    ]
    clustered_point_cloud_data_table[
      , rounded_mode_Y := rounded_mode_Y + min_y
    ]

    # Subset tree clusters with centers inside the core area of the
    # focal subplot and discard the clusters with centers in the buffer area
    clustered_point_cloud_data_table <- subset(
      clustered_point_cloud_data_table,
      rounded_mode_X >= core_min_x & rounded_mode_X < core_max_x
      & rounded_mode_Y >= core_min_y & rounded_mode_Y < core_max_y
    )

    # Collect the clustered point cloud in the results list
    return(clustered_point_cloud_data_table)
  }

  # Apply the mean shift wrapper function in parallel using pblapply to display
  # a progress bar
  res_list <- pbapply::pblapply(
    cl = my_cluster, X = point_clouds, FUN = mean_shift_buffered
  )

  # Bind all point clouds from the list in one large data.table
  res_data_table <- data.table::rbindlist(res_list)

  # Assign IDs to each cluster based on the rounded coordinates
  res_data_table[
    , ID := .GRP,
    by = .(rounded_mode_X, rounded_mode_Y, rounded_mode_Z)
  ]

  parallel::stopCluster(my_cluster)

  return(res_data_table)
}
