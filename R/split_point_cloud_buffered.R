#' Split point cloud into subsets with buffer areas around them
#'
#' The function splits one large point cloud into several smaller point clouds.
#' It allows to specify a buffer width around the core area where points are
#' included.
#'
#' @param point_cloud A data.table containing columns with x-, y-, and z-
#'   coordinates.
#' @param core_width Width of the core area in meters.
#' @param buffer_width Width of the buffer around the core area in meters.
#'
#' @return List of data.tables that each contains the coordinates of a point
#'   cloud subset together with a boolean column "Buffer" that labels core- and
#'   buffer-points.
#'
#' @importFrom data.table :=
#'
#' @export
split_point_cloud_buffered <- function(point_cloud, core_width, buffer_width) {

  # Convert to data.table
  point_cloud <- data.table::data.table(point_cloud)

  # Calculate the absolute lower left corner of the point cloud
  abs.llX <- plyr::round_any(
    min(point_cloud$X, na.rm = T), accuracy = core_width, f = floor
  )
  abs.llY <- plyr::round_any(
    min(point_cloud$Y, na.rm = T), accuracy = core_width, f = floor
  )

  # Calculate spatial indices for small subplots
  if (data.table::is.data.table(point_cloud)) {
    point_cloud[
      , sBPC_SpatID := calculate_plot_index(
        xcor = X, ycor = Y, res = core_width, minx = abs.llX, miny = abs.llY
      )
    ]
  }

  # Calculate coordinates of the lower left plot corners
  point_cloud[
    , sBPC_llX := plyr::round_any(X, accuracy = core_width, f = floor)
  ]
  point_cloud[
    , sBPC_llY := plyr::round_any(Y, accuracy = core_width, f = floor)
  ]

  # Make a unique data.table with all coordinate corner combinations and spatial indices
  coord.dt <- unique(subset(
    point_cloud, select = c("sBPC_SpatID", "sBPC_llX", "sBPC_llY")
  ))
  data.table::setorderv(
    coord.dt, cols = c("sBPC_SpatID", "sBPC_llX", "sBPC_llY")
  )
  data.table::setnames(
    coord.dt,
    old = names(coord.dt),
    new = c("sBPC_nSpatID", "sBPC_nllX", "sBPC_nllY")
  )

  # Create a copy of the point cloud which will become the result and to which
  # the buffers will be added successively
  result.dt <- data.table::copy(point_cloud)
  result.dt[, Buffer := 0]

  # Add the lower buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX]
  buf.dt[, sBPC_nllY := sBPC_llY + core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, Y < sBPC_nllY & Y >= sBPC_nllY - buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the lower left buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX + core_width]
  buf.dt[, sBPC_nllY := sBPC_llY + core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer_width & Y < sBPC_nllY & Y >= sBPC_nllY - buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the left buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX + core_width]
  buf.dt[, sBPC_nllY := sBPC_llY]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper left buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX + core_width]
  buf.dt[, sBPC_nllY := sBPC_llY - core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer_width & Y > sBPC_nllY + core_width & Y <= sBPC_nllY + core_width + buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX]
  buf.dt[, sBPC_nllY := sBPC_llY - core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, Y > sBPC_nllY + core_width & Y <= sBPC_nllY + core_width + buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper right buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX - core_width]
  buf.dt[, sBPC_nllY := sBPC_llY - core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + core_width & X <= sBPC_nllX + core_width + buffer_width & Y > sBPC_nllY + core_width & Y <= sBPC_nllY + core_width + buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the right buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX - core_width]
  buf.dt[, sBPC_nllY := sBPC_llY]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + core_width & X <= sBPC_nllX + core_width + buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the lower right buffer points
  buf.dt <- data.table::copy(point_cloud)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX - core_width]
  buf.dt[, sBPC_nllY := sBPC_llY + core_width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by = c("sBPC_nllX", "sBPC_nllY"), all.x = T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + core_width & X <= sBPC_nllX + core_width + buffer_width & Y < sBPC_nllY & Y >= sBPC_nllY - buffer_width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select = names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Remove all rows with NA as spatial index
  result.dt <- subset(result.dt, !is.na(sBPC_SpatID))

  # Split into a list of separate point clouds with buffers based on spatial index
  data.table::setorderv(result.dt, cols = c("sBPC_SpatID", "Buffer", "X", "Y"))
  result.list <- split(result.dt, by = "sBPC_SpatID")

  return(result.list)
}
