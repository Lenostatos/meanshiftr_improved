#include "littleFunctionsImproved.h"

#include <Rcpp.h>
#include <cmath>


//' Mean shift clustering
//'
//' Adaptive mean shift clustering to delineate tree crowns from lidar point
//' clouds.
//'
//' @param pointCloud Point cloud data in a three-column matrix where the
//'   columns represent X, Y and Z coordinates and each row represents one
//'   point.
//' @param crownDiameter2TreeHeight Numeric scalar. Ratio of crown diameter
//'   to tree height. Determines kernel diameter based on the height of its
//'   center.
//' @param crownHeight2TreeHeight Numeric scalar. Ratio of crown height to tree
//'   height. Determines kernel height based on the height of its center.
//' @param maxNumCentroidsPerMode Integer scalar. Maximum number of
//'   iterations, i.e. steps that the kernel can move for each point. If no mode
//'   is found after \code{maxNumCentroidsPerMode} iterations, the centroid
//'   that was calculated last is treated as the mode.
//'
//' @return A data.frame with the coordinates in \code{pointCloud} and three
//'   additional columns with the coordinates of the calculated modes.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame meanShiftClassicImproved(
    Rcpp::NumericMatrix pointCloud,
    double crownDiameter2TreeHeight, double crownHeight2TreeHeight,
    int maxNumCentroidsPerMode = 200
){
  // Create three vectors that all have the length of the incoming point cloud.
  // These vectors will store the coordinates of the calculated modes.
  int nrows{ pointCloud.nrow() };

  Rcpp::NumericVector modesX(nrows);
  Rcpp::NumericVector modesY(nrows);
  Rcpp::NumericVector modesZ(nrows);

  // Process one point after the other.
  for(int i{ 0 }; i < nrows; i++){

    // Get the current point's coordinates
    double centroidX{ pointCloud(i, 0) };
    double centroidY{ pointCloud(i, 1) };
    double centroidZ{ pointCloud(i, 2) };

    // Declare variables for storing the centroid of the previous iteration
    double oldX;
    double oldY;
    double oldZ;

    // Keep iterating while neither the mode nor the maximum number of
    // iterations are reached
    int numIterations{ 0 };

    do {
      // Increase the iteration counter
      numIterations += 1;

      // Initialize helper variables for the calculation of kernel centroids
      double sumX{ 0 };
      double sumY{ 0 };
      double sumZ{ 0 };
      double sumWeights{ 0 };

      // Remember the centroid of the previous iteration.
      oldX = centroidX;
      oldY = centroidY;
      oldZ = centroidZ;

      // Calculate cylinder dimensions based on point height
      double cylinderRadius = crownDiameter2TreeHeight * centroidZ * 0.5;
      double cylinderHeight = crownHeight2TreeHeight * centroidZ * 0.75;
      double cylinderMiddleZ = centroidZ + cylinderHeight * 1.0/6.0;

      // Loop through all points to identify the neighbors of the current point
      for(int j{ 0 }; j < nrows; j++) {

        double neighborX{ pointCloud(j, 0) };
        double neighborY{ pointCloud(j, 1) };
        double neighborZ{ pointCloud(j, 2) };

        if (
          intersectsCylinder(
            neighborX, neighborY, neighborZ,
            cylinderRadius, cylinderHeight,
            centroidX, centroidY, cylinderMiddleZ
          )
        ) {
          // Calculate the centroid by multiplying all coodinates by their
          // weights, depending on their relative position within the cylinder,
          // summing up the products and dividing by the sum of all weights
          double verticalweight{ calculateVerticalWeight(
            neighborZ, cylinderMiddleZ, cylinderHeight
          ) };
          double horizontalweight{ calculateHorizontalWeight(
              neighborX, neighborY, cylinderRadius, centroidX, centroidY
          ) };
          double weight{ verticalweight * horizontalweight };

          sumX += weight * neighborX;
          sumY += weight * neighborY;
          sumZ += weight * neighborZ;
          sumWeights += weight;
        }
      }

      centroidX = sumX / sumWeights;
      centroidY = sumY / sumWeights;
      centroidZ = sumZ / sumWeights;

      // If the new position is very close to the previous position (kernel
      // stopped moving), or if the maximum number of iterations is reached, stop
      // the iterations.
    } while (
      std::sqrt(
          std::pow(centroidX - oldX, 2.0)
        + std::pow(centroidY - oldY, 2.0)
        + std::pow(centroidZ - oldZ, 2.0)
      ) > 0.01
      && numIterations < maxNumCentroidsPerMode
    );

    // Store the found position as the mode position for the current point
    modesX[i] = centroidX;
    modesY[i] = centroidY;
    modesZ[i] = centroidZ;
  }

  // Return the result as a data.frame with XYZ-coordinates of all points and
  // their corresponding modes
  return Rcpp::DataFrame::create(
    Rcpp::Named("X") = pointCloud(Rcpp::_, 0),
    Rcpp::Named("Y") = pointCloud(Rcpp::_, 1),
    Rcpp::Named("Z") = pointCloud(Rcpp::_, 2),
    Rcpp::Named("modeX") = modesX,
    Rcpp::Named("modeY") = modesY,
    Rcpp::Named("modeZ") = modesZ
  );
}
