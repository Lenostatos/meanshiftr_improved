#include <Rcpp.h>
#include <cmath>
#include "LittleFunctionsCollection.h"
using namespace Rcpp;


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
DataFrame meanShiftClassic(
    NumericMatrix pointCloud,
    double crownDiameter2TreeHeight, double crownHeight2TreeHeight,
    int maxNumCentroidsPerMode = 200
){
  // Create three vectors that all have the length of the incoming point cloud.
  // These vectors will store the coordinates of the calculated modes.
  int nrows = pointCloud.nrow();

  NumericVector modesX(nrows);
  NumericVector modesY(nrows);
  NumericVector modesZ(nrows);

  // Process one point after the other.
  for(int i = 0; i < nrows; i++){

    // Initialize variables to store the mean coordinates of all neighbors with
    // the actual coordinates of the current point from where the kernel starts
    // moving.
    double centroidX = (double) pointCloud(i, 0);
    double centroidY = (double) pointCloud(i, 1);
    double centroidZ = (double) pointCloud(i, 2);

    // Declare variables for storing the centroid of the previous iteration.
    double oldX = -100.0;
    double oldY = -100.0;
    double oldZ = -100.0;

    // Initialize helper variables for the calculation of kernel centroids.
    double sumX = 0.0;
    double sumY = 0.0;
    double sumZ = 0.0;
    double sumP = 0.0;

    double verticalWeight = 0.0;
    double horizontalWeight = 0.0;
    double weight = 0.0;

    // Keep iterating while neither the mode nor the maximum number of
    // iterations are reached.
    int numIterations = 0;
    do {

      sumX = 0.0;
      sumY = 0.0;
      sumZ = 0.0;
      sumP = 0.0;

      // Increase the iteration counter
      numIterations = numIterations + 1;

      // Calculate cylinder dimensions based on point height
      double radius = crownDiameter2TreeHeight * centroidZ * 0.5;
      double height = crownHeight2TreeHeight * centroidZ;

      // Remember the centroid of the previous iteration.
       oldX = centroidX;
       oldY = centroidY;
       oldZ = centroidZ;

      // Loop through all points to identify the neighbors of the current point
      for(int j = 0; j < nrows; j++) {

        double pointX = (double) pointCloud(j, 0);
        double pointY = (double) pointCloud(j, 1);
        double pointZ = (double) pointCloud(j, 2);

        if (InCylinder(
          pointX, pointY, pointZ,
          radius, height,
          centroidX, centroidY, centroidZ
        ) == true) {
          // Calculate the mode by multiplying
          // all coodinates by their weights, depending on their relative
          // position within the cylinder, summing up the products and dividing
          // by the sum of all weights
          verticalWeight = EpanechnikovFunction(height, centroidZ, pointZ);
          horizontalWeight = GaussFunction(
            radius, centroidX, centroidY, pointX, pointY
          );
          weight = verticalWeight * horizontalWeight;

          sumX = sumX + weight * pointX;
          sumY = sumY + weight * pointY;
          sumZ = sumZ + weight * pointZ;
          sumP = sumP + weight;
        }
      }

      centroidX = sumX / sumP;
      centroidY = sumY / sumP;
      centroidZ = sumZ / sumP;

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
  return DataFrame::create(
    _("X") = pointCloud(_, 0),
    _("Y") = pointCloud(_, 1),
    _("Z") = pointCloud(_, 2),
    _("modeX") = modesX,
    _("modeY") = modesY,
    _("modeZ") = modesZ
  );
}
