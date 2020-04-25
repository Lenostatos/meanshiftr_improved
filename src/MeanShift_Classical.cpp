#include "LittleFunctionsCollection.h"

#include <Rcpp.h>
#include <cmath>


//' Mean shift clustering
//'
//' Adaptive mean shift clustering to delineate tree crowns from lidar point
//' clouds.
//'
//' @param pc Point cloud data in a three-column matrix where the columns
//'   represent X, Y and Z coordinates and each row represents one point.
//' @param CD2TH_fac Numeric scalar. Ratio of crown diameter to tree height.
//'   Determines kernel diameter based on the height of its center.
//' @param CH2TH_fac Numeric scalar. Ratio of crown height to tree height.
//'   Determines kernel height based on the height of its center.
//' @param MaxIter Integer scalar. Maximum number of iterations, i.e. steps that
//'   the kernel can move for each point. If no mode is found after
//'   \code{MaxIter} iterations, the centroid that was calculated last is
//'   treated as the mode.
//' @param UniformKernel Boolean scalar. Set to \code{TRUE} in order to turn off
//'   distance weighting within the kernel (defaults to \code{FALSE}).
//'
//' @return A data.frame with the coordinates in \code{pc} and three additional
//'   columns with the coordinates of the calculated modes.
//'
//' @export
// [[Rcpp::export]]
DataFrame MeanShift_Classical(
    NumericMatrix pc,
    double CD2TH_fac, double CH2TH_fac,
    int MaxIter = 200, bool UniformKernel = false
){
  // Create three vectors that all have the length of the incoming point cloud.
  // These vectors will store the coordinates of the calculated modes.
  int nrows{ pc.nrow() };
  NumericVector modesX(nrows);
  NumericVector modesY(nrows);
  NumericVector modesZ(nrows);

  // Process one point after the other.
  for(int i{ 0 }; i < nrows; i++){

    // Initialize variables to store the mean coordinates of all neighbors with
    // the actual coordinates of the current point from where the kernel starts
    // moving.
    double centroidX{ pc(i, 0) };
    double centroidY{ pc(i, 1) };
    double centroidZ{ pc(i, 2) };

    // Declare variables for storing the centroid of the previous iteration.
    double oldX;
    double oldY;
    double oldZ;

    // Initialize helper variables for the calculation of kernel centroids.
    double sumx{ 0 };
    double sumy{ 0 };
    double sumz{ 0 };
    double sump{ 0 };

    double verticalweight{ 0 };
    double horizontalweight{ 0 };
    double weight{ 0 };

    // Keep iterating while neither the mode nor the maximum number of
    // iterations are reached.
    int IterCounter{ 0 };
    do {

      sumx = 0;
      sumy = 0;
      sumz = 0;
      sump = 0;

      // Increase the iteration counter
      IterCounter += 1;

      // Calculate cylinder dimensions based on point height
      double radius = CD2TH_fac * centroidZ * 0.5;
      double height = CH2TH_fac * centroidZ;

      // Remember the centroid of the previous iteration.
       oldX = centroidX;
       oldY = centroidY;
       oldZ = centroidZ;

      // Loop through all points to identify the neighbors of the current point
      for(int j{ 0 }; j < nrows; j++) {

        double jx{ pc(i, 0) };
        double jy{ pc(i, 1) };
        double jz{ pc(i, 2) };

        if (InCylinder(
          jx, jy, jz,
          radius, height,
          centroidX, centroidY, centroidZ
        )) {
          // If the option of a uniform kernel is set to true calculate the mode
          // by summing up all coodinates and dividing by the number of points
          if (UniformKernel) {
            sumx += jx;
            sumy += jy;
            sumz += jz;
            sump += 1;
          }
          // Else calculate the mode by multiplying all coodinates by their
          // weights, depending on their relative position within the cylinder,
          // summing up the products and dividing by the sum of all weights
          else {
            verticalweight = EpanechnikovFunction(height, centroidZ, jz);
            horizontalweight = GaussFunction(
              radius, centroidX, centroidY, jx, jy
            );
            weight = verticalweight * horizontalweight;

            sumx += weight * jx;
            sumy += weight * jy;
            sumz += weight * jz;
            sump += weight;
          }
        }
      }

      centroidX = sumx / sump;
      centroidY = sumy / sump;
      centroidZ = sumz / sump;

    // If the new position is very close to the previous position (kernel
    // stopped moving), or if the maximum number of iterations is reached, stop
    // the iterations.
    } while (
      std::abs(centroidX - oldX) < 0.01
      && std::abs(centroidY - oldY) < 0.01
      && std::abs(centroidZ - oldZ) < 0.01
        && IterCounter < MaxIter
    );

    // Store the found position as the mode position for the current point
    modesX[i] = centroidX;
    modesY[i] = centroidY;
    modesZ[i] = centroidZ;
  }

  // Return the result as a data.frame with XYZ-coordinates of all points and
  // their corresponding modes
  return Rcpp::DataFrame::create(
    Rcpp::Named("X") = pc(_, 0),
    Rcpp::Named("Y") = pc(_, 1),
    Rcpp::Named("Z") = pc(_, 2),
    Rcpp::Named("modeX") = modesX,
    Rcpp::Named("modeY") = modesY,
    Rcpp::Named("modeZ") = modesZ
  );
}
