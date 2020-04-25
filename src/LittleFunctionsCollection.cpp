#include "LittleFunctionsCollection.h"

#include <Rcpp.h>
#include <cmath>  // for std::pow, std::abs, std::min, std::hypot, std::exp


bool InCylinder(
  double x, double y, double z,
  double radius, double height,
  double centerX, double centerY, double centerZ
){
  return (std::pow((x - centerX), 2.0) + std::pow((y - centerY), 2.0)
            <= std::pow(radius, 2.0))
    && (z >= (centerZ - (0.5 * height)))
    && (z <= (centerZ + (0.5 * height)));
}


double VerticalDistance(double height, double centerZ, double pointZ) {
  // Calculate the relative distances to top and bottom...
  double BottomDistance{
    std::abs((centerZ - height/4.0 - pointZ) / (height * 3.0/8.0))
  };
  double TopDistance{
    std::abs((centerZ + height/2 - pointZ) / (height * 3.0/8.0))
  };

  // ...and return the smaller of the two.
  return std::min(BottomDistance, TopDistance);
}

//Equivalent R code
//distx <- function(h, CtrZ, PointZ){
//  bottomdist <- abs((CtrZ-h/4-PointZ)/(3*h/8))
//  topdist <- abs((CtrZ+h/2-PointZ)/(3*h/8))
//  mindist <- pmin(bottomdist, topdist)
//  return(mindist)
//}


short VerticalMask(double height, double centerZ, double pointZ) {
  if((pointZ >= centerZ - height/4.0) && (pointZ <= centerZ + height/2.0)) {
    return 1;
  } else {
    return 0;
  }
}

//Equivalent R code
//maskx <- function(h, CtrZ, PointZ){
//  maskvec <- ifelse(PointZ >= CtrZ-h/4 & PointZ <= CtrZ+h/2, 1, 0)
//  return(maskvec)
//}


// Epanechnikov function for vertical filter
double EpanechnikovFunction(double height, double centerZ, double pointZ) {
  return VerticalMask(height, centerZ, pointZ)
    * (1 - std::pow((1 - VerticalDistance(height, centerZ, pointZ)), 2.0));
}

//Equivalent R code
//Epanechnikov <- function(h, CtrZ, PointZ){
//  output <- maskx(h, CtrZ, PointZ)*(1-(1-distx(h, CtrZ, PointZ))^2)
//  return(output)
//}


// Gauss function for horizontal filter
double GaussFunction(
  double radius, double centerX, double centerY, double pointX, double pointY
) {
  double distance{ std::hypot(pointX - centerX, pointY - centerY) };
  double normDistance{ distance / radius };
  return std::exp(-5.0 * std::pow(normDistance, 2.0));
}

//Equivalent R code
//gauss <- function(w, CtrX, CtrY, PointX, PointY){
//  distance <- ((PointX-CtrX)^2+(PointY-CtrY)^2)^0.5
//  norm.distance <- distance/w
//  output <- exp(-5*norm.distance^2)
//  return(output)
//}
