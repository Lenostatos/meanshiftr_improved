#include <Rcpp.h>
using namespace Rcpp;

// Declarations of all the little functions used by the main functions

//' Is a point within a cylinder?
//'
//' Indicates whether a point \code{[x, y, z]} lies within a vertical cylider
//' that is defined by a radius, a height, and a center
//' (\code{[centerX, centerY, centerZ]}).
//'
//' @param x,y,z The coordinates of the point.
//' @param radius The radius of the cylinder.
//' @param height The height of the cylinder.
//' @param centerX,centerY,centerZ The center of the cylinder.
//'
//' @return \code{true} if point lies within cylinder, \code{false} otherwise.
bool InCylinder(
    double x, double y, double z,
    double radius, double height,
    double centerX, double centerY, double centerZ
);

//' Relative vertical distance to cylinder boundary
//'
//' Calculates the relative vertical position of a point between the center and
//' the closest outer boundary of a cylinder that is defined as the upper three
//' quarters of a vertical cylinder with height \code{height} and a vertical
//' center at \code{centerZ}.
double VerticalDistance(double height, double centerZ, double pointZ);

//' 1-0 mask for the upper three quarters of a cylinder
//'
//' Returns 1 if point lies within the upper three quarters of the vertical
//' cylinder defined by its \code{height} and vertical center \code{centerZ},
//' 0 otherwise.
short VerticalMask(double height, double centerZ, double pointZ);

//' Epanechnikov weighting of a point's vertical position in a cylinder
//'
//' Returns \eqn{ 1 - x^2 } where x is the vertical distance of \code{pointZ} to
//' the center of a cylinder that is the three upper quarters of a vertical
//' cylinder with height \code{height} and vertical center \code{centerZ}. The
//' distance is normalized with half of the three-quarter cylinder's height.
//' Returns 0, if point does not lie within the three upper quarters.
double EpanechnikovFunction(double height, double centerZ, double pointZ);

//' Gauss weighting of a point's horizontal position in a cylinder
//'
//' Returns \eqn{ 5 * x^2 } where x is the horizontal distance of point
//' (\code{pointX, pointY}) to the center (\code{centerX, centerY}) relative to
//' the cylinder's radius.
double GaussFunction(
    double radius, double centerX, double centerY, double pointX, double pointY
);
