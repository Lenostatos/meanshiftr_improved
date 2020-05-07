#include "littleFunctionsImproved.h"

#include <cmath>  // for std::pow, std::abs, std::hypot, std::exp


bool intersectsCylinder(
    const double pointX, const double pointY, const double pointZ,
    const double cylinderRadius, const double cylinderHeight,
    const double centerX, const double centerY, const double centerZ
) {
  double topHeight{ centerZ + (0.5 * cylinderHeight) };
  double bottomHeight{ topHeight - cylinderHeight };

  return std::pow((pointX - centerX), 2.0) + std::pow((pointY - centerY), 2.0)
            <= std::pow(cylinderRadius, 2.0)
    && bottomHeight <= pointZ && pointZ <= topHeight;
}


double calculateVerticalWeight(
    const double pointZ, const double cylinderMiddleZ, const double cylinderHeight
) {
  double relativeVerticalDistanceToCenter{
    std::abs(cylinderMiddleZ - pointZ) / (cylinderHeight * 0.5)
  };
  return epanechnikov(relativeVerticalDistanceToCenter);
}

double calculateHorizontalWeight(
    const double pointX, const double pointY,
    const double cylinderRadius, const double centerX, const double centerY
) {
  double relativeHorizontalDistanceToCenter{
    std::hypot(centerX - pointX, centerY - pointY) / cylinderRadius
  };
  return gauss(relativeHorizontalDistanceToCenter);
}


/** The gaussian function f(x) = exp(gaussian_gamma * x^2).
 *
 *  Analogous to equation (11) in Ferraz et al. 2012.
 */
double gauss(const double x) {
  return std::exp(-5.0 * std::pow(x, 2));
}

/** The epanechnikov distribution function f(x) = 1 - x^2 .
 *
 *  Analogous to parts of equation (14) in Ferraz et al. 2012.
 */
double epanechnikov(const double x) {
  return 1 - std::pow(x, 2);
}
