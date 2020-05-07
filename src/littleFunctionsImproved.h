#ifndef LITTLE_FUNCTIONS_IMPROVED_H
#define LITTLE_FUNCTIONS_IMPROVED_H

bool intersectsCylinder(
    const double pointX, const double pointY, const double pointZ,
    const double cylinderRadius, const double cylinderHeight,
    const double centerX, const double centerY, const double centerZ
);


double gauss(const double x);

double epanechnikov(const double x);


double calculateVerticalWeight(
    const double pointZ, const double cylinderMiddleZ, const double cylinderHeight
);

double calculateHorizontalWeight(
    const double pointX, const double pointY,
    const double cylinderRadius, const double centerX, const double centerY
);

#endif  // define LITTLE_FUNCTIONS_IMPROVED_H
