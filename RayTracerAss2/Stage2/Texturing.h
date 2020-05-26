#ifndef __TEXTURING_H
#define __TEXTURING_H

#include "Scene.h"
#include "Intersection.h"

// apply computed checkerboard texture
Colour8 applyCheckerboard(const Intersection* intersect);

// apply computed turbulence texture
Colour8 applyCircles(const Intersection* intersect);

// apply computed wood texture
Colour8 applyWood(const Intersection* intersect);

#endif // __TEXTURING_H
