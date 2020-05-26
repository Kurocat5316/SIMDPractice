/*  The following code is a VERY heavily modified from code originally sourced from:
	Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
	It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#ifndef __LIGHTING_H
#define __LIGHTING_H

#include "Scene.h"
#include "Intersection.h"

// test to see if light ray collides with any of the scene's objects
__m256 isInShadow(const Scene* scene, const Ray* lightRay, const float lightDist);

// apply diffuse lighting with respect to material's colouring
Colour8 applyDiffuse(const Ray* lightRay, const Light* currentLight, const Intersection* intersect);

// apply specular lighting using Blinn
Colour8 applySpecular(const Ray8* lightRay, const Light8* currentLight, const __m256 fLightProjection, const Ray8* viewRay, const Intersection* intersect);

// apply diffuse and specular lighting contributions for all lights in scene taking shadowing into account
Colour applyLighting(const Scene* scene, const Ray* viewRay, const Intersection* intersect); 


#endif // __LIGHTING_H
