/*  The following code is a VERY heavily modified from code originally sourced from:
Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#ifndef __SCENE_H
#define __SCENE_H

#include "SceneObjects.h"
#include <immintrin.h>

// description of a single static scene
typedef struct Scene
{
	Point cameraPosition;					// camera location
	float cameraRotation;					// direction camera points
	float cameraFieldOfView;				// field of view for the camera

	float exposure;							// image exposure

	unsigned int skyboxMaterialId;

	// scene object counts
	unsigned int numMaterials;
	unsigned int numSpheres;
	unsigned int numTriangles;
	unsigned int numLights;

	// scene objects
	Material* materialContainer;
	Sphere* sphereContainer;
	Triangle* triangleContainer;
	Light* lightContainer;

	// SIMD spheres
	unsigned int numSpheresSIMD;
	__m256* spherePosX;
	__m256* spherePosY;
	__m256* spherePosZ;
	__m256* sphereSize;
	__m256i* sphereMaterialId;

	// SIMD triangles
	unsigned int numTriganleSIMD;
	__m256* p1X;
	__m256* p1Y;
	__m256* p1Z;
	__m256* p2X;
	__m256* p2Y;
	__m256* p2Z;
	__m256* p3X;
	__m256* p3Y;
	__m256* p3Z;
	__m256* normalX;
	__m256* normalY;
	__m256* normalZ;
	__m256i* triangleMaterialId;
} Scene;

bool init(const char* inputName, Scene& scene);

#endif // __SCENE_H
