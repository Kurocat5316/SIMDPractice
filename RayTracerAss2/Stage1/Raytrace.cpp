/*  The following code is a VERY heavily modified from code originally sourced from:
Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#define TARGET_WINDOWS

#pragma warning(disable: 4996)
#include "Timer.h"
#include "Primitives.h"
#include "Scene.h"
#include "Lighting.h"
#include "Intersection.h"
#include "ImageIO.h"

unsigned int buffer[MAX_WIDTH * MAX_HEIGHT];

// reflect the ray from an object
Ray calculateReflection(const Ray* viewRay, const Intersection* intersect)
{
	// reflect the viewRay around the object's normal
	Ray newRay = { intersect->pos, viewRay->dir - (intersect->normal * intersect->viewProjection * 2.0f) };

	return newRay;
}


// refract the ray through an object
Ray calculateRefraction(const Ray* viewRay, const Intersection* intersect, float* currentRefractiveIndex)
{
	// change refractive index depending on whether we are in an object or not
	float oldRefractiveIndex = *currentRefractiveIndex;
	*currentRefractiveIndex = intersect->insideObject ? DEFAULT_REFRACTIVE_INDEX : intersect->material->density;

	// calculate refractive ratio from old index and current index
	float refractiveRatio = oldRefractiveIndex / *currentRefractiveIndex;

	// Here we take into account that the light movement is symmetrical from the observer to the source or from the source to the oberver.
	// We then do the computation of the coefficient by taking into account the ray coming from the viewing point.
	float fCosThetaT;
	float fCosThetaI = fabsf(intersect->viewProjection);

	// glass-like material, we're computing the fresnel coefficient.
	if (fCosThetaI >= 1.0f)
	{
		// In this case the ray is coming parallel to the normal to the surface
		fCosThetaT = 1.0f;
	}
	else
	{
		float fSinThetaT = refractiveRatio * sqrtf(1 - fCosThetaI * fCosThetaI);

		// Beyond the angle (1.0f) all surfaces are purely reflective
		fCosThetaT = (fSinThetaT * fSinThetaT >= 1.0f) ? 0.0f : sqrtf(1 - fSinThetaT * fSinThetaT);
	}

	// Here we compute the transmitted ray with the formula of Snell-Descartes
	Ray newRay = { intersect->pos, (viewRay->dir + intersect->normal * fCosThetaI) * refractiveRatio - (intersect->normal * fCosThetaT) };

	return newRay;
}


// follow a single ray until it's final destination (or maximum number of steps reached)
Colour traceRay(const Scene* scene, Ray viewRay)
{
	Colour output(0.0f, 0.0f, 0.0f); 								// colour value to be output
	float currentRefractiveIndex = DEFAULT_REFRACTIVE_INDEX;		// current refractive index
	float coef = 1.0f;												// amount of ray left to transmit
	Intersection intersect;											// properties of current intersection

																	// loop until reached maximum ray cast limit (unless loop is broken out of)
	for (int level = 0; level < MAX_RAYS_CAST; ++level)
	{
		// check for intersections between the view ray and any of the objects in the scene
		// exit the loop if no intersection found
		if (!objectIntersection(scene, &viewRay, &intersect)) break;

		// calculate response to collision: ie. get normal at point of collision and material of object
		calculateIntersectionResponse(scene, &viewRay, &intersect);

		// apply the diffuse and specular lighting 
		if (!intersect.insideObject) output += coef * applyLighting(scene, &viewRay, &intersect);

		// if object has reflection or refraction component, adjust the view ray and coefficent of calculation and continue looping
		if (intersect.material->reflection)
		{
			viewRay = calculateReflection(&viewRay, &intersect);
			coef *= intersect.material->reflection;
		}
		else if (intersect.material->refraction)
		{
			viewRay = calculateRefraction(&viewRay, &intersect, &currentRefractiveIndex);
			coef *= intersect.material->refraction;
		}
		else
		{
			// if no reflection or refraction, then finish looping (cast no more rays)
			return output;
		}
	}

	// if the calculation coefficient is non-zero, read from the environment map
	if (coef > 0.0f)
	{
		Material& currentMaterial = scene->materialContainer[scene->skyboxMaterialId];

		output += coef * currentMaterial.diffuse;
	}

	return output;
}

// render a section of the scene at given width and height and anti-aliasing level
void renderSection(Scene* scene, const int width, const int height, const int aaLevel, const int blockSize, unsigned int* out, const unsigned int colourMask, unsigned int* currentBlockShared)
{
	// angle between each successive ray cast (per pixel, anti-aliasing uses a fraction of this)
	const float dirStepSize = 1.0f / (0.5f * width / tanf(PIOVER180 * 0.5f * scene->cameraFieldOfView));

	// calculate exactly how many blocks are needed (and deal with cases where the blockSize doesn't exactly divide)
	unsigned int blocksWide = (width - 1) / blockSize + 1;
	unsigned int blocksHigh = (height - 1) / blockSize + 1;
	unsigned int blocksTotal = blocksWide * blocksHigh;

	// current block index
	unsigned int currentBlock;

	while ((currentBlock = InterlockedIncrement(currentBlockShared)) < blocksTotal)
	{
		// block x,y position
		const int bx = currentBlock % blocksWide;
		const int by = currentBlock / blocksWide;

		// block coordinates (making sure not to exceed image bounds with non-divisible block sizes)
		const int xMin = bx * blockSize - width / 2;
		const int xMax = (std::min)(xMin + blockSize, width / 2);
		const int yMin = by * blockSize - height / 2;
		const int yMax = (std::min)(yMin + blockSize, height / 2);

		// calculate the array location of the start of the block
		unsigned int* outBlock = buffer + width * by * blockSize + bx * blockSize;

		// jump required to get to the start of the next line of the block
		unsigned int outJump = width - (xMax - xMin);

		// loop through all the pixels
		for (int y = yMin; y < yMax; ++y)
		{
			for (int x = xMin; x < xMax; ++x)
			{
				Colour output(0.0f, 0.0f, 0.0f);

				// calculate multiple samples for each pixel
				const float sampleStep = 1.0f / aaLevel, sampleRatio = 1.0f / (aaLevel * aaLevel);

				// loop through all sub-locations within the pixel
				for (float fragmentx = float(x); fragmentx < x + 1.0f; fragmentx += sampleStep) //1.0f / aaLevel)
				{
					for (float fragmenty = float(y); fragmenty < y + 1.0f; fragmenty += sampleStep) //1.0f / aaLevel)
					{
						// direction of default forward facing ray
						Vector dir = { fragmentx * dirStepSize, fragmenty * dirStepSize, 1.0f };

						// rotated direction of ray
						Vector rotatedDir = {
							dir.x * cosf(scene->cameraRotation) - dir.z * sinf(scene->cameraRotation),
							dir.y,
							dir.x * sinf(scene->cameraRotation) + dir.z * cosf(scene->cameraRotation) };

						// view ray starting from camera position and heading in rotated (normalised) direction
						Ray viewRay = { scene->cameraPosition, normalise(rotatedDir) };

						// follow ray and add proportional of the result to the final pixel colour
						output += sampleRatio * traceRay(scene, viewRay);
					}
				}

				// colour the pixel
				output.colourise(colourMask);

				// store saturated final colour value in image buffer
				*outBlock++ = output.convertToPixel(scene->exposure);
			}

			// move to the start of the next line of the block
			outBlock += outJump;
		}
	}
}


// data for threads
struct ThreadParams
{
	Scene* scene;
	int width;
	int height;
	int aaLevel;
	int blockSize;
	unsigned int* out;
	unsigned int colourMask;
	unsigned int* currentBlockShared;
};


// thread callback for rendering
DWORD __stdcall renderSectionThread(LPVOID inData)
{
	// cast the void* structure to something useful
	ThreadParams* params = (ThreadParams*)inData;

	// call the real render function
	renderSection(params->scene, params->width, params->height, params->aaLevel, params->blockSize, params->out, params->colourMask, params->currentBlockShared);

	// exit with success
	ExitThread(NULL);
}


// render scene at given width and height and anti-aliasing level using a specified number of threads
void render(Scene* scene, const int width, const int height, const int aaLevel, const unsigned int threadCount, const int blockSize, const bool colourise)
{
	// reserve space for threads and their parameters
	HANDLE* threads = new HANDLE[threadCount];
	ThreadParams* params = new ThreadParams[threadCount];

	// one less than the current block to render (shared between threads)
	unsigned int currentBlockShared = -1;

	// loop through all the squares
	for (unsigned int i = 0; i < threadCount; ++i)
	{
		// debug
		//printf("thread rendering: [%d,%d] (%d,%d)->(%d,%d) => %x\n", bx, by, xMin, yMin, xMax, yMax, out);

		// set up thread parameters
		params[i] = { scene, width, height, aaLevel, blockSize, buffer, colourise ? (i % 8) : 7, &currentBlockShared };

		// start thread
		threads[i] = CreateThread(NULL, 0, renderSectionThread, (LPVOID)&params[i], 0, NULL);
	}

	// wait until all the threads are done
	if (threadCount <= 64)
	{
		WaitForMultipleObjects(threadCount, threads, TRUE, INFINITE);
	}
	else
	{
		for (unsigned int i = 0; i < threadCount; i++) {
			WaitForSingleObject(threads[i], INFINITE);
		}
	}

	// clean up thread and param storage
	delete[] params;
	delete[] threads;
}



// allocate space fro SoA, and copy values from AoS to SoA 
void simdifySceneContainers(Scene& scene)
{
	// helper size (so we don't just have 8 everywhere)
	unsigned int valuesPerVector = sizeof(__m256) / sizeof(float);

	// make SoA SIMD copies of spheres (if there are any)
	if (scene.numSpheres == 0)
	{
		scene.numSpheresSIMD = 0;
	}
	else
	{
		// mathemagical way of calculating ceilf(scene.numSpheres / 8.0f)
		scene.numSpheresSIMD = (((int)scene.numSpheres) - 1) / valuesPerVector + 1;

		// allocate the correct amount of space at the correct alignment for SIMD operations
		scene.spherePosX = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numSpheresSIMD, 32);
		scene.spherePosY = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numSpheresSIMD, 32);
		scene.spherePosZ = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numSpheresSIMD, 32);
		scene.sphereSize = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numSpheresSIMD, 32);
		scene.sphereMaterialId = (__m256i*) _aligned_malloc(sizeof(__m256i) * scene.numSpheresSIMD, 32);

		// initialise SoA structures
		for (unsigned int i = 0; i < scene.numSpheresSIMD * valuesPerVector; ++i)
		{
			// don't let the source index extend out of the AoS array
			// i.e. copy the last value into the extra array slots when numSpheres isn't exactly divisible by 8
			// pretty lazy way to fix this, but it works
			int sourceIndex = i < scene.numSpheres ? i : scene.numSpheres - 1;

			scene.spherePosX[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.sphereContainer[sourceIndex].pos.x;
			scene.spherePosY[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.sphereContainer[sourceIndex].pos.y;
			scene.spherePosZ[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.sphereContainer[sourceIndex].pos.z;
			scene.sphereSize[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.sphereContainer[sourceIndex].size;
			scene.sphereMaterialId[i / valuesPerVector].m256i_i32[i % valuesPerVector] = scene.sphereContainer[sourceIndex].materialId;
		}
	}


	if (scene.numTriangles == 0)
		scene.numTriganleSIMD = 0;
	else
	{
		scene.numTriganleSIMD = (((int)scene.numTriangles) - 1) / valuesPerVector + 1;

		// allocate the correct amount of space at the correct alignment for SIMD operations
		scene.p1X = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p1Y = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p1Z = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p2X = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p2Y = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p2Z = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p3X = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p3Y = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.p3Z = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.normalX = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.normalY = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.normalZ = (__m256*) _aligned_malloc(sizeof(__m256) * scene.numTriganleSIMD, 32);
		scene.triangleMaterialId = (__m256i*) _aligned_malloc(sizeof(__m256i) * scene.numTriganleSIMD, 32);

		for (unsigned int i = 0; i < scene.numTriganleSIMD * valuesPerVector; ++i)
		{
			// don't let the source index extend out of the AoS array
			// i.e. copy the last value into the extra array slots when numSpheres isn't exactly divisible by 8
			// pretty lazy way to fix this, but it works
			int sourceIndex = i < scene.numTriangles ? i : scene.numTriangles - 1;

			scene.p1X[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p1.x;
			scene.p1Y[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p1.y;
			scene.p1Z[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p1.z;

			scene.p2X[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p2.x;
			scene.p2Y[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p2.y;
			scene.p2Z[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p2.z;

			scene.p3X[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p3.x;
			scene.p3Y[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p3.y;
			scene.p3Z[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].p3.z;

			scene.normalX[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].normal.x;
			scene.normalY[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].normal.y;
			scene.normalZ[i / valuesPerVector].m256_f32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].normal.z;

			scene.triangleMaterialId[i / valuesPerVector].m256i_i32[i % valuesPerVector] = scene.triangleContainer[sourceIndex].materialId;
		}
	}

}



// read command line arguments, render, and write out BMP file
int main(int argc, char* argv[])
{
	int width = 256;
	int height = 256;
	int samples = 1;

	// rendering options
	int times = 5;
	unsigned int threads = 8;			
	bool colourise = false;				
	unsigned int blockSize = 8;		

	// default input / output filenames
	const char* inputFilename = "../Scenes/bunny10k.txt";
	char outputFilenameBuffer[1000];
	char* outputFilename = outputFilenameBuffer;

	// do stuff with command line args
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-size") == 0)
		{
			width = atoi(argv[++i]);
			height = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-samples") == 0)
		{
			samples = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-input") == 0)
		{
			inputFilename = argv[++i];
		}
		else if (strcmp(argv[i], "-output") == 0)
		{
			outputFilename = argv[++i];
		}
		else if (strcmp(argv[i], "-runs") == 0)
		{
			times = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-threads") == 0)
		{
			threads = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-colourise") == 0)
		{
			colourise = true;
		}
		else if (strcmp(argv[i], "-blockSize") == 0)
		{
			blockSize = atoi(argv[++i]);
		}
		else
		{
			fprintf(stderr, "unknown argument: %s\n", argv[i]);
		}
	}

	// nasty (and fragile) kludge to make an ok-ish default output filename (can be overriden with "-output" command line option)
	sprintf(outputFilenameBuffer, "../Outputs/%s_%dx%dx%d_%s.bmp", (strrchr(inputFilename, '/') + 1), width, height, samples, (strrchr(argv[0], '\\') + 1));

	// read scene file
	Scene scene;
	if (!init(inputFilename, scene))
	{
		fprintf(stderr, "Failure when reading the Scene file.\n");
		return -1;
	}

	// do the SoA things
	simdifySceneContainers(scene);

	// total time taken to render all runs (used to calculate average)
	int totalTime = 0;
	for (int i = 0; i < times; i++)
	{
		Timer timer;															// create timer
		render(&scene, width, height, samples, threads, blockSize, colourise);	// raytrace scene
		timer.end();															// record end time
		totalTime += timer.getMilliseconds();									// record total time taken
	}

	// output timing information (times run and average)
	printf("average time taken (%d run(s)): %ums\n", times, totalTime / times);

	// output BMP file
	write_bmp(outputFilename, buffer, width, height, width);
}
