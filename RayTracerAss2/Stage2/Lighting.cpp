/*  The following code is a VERY heavily modified from code originally sourced from:
Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#include "Lighting.h"
#include "Colour.h"
#include "Intersection.h"
#include "Texturing.h"
#include "PrimitivesSIMD.h"

// test to see if light ray collides with any of the scene's objects
// short-circuits when first intersection discovered, because no matter what the object will be in shadow
__m256 isInShadow(const Scene* scene, const Ray8* lightRay, const __m256 lightDist)
{
	//initial record to all true (-nan)
	__m256 tmp = _mm256_set1_ps(true);
	__m256 tmp2 = _mm256_set1_ps(false);
	__m256 record = _mm256_cmp_ps(tmp, tmp2, _CMP_GT_OQ);

	for (int z = 0; z < 8; z++) {
		Ray lightRayT;
		lightRayT.start.x = lightRay->start.x.m256_f32[z];
		lightRayT.start.y = lightRay->start.y.m256_f32[z];
		lightRayT.start.z = lightRay->start.z.m256_f32[z];

		lightRayT.dir.x = lightRay->dir.xs.m256_f32[z];
		lightRayT.dir.y = lightRay->dir.ys.m256_f32[z];
		lightRayT.dir.z = lightRay->dir.zs.m256_f32[z];

		const Ray* lightRay2 = &lightRayT;

		float t = lightDist.m256_f32[z];

		// search for sphere collision
		if (isSphereIntersected(scene, lightRay2, t)) record.m256_f32[z] = false;

		// search for triangle collision
		for (unsigned int i = 0; i < scene->numTriangles; ++i)
		{
			if (isTriangleIntersected(&scene->triangleContainer[i], lightRay2, &t))
			{
				record.m256_f32[z] = false;
			}
		}
	}

	// not in shadow
	return record;
}


// apply diffuse lighting with respect to material's colouring
Colour8 applyDiffuse(const Ray8* lightRay, const Light8* currentLight, const Intersection* intersect)
{
	Colour8 output;

	

	//double i = Material::CHECKERBOARD;

	for (int i = 0; i < 8; i++) {

			if (intersect->material->type == Material::GOURAUD) {
				output = Colour8(intersect->material->diffuse.red, intersect->material->diffuse.green, intersect->material->diffuse.blue);
				break;
			}
			if (intersect->material->type == Material::CHECKERBOARD) {
				output = applyCheckerboard(intersect);
				break;
			}
			if (intersect->material->type == Material::CIRCLES) {
				output = applyCircles(intersect);
				break;
			}
			if (intersect->material->type == Material::WOOD) {
				output = applyWood(intersect);
				break;
			}

	}

	Vector8 normal;
	normal.xs = _mm256_set1_ps(intersect->normal.x);
	normal.ys = _mm256_set1_ps(intersect->normal.y);
	normal.zs = _mm256_set1_ps(intersect->normal.z);

	Vector8 tmp = lightRay->dir * normal;
	__m256 lambert = tmp.xs + tmp.ys + tmp.zs;

	return lambert * currentLight->intensity * output;
}


// Blinn 
// The direction of Blinn is exactly at mid point of the light ray and the view ray. 
// We compute the Blinn vector and then we normalize it then we compute the coeficient of blinn
// which is the specular contribution of the current light.
Colour8 applySpecular(const Ray8* lightRay, const Light8* currentLight, const __m256 fLightProjection, const Ray* viewRay, const Intersection* intersect)
{
	Vector8 viewRayDir(viewRay->dir.x, viewRay->dir.y, viewRay->dir.z);

	//Vector blinnDir = lightRay->dir - viewRay->dir;
	Vector8 blinnDir = lightRay->dir - viewRayDir;

	//float blinn = invsqrtf(blinnDir.dot()) * std::max(fLightProjection - intersect->viewProjection, 0.0f);
	__m256 blinn = invsqrtf(blinnDir.dot()) * _mm256_max_ps(fLightProjection - _mm256_set1_ps(intersect->viewProjection), _mm256_set1_ps(0.0f));


	//blinn = powf(blinn, intersect->material->power);
	for (int i = 0; i < 8; i++) {
		blinn.m256_f32[i] = powf(blinn.m256_f32[i], intersect->material->power);
	}
	
	float test = 8;

	Colour8 specular(intersect->material->specular.red, intersect->material->specular.green, intersect->material->specular.blue);
	return blinn * specular * currentLight->intensity;
}


// apply diffuse and specular lighting contributions for all lights in scene taking shadowing into account
Colour applyLighting(const Scene* scene, const Ray* viewRay, const Intersection* intersect)
{
	// colour to return (starts as black)
	Colour output(0.0f, 0.0f, 0.0f);

	// same starting point for each light ray
	//Ray lightRay = { intersect->pos };
	Point8 intersectPos;
	intersectPos.x = _mm256_set1_ps(intersect->pos.x);
	intersectPos.y = _mm256_set1_ps(intersect->pos.y);
	intersectPos.z = _mm256_set1_ps(intersect->pos.z);

	Vector8 normal;
	normal.xs = _mm256_set1_ps(intersect->normal.x);
	normal.ys = _mm256_set1_ps(intersect->normal.y);
	normal.zs = _mm256_set1_ps(intersect->normal.z);

	Ray8 lightRay;
	lightRay.start = intersectPos;

	// loop through all the lights
	//for (unsigned int j = 0; j < scene->numLights; ++j)
	for(unsigned int j = 0; j < scene->numLightsSIMD; ++j)
	{
		Colour8 output2(0.0f, 0.0f, 0.0f);
		// get reference to current light
		//const Light* currentLight = &scene->lightContainer[j];
		Light8* currentLight1 = new Light8();
		currentLight1->pos.x = scene->lightPosX[j];
		currentLight1->pos.y = scene->lightPosY[j];
		currentLight1->pos.z = scene->lightPosZ[j];

		currentLight1->intensity = Colour8( *scene->intensityR, *scene->intensityG, *scene->intensityB);

		// light ray direction need to equal the normalised vector in the direction of the current light
		// as we need to reuse all the intermediate components for other calculations, 
		// we calculate the normalised vector by hand instead of using the normalise function
		//lightRay.dir = currentLight->pos - intersect->pos;
		lightRay.dir = currentLight1->pos - intersectPos;

		//float angleBetweenLightAndNormal = lightRay.dir * intersect->normal;
		Vector8 tmp = lightRay.dir * normal;
		__m256 angleBetweenLightAndNormal = tmp.xs + tmp.ys + tmp.zs;


		// skip this light if it's behind the object (ie. both light and normal pointing in the same direction)
		//if (angleBetweenLightAndNormal <= 0.0f)
		//{
		//	continue;
		//}

		__m256 judgement1 = _mm256_cmp_ps(angleBetweenLightAndNormal, _mm256_set1_ps(0.0f), _CMP_GT_OQ);

		// distance to light from intersection point (and it's inverse)
		//float lightDist = sqrtf(lightRay.dir.dot());
		__m256 lightDist = _mm256_sqrt_ps( lightRay.dir.dot());

		__m256 invLightDist = _mm256_set1_ps( 1.0f) / lightDist;

		// light ray projection
		__m256 lightProjection = invLightDist * angleBetweenLightAndNormal;

		// normalise the light direction
		lightRay.dir = lightRay.dir * Vector8(invLightDist);

		// only apply lighting from this light if not in shadow of some other object
		__m256 judgement2 = isInShadow(scene, &lightRay, lightDist);

		/*if (!isInShadow(scene, &lightRay, lightDist))
		{
			// add diffuse lighting from colour / texture
			output += applyDiffuse(&lightRay, currentLight, intersect);

			// add specular lighting
			output += applySpecular(&lightRay, currentLight, lightProjection, viewRay, intersect);
		}*/

		__m256 judgement = _mm256_and_ps(judgement1, judgement2);

		Colour8 tmpCol = applyDiffuse(&lightRay, currentLight1, intersect);
		Colour8 tmpCol2 = applySpecular(&lightRay, currentLight1, lightProjection, viewRay, intersect);

		output2.red = select(judgement1, output2.red + tmpCol.red, output2.red);
		output2.red = select(judgement, output2.red + tmpCol2.red, output2.red);

		output2.green = select(judgement, output2.green + tmpCol.green, output2.green);
		output2.green = select(judgement, output2.green + tmpCol2.green, output2.green);

		output2.blue = select(judgement, output2.blue + tmpCol.blue, output2.blue);
		output2.blue = select(judgement, output2.blue + tmpCol2.blue, output2.blue);

		for (int i = 0; i < 8; i++)
		{
			if (i + j * 8 >= scene->numLights)
				break;
			output.red += output2.red.m256_f32[i];
			output.green += output2.green.m256_f32[i];
			output.blue += output2.blue.m256_f32[i];
		}
	}

	

	return output;
}
