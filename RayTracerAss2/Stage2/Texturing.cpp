#include "Texturing.h"
#include "Colour.h"
#include "Intersection.h"
#include "PrimitivesSIMD.h"
#include "MathSIMD.h"

// apply computed checkerboard texture
Colour8 applyCheckerboard(const Intersection* intersect)
{
	//Point p = (intersect->pos - intersect->material->offset) / intersect->material->size;

	Point8 pos(intersect->pos.x, intersect->pos.y, intersect->pos.z);
	Vector8 offset(intersect->material->offset.x, intersect->material->offset.y, intersect->material->offset.z);
	__m256 size = _mm256_set1_ps( intersect->material->size);

	Point8 p = (pos - offset) / size;

	__m256 which = (_mm256_floor_ps(p.x) + _mm256_floor_ps(p.y) + _mm256_floor_ps(p.z)) & _mm256_set1_ps(1);

	//int which = (int(floorf(p.x)) + int(floorf(p.y)) + int(floorf(p.z))) & 1;

	Colour8 color;
	color.red = select(which, _mm256_set1_ps(intersect->material->diffuse.red), _mm256_set1_ps(intersect->material->diffuse2.red));
	color.green = select(which, _mm256_set1_ps(intersect->material->diffuse.green), _mm256_set1_ps(intersect->material->diffuse2.green));
	color.blue = select(which, _mm256_set1_ps(intersect->material->diffuse.blue), _mm256_set1_ps(intersect->material->diffuse2.blue));

	return color;
}


// apply computed circular texture
Colour8 applyCircles(const Intersection* intersect)
{
	//Point p = (intersect->pos - intersect->material->offset) / intersect->material->size;

	Point8 pos(intersect->pos.x, intersect->pos.y, intersect->pos.z);
	Vector8 offset(intersect->material->offset.x, intersect->material->offset.y, intersect->material->offset.z);
	__m256 size = _mm256_set1_ps(intersect->material->size);

	Point8 p = (pos - offset) / size;

	//int which = int(floorf(sqrtf(p.x*p.x + p.y*p.y + p.z*p.z))) & 1;

	__m256 which = (_mm256_floor_ps(_mm256_sqrt_ps(p.x*p.x + p.y*p.y + p.z*p.z))) & _mm256_set1_ps(1);

	//return (which ? intersect->material->diffuse : intersect->material->diffuse2);

	Colour8 color;
	color.red = select(which, _mm256_set1_ps(intersect->material->diffuse.red), _mm256_set1_ps(intersect->material->diffuse2.red));
	color.green = select(which, _mm256_set1_ps(intersect->material->diffuse.green), _mm256_set1_ps(intersect->material->diffuse2.green));
	color.blue = select(which, _mm256_set1_ps(intersect->material->diffuse.blue), _mm256_set1_ps(intersect->material->diffuse2.blue));

	return color;
}


// apply computed wood grain texture
Colour8 applyWood(const Intersection* intersect)
{
	//Point p = (intersect->pos - intersect->material->offset) / intersect->material->size;

	Point8 pos(intersect->pos.x, intersect->pos.y, intersect->pos.z);
	Vector8 offset(intersect->material->offset.x, intersect->material->offset.y, intersect->material->offset.z);
	__m256 size = _mm256_set1_ps(intersect->material->size);

	Point8 p = (pos - offset) / size;


	// squiggle up where the point is
	//p = { p.x * cosf(p.y * 0.996f) * sinf(p.z * 1.023f), cosf(p.x) * p.y * sinf(p.z * 1.211f), cosf(p.x * 1.473f) * cosf(p.y * 0.795f) * p.z };

	p.x = p.x * cos256_ps(p.y * _mm256_set1_ps(0.996f)) * sin256_ps(p.z * _mm256_set1_ps(1.023f));
	p.y = cos256_ps(p.x) * p.y * sin256_ps(p.z * _mm256_set1_ps(1.211f));
	p.z = cos256_ps(p.x * _mm256_set1_ps(1.473f)) * cos256_ps(p.y * _mm256_set1_ps(0.795f)) *p.z;

	//p = { p.x * cos256_ps(p.y * _mm256_set1_ps(0.996f)) * sin256_ps(p.z * _mm256_set1_ps(1.023f)) };

	//int which = int(floorf(sqrtf(p.x*p.x + p.y*p.y + p.z*p.z))) & 1;

	__m256 which = (_mm256_floor_ps(_mm256_sqrt_ps(p.x*p.x + p.y*p.y + p.z*p.z))) & _mm256_set1_ps(1);

	//return (which ? intersect->material->diffuse : intersect->material->diffuse2);

	Colour8 color;
	color.red = select(which, _mm256_set1_ps(intersect->material->diffuse.red), _mm256_set1_ps(intersect->material->diffuse2.red));
	color.green = select(which, _mm256_set1_ps(intersect->material->diffuse.green), _mm256_set1_ps(intersect->material->diffuse2.green));
	color.blue = select(which, _mm256_set1_ps(intersect->material->diffuse.blue), _mm256_set1_ps(intersect->material->diffuse2.blue));

	return color;
}
