/*  The following code is a VERY heavily modified from code originally sourced from:
	Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
	It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#include "Intersection.h"
#include <immintrin.h>
#include "PrimitivesSIMD.h"


// helper function to find "horizontal" minimum (and corresponding index value from another vector)
__forceinline void selectMinimumAndIndex(__m256 values, __m256i indexes, float* min, int* index)
{
	// find min of elements 1&2, 3&4, 5&6, and 7&8
	__m256 minNeighbours = _mm256_min_ps(values, _mm256_permute_ps(values, 0x31));
	// find min of min(1,2)&min(5,6) and min(3,4)&min(7,8)
	__m256 minNeighbours2 = _mm256_min_ps(minNeighbours, _mm256_permute2f128_ps(minNeighbours, minNeighbours, 0x05));
	// find final minimum 
	__m256 mins = _mm256_min_ps(minNeighbours2, _mm256_permute_ps(minNeighbours2, 0x02));

	// find all elements that match our minimum
	__m256i matchingTs = _mm256_castps_si256(_mm256_set1_ps(mins.m256_f32[0]) != values);
	// set all other elements to be MAX_INT (-1 but unsigned)
	__m256i matchingIndexes = matchingTs | indexes;

	// find minimum of remaining indexes (so smallest index will be chosen) using that same technique as above but with heaps of ugly casts
	__m256i minIndexNeighbours = _mm256_min_epu32(matchingIndexes, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(matchingIndexes), 0x31)));
	__m256i minIndexNeighbours2 = _mm256_min_epu32(minIndexNeighbours, _mm256_castps_si256(_mm256_permute2f128_ps(
		_mm256_castsi256_ps(minIndexNeighbours), _mm256_castsi256_ps(minIndexNeighbours), 0x05)));
	__m256i minIndex = _mm256_min_epu32(minIndexNeighbours2, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(minIndexNeighbours2), 0x02)));

	// "return" minimum and associated index through reference parameters
	*min = mins.m256_f32[0];
	*index = minIndex.m256i_i32[0];
}

// test to see if collision between ray and a plane happens before time t (equivalent to distance)
// updates closest collision time (/distance) if collision occurs
// see: http://en.wikipedia.org/wiki/Line-sphere_intersection
// see: http://www.codermind.com/articles/Raytracer-in-C++-Part-I-First-rays.html
// see: Step 8 of http://meatfighter.com/juggler/ 
// this code make heavy use of constant term removal due to ray always being a unit vector
bool isSphereIntersected(const Scene* scene, const Ray* r, float* t, int* index)
{
	float tInitial = *t;

	// ray start and direction
	Vector8 rStart(r->start.x, r->start.y, r->start.z);
	Vector8 rDir(r->dir.x, r->dir.y, r->dir.z);

	// constants
	const __m256 epsilons = _mm256_set1_ps(EPSILON);
	const __m256 zeros = _mm256_set1_ps(0.0f);
	const __m256i eights = _mm256_set1_epi32(8);

	// best ts found so far and associated triangle indexes
	__m256 ts = _mm256_set1_ps(tInitial);
	__m256i indexes = _mm256_set1_epi32(*index);

	// current corresponding index
	__m256i ijs = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

	// search for sphere collisions, storing closest one found
	for (unsigned int i = 0; i < scene->numSpheresSIMD; ++i)
	{
		//Sphere& sphere = scene.sphereContainer[i * 8 + j];
		Vector8 pos(scene->spherePosX[i], scene->spherePosY[i], scene->spherePosZ[i]);
		__m256 sizes = scene->sphereSize[i];

		// Vector dist = pos - r.start;
		Vector8 dist = pos - rStart;

		// float B = r.dir * dist;
		__m256 Bs = dot(rDir, dist);

		// float D = B * B - dist * dist + size * size;
		__m256 Ds = Bs * Bs - dot(dist, dist) + sizes * sizes;

		// if D < 0, no intersection, so don't try and calculate the point of intersection
		//if (D < 0.0f) continue;
		__m256 DLessThanZeros = Ds < zeros;

		// calculate both intersection times(/distances)
		//float t0 = B - sqrtf(D);
		//float t1 = B + sqrtf(D);
		__m256 sqrtDs = _mm256_sqrt_ps(Ds);
		__m256 t0s = Bs - sqrtDs;
		__m256 t1s = Bs + sqrtDs;

		// check to see if either of the two sphere collision points are closer than time parameter
		//if ((t1 > EPSILON) && (t1 < t))
		__m256 t1GreaterThanEpsilonAndSmallerThanTs = (t1s > epsilons) & (t1s < ts);
		//else if ((t0 > EPSILON) && (t0 < t))
		__m256 t0GreaterThanEpsilonAndSmallerThanTs = (t0s > epsilons) & (t0s < ts);

		// select best ts 
		__m256 temp = select(t1GreaterThanEpsilonAndSmallerThanTs, t1s, ts);
		__m256 temp2 = select(t0GreaterThanEpsilonAndSmallerThanTs, t0s, temp);
		ts = select(DLessThanZeros, ts, temp2); 

		// select best corresponding sphere indexes
		__m256i temp3 = select(_mm256_castps_si256(t1GreaterThanEpsilonAndSmallerThanTs), ijs, indexes);
		__m256i temp4 = select(_mm256_castps_si256(t0GreaterThanEpsilonAndSmallerThanTs), ijs, temp3);
		indexes = select(_mm256_castps_si256(DLessThanZeros), indexes, temp4);

		// increase the index counters
		ijs = _mm256_add_epi32(ijs, eights);
	}

	// extract the best t and corresponding triangle index
	selectMinimumAndIndex(ts, indexes, t, index);

	return *t < tInitial;
}


// short-circuiting version of sphere intersection test that only returns true/false
bool isSphereIntersected(const Scene* scene, const Ray* r, float t)
{
	// ray start and direction
	Vector8 rStart(r->start.x, r->start.y, r->start.z);
	Vector8 rDir(r->dir.x, r->dir.y, r->dir.z);

	// constants
	const __m256 epsilons = _mm256_set1_ps(EPSILON);
	const __m256 zeros = _mm256_set1_ps(0.0f);

	// starting t
	const __m256 ts = _mm256_set1_ps(t);

	// search for sphere collisions, storing closest one found
	for (unsigned int i = 0; i < scene->numSpheresSIMD; ++i)
	{
		//Sphere& sphere = scene.sphereContainer[i * 8 + j];
		Vector8 pos(scene->spherePosX[i], scene->spherePosY[i], scene->spherePosZ[i]);
		__m256 sizes = scene->sphereSize[i];

		// Vector dist = pos - r.start;
		Vector8 dist = pos - rStart; 

		// float B = r.dir * dist;
		__m256 Bs = dot(rDir, dist);
		
		// float D = B * B - dist * dist + size * size;
		__m256 Ds = Bs * Bs - dot(dist, dist) + sizes * sizes;

		// if D < 0, no intersection, so don't try and calculate the point of intersection
		//if (D < 0.0f) continue;
		__m256 DLessThanZeros = Ds < zeros;

		// calculate both intersection times(/distances)
		//float t0 = B - sqrtf(D);
		//float t1 = B + sqrtf(D);
		__m256 sqrtDs = _mm256_sqrt_ps(Ds);
		__m256 t0s = Bs - sqrtDs;
		__m256 t1s = Bs + sqrtDs;

		// check to see if either of the two sphere collision points are closer than time parameter
		//if ((t1 > EPSILON) && (t1 < t))
		__m256 t1GreaterThanEpsilonAndSmallerThanTs = (t1s > epsilons) & (t1s < ts);
		//else if ((t0 > EPSILON) && (t0 < t))
		__m256 t0GreaterThanEpsilonAndSmallerThanTs = (t0s > epsilons) & (t0s < ts);

		// combine all the success cases together
		__m256 success = _mm256_andnot_ps(DLessThanZeros, t0GreaterThanEpsilonAndSmallerThanTs | t1GreaterThanEpsilonAndSmallerThanTs);

		// if any are successful, short-circuit
		if (_mm256_movemask_ps(success)) return true;
	}

	return false;
}


// short-circuiting version of sphere intersection test that only returns true/false
// this version doesn't use helper functions from PrimitivesSIMD.h and is only included for reference
/*bool isSphereIntersected(const Scene* scene, const Ray* r, float t)
{
	// ray start and direction
	const __m256 rStartxs = _mm256_set1_ps(r->start.x);
	const __m256 rStartys = _mm256_set1_ps(r->start.y);
	const __m256 rStartzs = _mm256_set1_ps(r->start.z);
	const __m256 rDirxs = _mm256_set1_ps(r->dir.x);
	const __m256 rDirys = _mm256_set1_ps(r->dir.y);
	const __m256 rDirzs = _mm256_set1_ps(r->dir.z);

	// constants
	const __m256 epsilons = _mm256_set1_ps(EPSILON);
	const __m256 zeros = _mm256_set1_ps(0.0f);

	// starting t
	const __m256 ts = _mm256_set1_ps(t);

	// search for sphere collisions, storing closest one found
	for (unsigned int i = 0; i < scene->numSpheresSIMD; ++i)
	{
		//Sphere& sphere = scene.sphereContainer[i * 8 + j];
		__m256 posxs = scene->spherePosX[i];
		__m256 posys = scene->spherePosY[i];
		__m256 poszs = scene->spherePosZ[i];
		__m256 sizes = scene->sphereSize[i];

		// Vector dist = pos - r.start;
		__m256 distxs = _mm256_sub_ps(posxs, rStartxs);
		__m256 distys = _mm256_sub_ps(posys, rStartys);
		__m256 distzs = _mm256_sub_ps(poszs, rStartzs);

		// float B = r.dir * dist;
		__m256 Bs = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rDirxs, distxs), _mm256_mul_ps(rDirys, distys)), _mm256_mul_ps(rDirzs, distzs));

		// float D = B * B - dist * dist + size * size;
		__m256 distDot = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(distxs, distxs), _mm256_mul_ps(distys, distys)), _mm256_mul_ps(distzs, distzs));
		__m256 Ds = _mm256_add_ps(_mm256_sub_ps(_mm256_mul_ps(Bs, Bs), distDot), _mm256_mul_ps(sizes, sizes));

		// if D < 0, no intersection, so don't try and calculate the point of intersection
		//if (D < 0.0f) continue;
		__m256 DLessThanZeros = _mm256_cmp_ps(Ds, zeros, _CMP_LT_OQ);

		// calculate both intersection times(/distances)
		//float t0 = B - sqrtf(D);
		//float t1 = B + sqrtf(D);
		__m256 sqrtDs = _mm256_sqrt_ps(Ds);
		__m256 t0s = _mm256_sub_ps(Bs, sqrtDs);
		__m256 t1s = _mm256_add_ps(Bs, sqrtDs);

		// check to see if either of the two sphere collision points are closer than time parameter
		//if ((t1 > EPSILON) && (t1 < t))
		__m256 t0GreaterThanEpsilons = _mm256_cmp_ps(t0s, epsilons, _CMP_GT_OQ);
		__m256 t0SmallerThanTs = _mm256_cmp_ps(t0s, ts, _CMP_LT_OQ);
		__m256 t0GreaterThanEpsilonAndSmallerThanTs = _mm256_and_ps(t0GreaterThanEpsilons, t0SmallerThanTs);
		//else if ((t0 > EPSILON) && (t0 < t))
		__m256 t1GreaterThanEpsilons = _mm256_cmp_ps(t1s, epsilons, _CMP_GT_OQ);
		__m256 t1SmallerThanTs = _mm256_cmp_ps(t1s, ts, _CMP_LT_OQ);
		__m256 t1GreaterThanEpsilonAndSmallerThanTs = _mm256_and_ps(t1GreaterThanEpsilons, t1SmallerThanTs);

		// combine all the success cases together
		__m256 success = _mm256_andnot_ps(DLessThanZeros,
			_mm256_or_ps(t0GreaterThanEpsilonAndSmallerThanTs, t1GreaterThanEpsilonAndSmallerThanTs));

		// if any are successful, short-circuit
		if (_mm256_movemask_ps(success)) return true;
	}

	return false;
}*/


// test to see if collision between ray and a triangle happens before time t (equivalent to distance)
// updates closest collision time (/distance) if collision occurs
// based on: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
// explanation at: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
// another: http://hugi.scene.org/online/hugi25/hugi%2025%20-%20coding%20corner%20graphics,%20sound%20&%20synchronization%20ken%20ray-triangle%20intersection%20tests%20for%20dummies.htm
bool isTriangleIntersected(const Triangle* tri, const Ray* r, float* t)
{

	// two edges of the triangle (as world coordinates offsets)
	Vector e1 = tri->p2 - tri->p1;
	Vector e2 = tri->p3 - tri->p1;

	// vector perpendicular to the ray's direction and the second edge of the triangle
	Vector h = cross(r->dir, e2);

	// determinant of above vector and the first edge of the triangle
	float det = e1 * h;

	// no intersection if this value is small (i.e. ray is parallel with triangle surface)
	if (det > -EPSILON && det < EPSILON) return false;
	
	float invDet = 1.0f / det;

	// distance vector between start of ray and first point of triangle
	Vector s = r->start - tri->p1;

	// barycentric coord u (i.e. p2 in coordinate system based around triangle's extents)
	float u = invDet * (s * h);

	// no intersection if value is "before" or "after" the triangle
	if (u < 0.0f || u > 1.0f) return false;

	// barycentric coord v (i.e. p3)
	Vector q = cross(s, e1);
	float v = invDet * (q * r->dir);

	// no intersection if value is "before" or "after" the triangle
	if (v < 0.0f || u + v > 1.0f) return false;

	// find point of intersection
	float t0 = invDet * (e2 * q);

	// check to see if triangle collision point is closer than time parameter
	if (t0 > EPSILON && t0 < *t)
	{
		*t = t0;
		return true;
	}

	return false;
}

bool isTriangleIntersected(const Scene* scene, const Ray* r, float* t, int* num)
{
	__m256 dirX = _mm256_set1_ps(r->dir.x);
	__m256 dirY = _mm256_set1_ps(r->dir.y);
	__m256 dirZ = _mm256_set1_ps(r->dir.z);

	Vector8 dir(dirX, dirY, dirZ);

	__m256 startX = _mm256_set1_ps(r->start.x);
	__m256 startY = _mm256_set1_ps(r->start.y);
	__m256 startZ = _mm256_set1_ps(r->start.z);
	
	__m256 tt = _mm256_set1_ps(*t);

	__m256 zero = _mm256_setzero_ps();
	__m256 one = _mm256_set1_ps(1.0f);

	__m256 record = _mm256_setzero_ps();

	bool flag = false;

	for (unsigned int index = 0; index < scene->numTriganleSIMD; ++index) {

		// two edges of the triangle (as world coordinates offsets)
		//Vector e1;
		//e1.x = scene->p2X[index / 8].m256_f32[index % 8] - scene->p1X[index / 8].m256_f32[index % 8];
		//e1.y = scene->p2Y[index / 8].m256_f32[index % 8] - scene->p1Y[index / 8].m256_f32[index % 8];
		//e1.z = scene->p2Z[index / 8].m256_f32[index % 8] - scene->p1Z[index / 8].m256_f32[index % 8];
		//Vector e2;
		//e2.x = scene->p3X[index / 8].m256_f32[index % 8] - scene->p1X[index / 8].m256_f32[index % 8];
		//e2.y = scene->p3Y[index / 8].m256_f32[index % 8] - scene->p1Y[index / 8].m256_f32[index % 8];
		//e2.z = scene->p3Z[index / 8].m256_f32[index % 8] - scene->p1Z[index / 8].m256_f32[index % 8];
		float test = scene->p1X[index].m256_f32[0];

		__m256 x1 = _mm256_sub_ps(scene->p2X[index], scene->p1X[index]);
		__m256 y1 = _mm256_sub_ps(scene->p2Y[index], scene->p1Y[index]);
		__m256 z1 = _mm256_sub_ps(scene->p2Z[index], scene->p1Z[index]);
		Vector8 e1(x1, y1, z1);

		__m256 x2 = _mm256_sub_ps(scene->p3X[index], scene->p1X[index]);
		__m256 y2 = _mm256_sub_ps(scene->p3Y[index], scene->p1Y[index]);
		__m256 z2 = _mm256_sub_ps(scene->p3Z[index], scene->p1Z[index]);
		Vector8 e2(x2, y2, z2);


		// vector perpendicular to the ray's direction and the second edge of the triangle
		//Vector h = cross(r->dir, e2);
		

		Vector8 h = cross(dir, e2);


		// determinant of above vector and the first edge of the triangle
		//float det = e1 * h;

		__m256 det = VectorMul(e1, h);

		// no intersection if this value is small (i.e. ray is parallel with triangle surface)
		__m256 epstlon = _mm256_set1_ps(-EPSILON);
		__m256 epstlon2 = _mm256_set1_ps(EPSILON);

		//if (det > -EPSILON && det < EPSILON) return false;
		__m256 judgement1 = _mm256_cmp_ps(det, epstlon, _CMP_GT_OQ);
		__m256 judgement2 = _mm256_cmp_ps(det, epstlon, _CMP_LT_OQ);
		__m256 judgementDet = _mm256_and_ps(judgement1, judgement2);

		//float invDet = 1.0f / det;
		__m256 invDet = _mm256_div_ps(_mm256_set1_ps(1.0f), det);

		// distance vector between start of ray and first point of triangle
		//Point temp1{ scene->p1X[index / 8].m256_f32[index % 8], scene->p1Y[index / 8].m256_f32[index % 8], scene->p1Z[index / 8].m256_f32[index % 8] };

		//Vector s = r->start - temp1;

		

		__m256 sX = _mm256_sub_ps(startX, scene->p1X[index]);
		__m256 sY = _mm256_sub_ps(startY, scene->p1Y[index]);
		__m256 sZ = _mm256_sub_ps(startZ, scene->p1Z[index]);

		Vector8 s(sX, sY, sZ);

		// barycentric coord u (i.e. p2 in coordinate system based around triangle's extents)
		//float u = invDet * (s * h);
		__m256 u = _mm256_mul_ps(invDet, VectorMul(s, h));

		// no intersection if value is "before" or "after" the triangle
		//if (u < 0.0f || u > 1.0f) return false;

		__m256 judgement3 = _mm256_cmp_ps(u, zero, _CMP_LT_OQ);
		__m256 judgement4 = _mm256_cmp_ps(u, one, _CMP_GT_OQ);
		__m256 judgementU = _mm256_or_ps(judgement3, judgement4);


		// barycentric coord v (i.e. p3)
		//Vector q = cross(s, e1);
		//float v = invDet * (q * r->dir);

		Vector8 q = cross(s, e1);
		__m256 v = _mm256_mul_ps(invDet, VectorMul(q, dir));

		// no intersection if value is "before" or "after" the triangle
		//if (v < 0.0f || u + v > 1.0f) return false;
		__m256 judgement5 = _mm256_cmp_ps(v, zero, _CMP_LT_OQ);
		__m256 judgement6 = _mm256_cmp_ps(_mm256_add_ps(u, v), one, _CMP_GT_OQ);
		__m256 judgementV = _mm256_or_ps(judgement5, judgement6);

		// find point of intersection
		//float t0 = invDet * (e2 * q);
		__m256 t0 = _mm256_mul_ps(invDet, VectorMul(e2, q));

		// check to see if triangle collision point is closer than time parameter
		/*if (t0 > EPSILON && t0 < *t)
		{
			*t = t0;
			return true;
		}

		return false;*/
		

		//__m256 judgement7 = _mm256_cmp_ps(t0, epstlon2, _CMP_GT_OQ);
		//__m256 judgement8 = _mm256_cmp_ps(t0, tt, _CMP_LE_OQ);
		//__m256 judgementT0 = _mm256_and_ps(judgement7, judgement8);

		__m256 judgement7 = _mm256_cmp_ps(t0, epstlon2, _CMP_LE_OQ);
		__m256 judgement8 = _mm256_cmp_ps(t0, tt, _CMP_GE_OQ);
		__m256 judgementT0 = _mm256_or_ps(judgement7, judgement8);

		__m256 cond1 = _mm256_or_ps(judgementDet, judgementU);
		__m256 cond2 = _mm256_or_ps(cond1, judgementV);
		__m256 cond3 = _mm256_or_ps(cond2, judgementT0);

		__m256 cond4 = _mm256_andnot_ps(cond3, t0);
		__m256 cond5 = _mm256_and_ps(cond3, tt);
		tt = _mm256_or_ps(cond4, cond5);

		__m256 cond6 = _mm256_and_ps(cond3, record);
		__m256 cond7 = _mm256_andnot_ps(cond3, _mm256_set1_ps(index));

		record = _mm256_or_ps(cond6, cond7);

		//get min
		__m256 v0 = tt;                /* _mm256_shuffle_ps instead of _mm256_permute_ps is also possible, see Peter Cordes' comments */
		__m256 v1 = _mm256_permute_ps(v0, 0b10110001); /* swap floats: 0<->1, 2<->3, 4<->5, 6<->7                         */
		__m256 v2 = _mm256_min_ps(v0, v1);
		__m256 v3 = _mm256_permute_ps(v2, 0b01001110); /* swap floats                                                     */
		__m256 v4 = _mm256_min_ps(v2, v3);
		__m256 v5 = _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(v4), 0b01001110)); /* swap 128-bit lanes      */
		__m256 v_min = _mm256_min_ps(v4, v5);
		__m256 mask = _mm256_cmp_ps(v0, v_min, 0);
		int indx = _tzcnt_u32(_mm256_movemask_ps(mask));

		//return min *t and index 
		if (v_min.m256_f32[0] < *t) {
			*t = v_min.m256_f32[0];
			int indx = _tzcnt_u32(_mm256_movemask_ps(mask));
			*num = indx + (int)record.m256_f32[indx] * 8;
			*num = *num < scene->numTriangles ? *num : scene->numTriangles - 1;
			flag = true;
		}
	}

	if (flag)
		return true;

	return false;

}


// calculate collision normal, viewProjection, object's material, and test to see if inside collision object
void calculateIntersectionResponse(const Scene* scene, const Ray* viewRay, Intersection* intersect)
{
	switch (intersect->objectType)
	{
	case Intersection::SPHERE:
		intersect->normal = normalise(intersect->pos - intersect->sphere->pos);
		intersect->material = &scene->materialContainer[intersect->sphere->materialId];
		break;
	case Intersection::TRIANGLE:
		intersect->normal = intersect->triangle->normal;
		intersect->material = &scene->materialContainer[intersect->triangle->materialId];
	}

	// calculate view projection
	intersect->viewProjection = viewRay->dir * intersect->normal; 

	// detect if we are inside an object (needed for refraction)
	intersect->insideObject = (intersect->normal * viewRay->dir > 0.0f);

	// if inside an object, reverse the normal
    if (intersect->insideObject)
    {
        intersect->normal = intersect->normal * -1.0f;
    }
}


// test to see if collision between ray and any object in the scene
// updates intersection structure if collision occurs
bool objectIntersection(const Scene* scene, const Ray* viewRay, Intersection* intersect)
{
	// set default distance to be a long long way away
    float t = MAX_RAY_DISTANCE;

	// no intersection found by default
	intersect->objectType = Intersection::NONE;

	// search for sphere collisions, storing closest one found
	int index = -1;
	if (isSphereIntersected(scene, viewRay, &t, &index))
	{
		intersect->objectType = Intersection::SPHERE;
		intersect->sphere = &scene->sphereContainer[index];
	}

	//unsigned int record = 1;
	//float recordT;
	//int flag = false;
	// search for triangle collisions, storing closest one found
	//for (unsigned int i = 0; i < scene->numTriangles; ++i)
	//{
	//
	//	if (isTriangleIntersected(&scene->triangleContainer[i], viewRay, &t))
	//	//if (isTriangleIntersected(scene, viewRay, &t, i, intersect))
	//	{
	//		intersect->objectType = Intersection::TRIANGLE;
	//		intersect->triangle = &scene->triangleContainer[i];
	//		record = i;
	//		recordT = t;
	//		flag = !flag;
	//	}
	//}


	int i = 0;
	if (isTriangleIntersected(scene, viewRay, &t, &i))
	{
		//i++;
		intersect->objectType = Intersection::TRIANGLE;
		intersect->triangle = &scene->triangleContainer[i];
		
	}




	// nothing detected, return false
	if (intersect->objectType == Intersection::NONE)
	{
		return false;
	}

	// calculate the point of the intersection
	intersect->pos = viewRay->start + viewRay->dir * t;

	return true;
}
