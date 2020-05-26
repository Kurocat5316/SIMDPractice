// SIMD versions of (some of the) helper functions from primitives.h

//NOTE: the only code that uses this at the moment is the non-short-circuiting isTriangleIntersected function

#ifndef __PRIMITIVES_SIMD_H
#define __PRIMITIVES_SIMD_H

#include <immintrin.h>

// a bunch of operators to replace nasty instrinsics
__forceinline __m256 operator - (const __m256 x, const __m256 y) { return _mm256_sub_ps(x, y); }
__forceinline __m256 operator + (const __m256 x, const __m256 y) { return _mm256_add_ps(x, y); }
__forceinline __m256 operator * (const __m256 x, const __m256 y) { return _mm256_mul_ps(x, y); }
__forceinline __m256 operator / (const __m256 x, const __m256 y) { return _mm256_div_ps(x, y); }
__forceinline __m256 operator < (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_LT_OQ); }
__forceinline __m256 operator > (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_GT_OQ); }
__forceinline __m256 operator <= (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_LE_OQ); }
__forceinline __m256 operator >= (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_GE_OQ); }
__forceinline __m256 operator & (const __m256 x, const __m256 y) { return _mm256_and_ps(x, y); }
__forceinline __m256 operator | (const __m256 x, const __m256 y) { return _mm256_or_ps(x, y); }
__forceinline __m256 operator == (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_EQ_OQ); }
__forceinline __m256 operator != (const __m256 x, const __m256 y) { return _mm256_cmp_ps(x, y, _CMP_NEQ_OQ); }


__forceinline __m256i operator & (const __m256i x, const __m256i y) { return _mm256_and_si256(x, y); }
__forceinline __m256i operator | (const __m256i x, const __m256i y) { return _mm256_or_si256(x, y); }


// Represent 8 vectors in one struct
struct Vector8
{
	__m256 xs, ys, zs;

	__forceinline Vector8(float x, float y, float z)
	{
		xs = _mm256_set1_ps(x);
		ys = _mm256_set1_ps(y);
		zs = _mm256_set1_ps(z);
	}

	__forceinline Vector8(__m256 xsIn, __m256 ysIn, __m256 zsIn)
	{
		xs = xsIn;
		ys = ysIn;
		zs = zsIn;
	}

	__forceinline Vector8()
	{
		xs = ys = zs = _mm256_setzero_ps();
	}
};


// helper operators / functions for Vector8s
__forceinline Vector8 operator - (const Vector8& v1, const Vector8& v2) { return { v1.xs - v2.xs, v1.ys - v2.ys, v1.zs - v2.zs }; }
__forceinline Vector8 operator + (const Vector8& v1, const Vector8& v2) { return { v1.xs + v2.xs, v1.ys + v2.ys, v1.zs + v2.zs }; }
__forceinline Vector8 operator * (const Vector8& v1, const Vector8& v2) { return { v1.xs * v2.xs, v1.ys * v2.ys, v1.zs * v2.zs }; }

__forceinline Vector8 cross(const Vector8& v1, const Vector8& v2)
{
	return { v1.ys * v2.zs - v1.zs * v2.ys, v1.zs * v2.xs - v1.xs * v2.zs, v1.xs * v2.ys - v1.ys * v2.xs };
}

__forceinline __m256 dot(const Vector8& v1, const Vector8& v2)
{
	return v1.xs * v2.xs + v1.ys * v2.ys + v1.zs * v2.zs;
}

__forceinline __m256 select(__m256 cond, __m256 ifTrue, __m256 ifFalse)
{
	return _mm256_or_ps(_mm256_and_ps(cond, ifTrue), _mm256_andnot_ps(cond, ifFalse));
}

__forceinline __m256i select(__m256i cond, __m256i ifTrue, __m256i ifFalse)
{
	return _mm256_or_si256(_mm256_and_si256(cond, ifTrue), _mm256_andnot_si256(cond, ifFalse));
}

__forceinline __m256 VectorMul(Vector8 vector1, Vector8 vector2) {
	__m256 tmpX = _mm256_mul_ps(vector1.xs, vector2.xs);
	__m256 tmpY = _mm256_mul_ps(vector1.ys, vector2.ys);
	__m256 tmpZ = _mm256_mul_ps(vector1.zs, vector2.zs);
	__m256 distance = _mm256_add_ps(tmpX, _mm256_add_ps(tmpY, tmpZ));
	return distance;
}
#endif

