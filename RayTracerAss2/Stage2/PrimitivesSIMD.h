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

// points consist of three coordinates and represent a point in 3d space
typedef struct Point8
{
	__m256 x, y, z;

	inline Point8() {}

	inline Point8(float x1, float y1, float z1) { this->x = _mm256_set1_ps(x1);  this->y = _mm256_set1_ps(y1);  this->z = _mm256_set1_ps(z1); }

	inline Point8(__m256 x1, __m256 y1, __m256 z1) { x = x1; y = y1; z = z1; }

	inline __m256 length()
	{
		__m256 xx = _mm256_mul_ps(x, x);
		__m256 yy = _mm256_mul_ps(y, y);
		__m256 zz = _mm256_mul_ps(y, y);
		__m256 total = _mm256_add_ps(xx, _mm256_add_ps(yy, zz));
		//return sqrtf(x * x + y * y + z * z);
		return _mm256_sqrt_ps(total);
	}
} Point8;

// Represent 8 vectors in one struct
struct Vector8
{
	__m256 xs, ys, zs;

	__forceinline Vector8(__m256 x) {
		xs = x;
		ys = x;
		zs = x;
	}

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

	inline __m256 dot() const
	{
		return xs * xs + ys * ys + zs * zs;
	}
};

typedef struct Ray8
{
	Point8 start;
	Vector8 dir;

}Ray8;


// helper operators / functions for Vector8s
__forceinline Vector8 operator - (const Vector8& v1, const Vector8& v2) { return { v1.xs - v2.xs, v1.ys - v2.ys, v1.zs - v2.zs }; }
__forceinline Vector8 operator + (const Vector8& v1, const Vector8& v2) { return { v1.xs + v2.xs, v1.ys + v2.ys, v1.zs + v2.zs }; }
__forceinline Vector8 operator * (const Vector8& v1, const Vector8& v2) { return { v1.xs * v2.xs, v1.ys * v2.ys, v1.zs * v2.zs }; }

__forceinline Vector8 cross(const Vector8& v1, const Vector8& v2)
{
	return { v1.ys * v2.zs - v1.zs * v2.ys, v1.zs * v2.xs - v1.xs * v2.zs, v1.xs * v2.ys - v1.ys * v2.xs };
}

// point - vector (produces a point)
inline Vector8 operator - (const Point8& p1, const Point8& p2)
{
	Vector8 v = { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
	return v;
}

// point - vector (produces a point)
inline Point8 operator - (const Point8& p, const Vector8& v)
{
	Point8 p2( p.x - v.xs, p.y - v.ys, p.z - v.zs );
	return p2;
}

// point / float (each component of the point is divided by the float)
inline Point8 operator / (const Point8& p, __m256 c)
{

	Point8 p2( p.x / c, p.y / c, p.z / c );
	return p2;
}

// helper function
inline __m256 invsqrtf(const __m256& x)
{
	return _mm256_set1_ps( 1.0f ) / _mm256_sqrt_ps(x);
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


#endif

