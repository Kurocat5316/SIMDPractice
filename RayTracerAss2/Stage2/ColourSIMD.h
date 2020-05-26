/*  The following code is a VERY heavily modified from code originally sourced from:
Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#ifndef __COLOUR_S
#define __COLOUR_S

#include <algorithm>
#include "PrimitivesSIMD.h"

// a colour consists of three primary components (red, green, and blue)
struct Colour8
{
	__m256 red, green, blue;

	// create uninitialised colour
	inline Colour8() { }

	// create initialised colour
	inline Colour8(__m256 r, __m256 g, __m256 b) : red(r), green(g), blue(b)
	{
	}

	inline Colour8(float r, float g, float b) : red(_mm256_set1_ps(r)), green(_mm256_set1_ps(g)), blue(_mm256_set1_ps(b))
	{
	}

	// create colour from integer (pixel stored in 0x00BBGGRR format)
	inline Colour8(unsigned int value) :
		red(_mm256_set1_ps((value & 0xFF) / 255.0f)),
		green(_mm256_set1_ps(((value >> 8) & 0xFF) / 255.0f)),
		blue(_mm256_set1_ps(((value >> 16) & 0xFF) / 255.0f))
	{
	}

	// colourise the RGB values based on a mask 
	inline void colourise(unsigned int colourMask)
	{
		red = red * _mm256_set1_ps((colourMask & 4) ? 1.5f : 0.75f);
		green = green * _mm256_set1_ps((colourMask & 2) ? 1.5f : 0.75f);
		blue = blue *  _mm256_set1_ps((colourMask & 1) ? 1.5f : 0.75f);
	}



	// colour += colour
	inline Colour8& operator += (const Colour8& c2)
	{
		this->red = red + c2.red;
		this->green = green + c2.green;
		this->blue = blue + c2.blue;
		return *this;
	}
};

// colour * colour
inline Colour8 operator * (const Colour8& c1, const Colour8& c2)
{
	return Colour8(c1.red * c2.red, c1.green * c2.green, c1.blue * c2.blue);
}

// colour + colour
inline Colour8 operator + (const Colour8& c1, const Colour8& c2)
{
	return Colour8(c1.red + c2.red, c1.green + c2.green, c1.blue + c2.blue);
}

// float * colour (float multiplied to each colour channel)
inline Colour8 operator * (float coef, const Colour8& c1)
{
	__m256 coefTmp = _mm256_set1_ps(coef);
	return Colour8(c1.red * coefTmp, c1.green * coefTmp, c1.blue * coefTmp);
}

inline Colour8 operator * (__m256 coef, const Colour8& c1)
{
	return Colour8(c1.red * coef, c1.green * coef, c1.blue * coef);
}

#endif //__COLOUR_S

