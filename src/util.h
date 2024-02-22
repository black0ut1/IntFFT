#pragma once

#include <cmath>
#include <complex>
#include <cstdint>

using complex = std::complex<double>;
using namespace std::complex_literals;

constexpr double NEG2PI = -2 * M_PI;

inline complex cexp(double theta) {
	return {cos(theta), sin(theta)};
}

/** Normalize angle to interval (-pi, pi) */
double normalize(double theta) {
	theta = fmod(theta, 2 * M_PI);

	if (theta >= M_PI)
		theta -= 2 * M_PI;
	else if (theta <= -M_PI)
		theta += 2 * M_PI;

	return theta;
}

/** std::complex does not support integer types */
struct complex_int {
	int32_t r;
	int32_t i;

	complex_int operator+(complex_int that) const {
		return {r + that.r, i + that.i};
	}

	complex_int operator-(complex_int that) const {
		return {r - that.r, i - that.i};
	}

	void operator+=(complex_int that) {
		r += that.r;
		i += that.i;
	}

	complex_int times1i() {
		return {-i, r};
	}

	int32_t real() const {
		return r;
	}

	int32_t imag() const {
		return i;
	}
};
