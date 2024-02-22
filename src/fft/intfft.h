#pragma once

#include "../util.h"
#include <vector>

template<double (*Q)(double)>
complex_int lift(complex_int x, double theta) {
	if (theta == 0)
		return x;

	double c = cos(theta);
	double s = sin(theta);
	int32_t re = x.r;
	int32_t im = x.i;

	if (0 < theta && theta < M_PI) { // Figure 5(c)
		im = -im;
		re += (int32_t) Q(im * (s - 1) / c);
		im += (int32_t) Q(re * c);
		re += (int32_t) Q(im * (s - 1) / c);

		return {im, re};
	} else if (-M_PI < theta && theta < 0) { // Figure 5(d)
		im = -im;
		re += (int32_t) Q(im * (s + 1) / c);
		im += (int32_t) Q(re * -c);
		re += (int32_t) Q(im * (s + 1) / c);

		return {-im, -re};
	}

	throw std::runtime_error("Theta is outside interval (-pi, pi)");
}

//////////// DFT ////////////

template<double (*Q)(double)>
std::vector<complex_int> IntDFT(std::vector<complex_int> &x) {
	int N = (int) x.size();
	std::vector<complex_int> X(x.size());

	for (int k = 0; k < N; ++k) {

		complex_int sum{};
		for (int n = 0; n < N; ++n) {
			double theta = normalize(NEG2PI * k * n / N);

			sum += lift<Q>(x, theta);
		}

		X[k] = sum;
	}

	return X;
}

//////////// Radix 2 FFT ////////////

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_radix2_DIT(const std::vector<complex_int>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex_int> y(N / 2);
	std::vector<complex_int> z(N / 2);
	for (int k = 0; k < N / 2; ++k) {
		y[k] = x[2 * k	  ];
		z[k] = x[2 * k + 1];
	}

	auto Y = IntFFT_radix2_DIT<Q>(y);
	auto Z = IntFFT_radix2_DIT<Q>(z);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 2; ++k) {
		double theta = normalize(NEG2PI * k / N);
		complex_int a = lift<Q>(Z[k], theta); // tf * Z[k]

		X[k]		 = Y[k] + a;
		X[k + N / 2] = Y[k] - a;
	}

	return X;
}

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_radix2_DIF(const std::vector<complex_int>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex_int> y(N / 2);
	std::vector<complex_int> z(N / 2);
	for (int k = 0; k < N / 2; ++k) {
		double theta = normalize(NEG2PI * k / N);

		complex_int a = x[k] + x[k + N / 2];
		complex_int b = x[k] - x[k + N / 2];

		y[k] = a;
		z[k] = lift<Q>(b, theta);
	}

	auto Y = IntFFT_radix2_DIF<Q>(y);
	auto Z = IntFFT_radix2_DIF<Q>(z);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 2; ++k) {
		X[2 * k]	 = Y[k];
		X[2 * k + 1] = Z[k];
	}

	return X;
}

//////////// Radix 4 FFT ////////////

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_radix4_DIT(const std::vector<complex_int>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex_int> y(N / 4);
	std::vector<complex_int> z(N / 4);
	std::vector<complex_int> g(N / 4);
	std::vector<complex_int> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		y[k] = x[4 * k	  ];
		z[k] = x[4 * k + 1];
		g[k] = x[4 * k + 2];
		h[k] = x[4 * k + 3];
	}

	auto Y = IntFFT_radix4_DIT<Q>(y);
	auto Z = IntFFT_radix4_DIT<Q>(z);
	auto G = IntFFT_radix4_DIT<Q>(g);
	auto H = IntFFT_radix4_DIT<Q>(h);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 4; ++k) {
		double theta1 = normalize(	  NEG2PI * k / N);
		double theta2 = normalize(2 * NEG2PI * k / N);
		double theta3 = normalize(3 * NEG2PI * k / N);
		complex_int a = lift<Q>(Z[k], theta1);
		complex_int b = lift<Q>(G[k], theta2);
		complex_int c = lift<Q>(H[k], theta3);

		X[k] 			 = (Y[k] + b) + (a + c);
		X[k + N / 4]	 = (Y[k] - b) - (a - c).times1i();
		X[k + N / 2]	 = (Y[k] + b) - (a + c);
		X[k + 3 * N / 4] = (Y[k] - b) + (a - c).times1i();
	}

	return X;
}

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_radix4_DIF(const std::vector<complex_int>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex_int> y(N / 4);
	std::vector<complex_int> z(N / 4);
	std::vector<complex_int> g(N / 4);
	std::vector<complex_int> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		double theta1 = normalize(	  NEG2PI * k / N);
		double theta2 = normalize(2 * NEG2PI * k / N);
		double theta3 = normalize(3 * NEG2PI * k / N);

		complex_int a = (x[k] + x[k + N / 2]) + (x[k + N / 4] + x[k + 3 * N / 4]);
		complex_int b = (x[k] - x[k + N / 2]) - (x[k + N / 4] - x[k + 3 * N / 4]).times1i();
		complex_int c = (x[k] + x[k + N / 2]) - (x[k + N / 4] + x[k + 3 * N / 4]);
		complex_int d = (x[k] - x[k + N / 2]) + (x[k + N / 4] - x[k + 3 * N / 4]).times1i();

		y[k] = a;
		z[k] = lift<Q>(b, theta1);
		g[k] = lift<Q>(c, theta2);
		h[k] = lift<Q>(d, theta3);
	}

	auto Y = IntFFT_radix4_DIF<Q>(y);
	auto Z = IntFFT_radix4_DIF<Q>(z);
	auto G = IntFFT_radix4_DIF<Q>(g);
	auto H = IntFFT_radix4_DIF<Q>(h);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 4; ++k) {
		X[4 * k	   ] = Y[k];
		X[4 * k + 1] = Z[k];
		X[4 * k + 2] = G[k];
		X[4 * k + 3] = H[k];
	}

	return X;
}

//////////// Split radix FFT ////////////

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_splitRadix_DIT(const std::vector<complex_int>& x) {
	auto N = (int32_t) x.size();
	if (N == 1)
		return x;
	if (N == 2)
		return {x[0] + x[1], x[0] - x[1]};

	std::vector<complex_int> y(N / 2);
	std::vector<complex_int> z(N / 4);
	std::vector<complex_int> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		y[2 * k] 	 = x[4 * k	  ];
		z[k] 		 = x[4 * k + 1];
		y[2 * k + 1] = x[4 * k + 2];
		h[k] 		 = x[4 * k + 3];
	}

	auto Y = IntFFT_splitRadix_DIT<Q>(y);
	auto Z = IntFFT_splitRadix_DIT<Q>(z);
	auto H = IntFFT_splitRadix_DIT<Q>(h);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 4; ++k) {
		double theta1 = normalize(	  NEG2PI * k / N);
		double theta3 = normalize(3 * NEG2PI * k / N);
		complex_int a = lift<Q>(Z[k], theta1);
		complex_int b = lift<Q>(H[k], theta3);

		X[k] 			 = Y[k] 		+ (a + b);
		X[k + N / 4] 	 = Y[k + N / 4] - (a - b).times1i();
		X[k + N / 2] 	 = Y[k] 		- (a + b);
		X[k + 3 * N / 4] = Y[k + N / 4] + (a - b).times1i();
	}

	return X;
}

template<double (*Q)(double)>
std::vector<complex_int> IntFFT_splitRadix_DIF(const std::vector<complex_int>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;
	if (N == 2)
		return {x[0] + x[1], x[0] - x[1]};

	std::vector<complex_int> y(N / 2);
	std::vector<complex_int> z(N / 4);
	std::vector<complex_int> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		double theta1 = normalize(	  NEG2PI * k / N);
		double theta3 = normalize(3 * NEG2PI * k / N);

		complex_int a =  (x[k] + x[k + N / 2]);
		complex_int b = ((x[k] - x[k + N / 2]) - (x[k + N / 4] - x[k + 3 * N / 4]).times1i());
		complex_int c = 						 (x[k + N / 4] + x[k + 3 * N / 4]);
		complex_int d = ((x[k] - x[k + N / 2]) + (x[k + N / 4] - x[k + 3 * N / 4]).times1i());

		y[k] 		 = a;
		z[k] 		 = lift<Q>(b, theta1);
		y[k + N / 4] = c;
		h[k] 		 = lift<Q>(d, theta3);
	}

	auto Y = IntFFT_splitRadix_DIF<Q>(y);
	auto Z = IntFFT_splitRadix_DIF<Q>(z);
	auto H = IntFFT_splitRadix_DIF<Q>(h);

	std::vector<complex_int> X(N);
	for (int k = 0; k < N / 4; ++k) {
		X[4 * k    ] = Y[2 * k];
		X[4 * k + 1] = Z[k];
		X[4 * k + 2] = Y[2 * k + 1];
		X[4 * k + 3] = H[k];
	}

	return X;
}
