#pragma once

#include "../util.h"
#include <vector>

//////////// DFT ////////////

std::vector<complex> DFT(const std::vector<complex> &x) {
	int N = (int) x.size();
	std::vector<complex> X(x.size());

	for (int k = 0; k < N; ++k) {

		complex sum;
		for (int n = 0; n < N; ++n) {
			complex tf = cexp(NEG2PI * k * n / N);

			sum += x[n] * tf;
		}

		X[k] = sum;
	}

	return X;
}

//////////// Radix 2 FFT ////////////

std::vector<complex> FFT_radix2_DIT(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex> y(N / 2);
	std::vector<complex> z(N / 2);
	for (int k = 0; k < N / 2; ++k) {
		y[k] = x[2 * k	  ];
		z[k] = x[2 * k + 1];
	}

	auto Y = FFT_radix2_DIT(y);
	auto Z = FFT_radix2_DIT(z);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 2; ++k) {
		complex tf = cexp(NEG2PI * k / N);

		X[k]		 = Y[k] + tf * Z[k];
		X[k + N / 2] = Y[k] - tf * Z[k];
	}

	return X;
}

std::vector<complex> FFT_radix2_DIF(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex> y(N / 2);
	std::vector<complex> z(N / 2);
	for (int k = 0; k < N / 2; ++k) {
		complex tf = cexp(NEG2PI * k / N);

		y[k] = 		(x[k] + x[k + N / 2]);
		z[k] = tf * (x[k] - x[k + N / 2]);
	}

	auto Y = FFT_radix2_DIF(y);
	auto Z = FFT_radix2_DIF(z);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 2; ++k) {
		X[2 * k]	 = Y[k];
		X[2 * k + 1] = Z[k];
	}

	return X;
}

//////////// Radix 4 FFT ////////////

std::vector<complex> FFT_radix4_DIT(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex> y(N / 4);
	std::vector<complex> z(N / 4);
	std::vector<complex> g(N / 4);
	std::vector<complex> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		y[k] = x[4 * k	  ];
		z[k] = x[4 * k + 1];
		g[k] = x[4 * k + 2];
		h[k] = x[4 * k + 3];
	}

	auto Y = FFT_radix4_DIT(y);
	auto Z = FFT_radix4_DIT(z);
	auto G = FFT_radix4_DIT(g);
	auto H = FFT_radix4_DIT(h);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 4; ++k) {
		complex tf1 = cexp(    NEG2PI * k / N);
		complex tf2 = cexp(2 * NEG2PI * k / N);
		complex tf3 = cexp(3 * NEG2PI * k / N);

		X[k] 			 = (Y[k] + tf2 * G[k]) + 	  (tf1 * Z[k] + tf3 * H[k]);
		X[k + N / 4]	 = (Y[k] - tf2 * G[k]) - 1i * (tf1 * Z[k] - tf3 * H[k]);
		X[k + N / 2]	 = (Y[k] + tf2 * G[k]) - 	  (tf1 * Z[k] + tf3 * H[k]);
		X[k + 3 * N / 4] = (Y[k] - tf2 * G[k]) + 1i * (tf1 * Z[k] - tf3 * H[k]);
	}

	return X;
}

std::vector<complex> FFT_radix4_DIF(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;

	std::vector<complex> y(N / 4);
	std::vector<complex> z(N / 4);
	std::vector<complex> g(N / 4);
	std::vector<complex> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		complex tf1 = cexp(    NEG2PI * k / N);
		complex tf2 = cexp(2 * NEG2PI * k / N);
		complex tf3 = cexp(3 * NEG2PI * k / N);

		y[k] = 		  (x[k] + x[k + N / 2]) + 	   (x[k + N / 4] + x[k + 3 * N / 4]);
		z[k] = tf1 * ((x[k] - x[k + N / 2]) - 1i * (x[k + N / 4] - x[k + 3 * N / 4]));
		g[k] = tf2 * ((x[k] + x[k + N / 2]) - 	   (x[k + N / 4] + x[k + 3 * N / 4]));
		h[k] = tf3 * ((x[k] - x[k + N / 2]) + 1i * (x[k + N / 4] - x[k + 3 * N / 4]));
	}

	auto Y = FFT_radix4_DIF(y);
	auto Z = FFT_radix4_DIF(z);
	auto G = FFT_radix4_DIF(g);
	auto H = FFT_radix4_DIF(h);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 4; ++k) {
		X[4 * k	   ] = Y[k];
		X[4 * k + 1] = Z[k];
		X[4 * k + 2] = G[k];
		X[4 * k + 3] = H[k];
	}

	return X;
}

//////////// Split radix FFT ////////////

std::vector<complex> FFT_splitRadix_DIT(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;
	if (N == 2)
		return {x[0] + x[1], x[0] - x[1]};

	std::vector<complex> y(N / 2);
	std::vector<complex> z(N / 4);
	std::vector<complex> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		y[2 * k] 	 = x[4 * k	  ];
		z[k] 		 = x[4 * k + 1];
		y[2 * k + 1] = x[4 * k + 2];
		h[k] 		 = x[4 * k + 3];
	}

	auto Y = FFT_splitRadix_DIT(y);
	auto Z = FFT_splitRadix_DIT(z);
	auto H = FFT_splitRadix_DIT(h);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 4; ++k) {
		complex tf1 = cexp(    NEG2PI * k / N);
		complex tf3 = cexp(3 * NEG2PI * k / N);

		X[k] 			 = Y[k] 		+ 	   (tf1 * Z[k] + tf3 * H[k]);
		X[k + N / 4] 	 = Y[k + N / 4] - 1i * (tf1 * Z[k] - tf3 * H[k]);
		X[k + N / 2] 	 = Y[k] 		- 	   (tf1 * Z[k] + tf3 * H[k]);
		X[k + 3 * N / 4] = Y[k + N / 4] + 1i * (tf1 * Z[k] - tf3 * H[k]);
	}

	return X;
}

std::vector<complex> FFT_splitRadix_DIF(const std::vector<complex>& x) {
	int N = (int) x.size();
	if (N == 1)
		return x;
	if (N == 2)
		return {x[0] + x[1], x[0] - x[1]};

	std::vector<complex> y(N / 2);
	std::vector<complex> z(N / 4);
	std::vector<complex> h(N / 4);
	for (int k = 0; k < N / 4; ++k) {
		complex tf1 = cexp(    NEG2PI * k / N);
		complex tf3 = cexp(3 * NEG2PI * k / N);

		y[k] 		 = 		  (x[k] + x[k + N / 2]);
		z[k] 		 = tf1 * ((x[k] - x[k + N / 2]) - 1i * (x[k + N / 4] - x[k + 3 * N / 4]));
		y[k + N / 4] = 		  							   (x[k + N / 4] + x[k + 3 * N / 4]);
		h[k] 		 = tf3 * ((x[k] - x[k + N / 2]) + 1i * (x[k + N / 4] - x[k + 3 * N / 4]));
	}

	auto Y = FFT_splitRadix_DIT(y);
	auto Z = FFT_splitRadix_DIT(z);
	auto H = FFT_splitRadix_DIT(h);

	std::vector<complex> X(N);
	for (int k = 0; k < N / 4; ++k) {
		X[4 * k    ] = Y[2 * k];
		X[4 * k + 1] = Z[k];
		X[4 * k + 2] = Y[2 * k + 1];
		X[4 * k + 3] = H[k];
	}

	return X;
}
