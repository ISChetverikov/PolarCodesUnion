#include <cmath>

#include "../include/Exceptions.h"
#include "../include/GaussianApproximation.h"


GaussianApproximation::GaussianApproximation(double sigma) {
	_sigma = sigma;
}

double GaussianApproximation::phi(double x) {
	if (x > 10.0)
		return PI_HALF_SQRT * sqrt(2.0 / x) * (1.0 - 10.0 / 7.0 / x) * exp(-x / 4);

	if (x < 0)
		throw ArgumentOutOfRangeException("Argument of function phi: " + std::to_string(x) + " is less than 0");

	return exp(-0.4527*pow(x, 0.86) + 0.0218);
}

double GaussianApproximation::phiDerivative(double x) {
	if (x < 0)
		throw ArgumentOutOfRangeException("Argument of derivative of function phi: " + std::to_string(x) + " is less than 10");

	return PI_HALF_SQRT * sqrt(2) * (-0.25 - 1.0 / 7 / x + 15.0 / 7.0 / x / x) * exp(-x / 4) / sqrt(x);
}

double GaussianApproximation::NumericalSolve(double y) {

	double e = 10e-14;
	int maxIter = 50;

	double x0 = 0;
	double x1 = 10 + e;

	size_t i = 0;
	while (abs(x0 - x1) >= e && i < maxIter) {
		double d = phiDerivative(x1);
		if (d > -e)
			d = -e;

		x0 = x1;
		x1 = x1 - (phi(x1) - y) / d;
		if (x1 < 10)
			x1 = 10;
		i++;
	}

	return x1;
}

double GaussianApproximation::phiInv(double y) {
	if (y < PHI_MIN || y > PHI_MAX)
		throw ArgumentOutOfRangeException("Argument of inversion of function phi: " + std::to_string(y) + " is out of range");

	if (y > PHI_GAP_MAX)
		return pow((0.0218 - log(y)) / 0.4527, 1 / 0.86);

	if (y > PHI_GAP_MIN)
		return 10;

	return NumericalSolve(y);
}

double GaussianApproximation::GetMu(size_t i, size_t n) {
	if (n == 1)
		return 2 / _sigma / _sigma;

	if (i % 2 == 1) {
		size_t j = i / 2 + 1;
		double temp = 1 - phi(GetMu(j, n / 2));

		return phiInv(1 - temp * temp);
	}
	else {
		size_t j = i / 2;
		return 2 * GetMu(j, n / 2);
	}
}

double GaussianApproximation::GetChannelErrorProbability(size_t j, size_t n) {
	double mu = GetMu(j, n);

	return 1.0 / 2 * erfc(sqrt(mu) / 2);
}