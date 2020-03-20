#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <future>
#include <thread>

static size_t myPow(size_t x, size_t p)
{
	if (p == 0) return 1;
	if (p == 1) return x;
	size_t tmp = myPow(x, p / 2);
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}

class Fourier
{
	size_t binExp{ 2 };
	size_t numberOfElems{ 4 };
	std::vector<size_t> sigma;
	std::vector<std::complex<double>> rootOfUnity;
	bool bitReversalDone{ false };

	void BitReversal() {
		for (size_t k = 0; k < numberOfElems; ++k) {
			size_t l{ k };
			size_t krev{ 0 };
			for (size_t j = 0; j < binExp; ++j) {
				krev = 2 * krev + l % 2;
				l = l / 2;
			}
			sigma.at(k) = krev;
			bitReversalDone = true;
		}
	}

	void AssignRootOfUnity() {
		if (!bitReversalDone) {
			return;
		}
		std::complex<double> zeta{ std::exp(std::complex<double>(0, -2 * M_PI / numberOfElems)) };
		std::complex<double> zpow{ zeta };
		rootOfUnity.at(0) = std::complex<double>(1.0, 0.0);
		for (size_t k = 1; k < numberOfElems; ++k) {
			rootOfUnity.at(sigma.at(k)) = zpow;
			zpow = zpow * zeta;
		}
	}

public:
	Fourier(size_t _binExp) :binExp{ _binExp }, numberOfElems{ myPow(2,binExp) } {
		sigma = std::vector<size_t>(numberOfElems);
		rootOfUnity = std::vector<std::complex<double>>(numberOfElems);
		BitReversal();
		AssignRootOfUnity();
	}

	std::vector<std::complex<double>> Transform(std::vector<std::complex<double>> _signal) {
		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<double>>(0);
		}
		auto _tmp = std::vector<std::complex<double>>(numberOfElems);
		size_t M = numberOfElems / 2;
		size_t K = 2;

		for (size_t gen = 1; gen <= binExp; gen++) {
			size_t l = 0;
			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; k++) {
					std::complex<double> u = _signal.at(k);
					std::complex<double> v = _signal.at(M + k) * rootOfUnity.at(j);
					_signal[k] = u + v;
					_signal[M + k] = u - v;
				}
				l = l + 2 * M;
				if (l >= numberOfElems) {
					break;
				}
			}
			M = M / 2;
			K = 2 * K;
		}

		for (size_t j = 0; j < numberOfElems; j++) {
			_tmp.at(j) = _signal.at(sigma.at(j));
		}

		return _tmp;
	}
};

