#include "Fourier.h"
#include <numeric>

int main() {
	Fourier f4(2);
	Fourier f8(3);
	Fourier f16(4);

	{
		unsigned int n = std::thread::hardware_concurrency();
		std::cout << n << " concurrent threads are supported.\n";

		size_t N{ 3 };
		Fourier f(N);

		std::vector<std::complex<double>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); ++j) {
			sig.at(j) = 1.0;
		}
		std::vector<std::complex<double>> res = f.Transform(sig);

		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
		}
		for (size_t j = 0; j < sig.size(); j += 4) {
			sig.at(1 + j) = 2.0;
			sig.at(3 + j) = 0.0;
		}
		res = f.Transform(sig);

		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		res = f.Transform(sig);
	}
	
	{
		size_t N{ 2 };
		Fourier f(N);

		std::vector<std::complex<double>> sig(myPow(2, N));
		auto sz = sig.size();
		for (auto j = 0; j < sig.size(); j+=2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		std::vector<std::complex<double>> res = f.Transform(sig);
	}

	std::complex<double> z(1.0, 1.0);
	double d = 5.0;
	z = z * d;
}