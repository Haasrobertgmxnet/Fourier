#include <iostream>
#include <chrono>
#include <random>
#include <mutex>
std::mutex coutMutex;

//#define DEBUG 1

#include "Fourier.h"
static const unsigned int NUM = 21;

int main() {
	unsigned int n = std::thread::hardware_concurrency();
	std::cout << n << " concurrent threads are supported.\n";
	{
		size_t ex;
		std::cout << "Exponent zur Basis zwei eingeben: ";
		std::cin >> ex;
		//std::random_device seed;
		//std::mt19937 engine(seed());
		std::vector<std::complex<double>> sig(myPow(2, ex));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j+=2) {
			sig.at(j) = 1.0;
		}
		{
			auto start = std::chrono::system_clock::now();
			Fourier f(ex, Mode::SingleThreaded);
			std::vector<std::complex<double>> res = f.Transform(sig);
			std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
			std::cout << "time single-threaded: " << dur.count() << " seconds" << std::endl;
		}
		{
			auto start = std::chrono::system_clock::now();
			Fourier f(ex, Mode::MultiThreaded);
			std::vector<std::complex<double>> res = f.Transform(sig);
			std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
			std::cout << "time multi-threaded: " << dur.count() << " seconds" << std::endl;
		}
	}

	{
		size_t N{ 3 };
		Fourier f(N, Mode::MultiThreaded);

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
		Fourier f(N, Mode::MultiThreaded);

		std::vector<std::complex<double>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j+=2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		std::vector<std::complex<double>> res = f.Transform(sig);
	}

	{
		size_t N{ 2 };
		Fourier f(N, Mode::MultiThreaded);

		std::vector<std::complex<double>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		std::vector<std::complex<double>> res = f.Transform(sig);
	}

	std::cout << "Bitte Zeichen eingeben: ";
	char ch;
	std::cin >> ch;
}