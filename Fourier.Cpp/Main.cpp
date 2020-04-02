#include <iostream>
#include <chrono>
#include <random>
#include <mutex>
std::mutex coutMutex;

//#define DEBUG 1

#include "Fourier.h"

int main() {
	unsigned int n = std::thread::hardware_concurrency();
	{
		size_t exNoOfPoints;
		size_t exNoOfThreads;
		std::cout << "Anzahl der Punkte: Exponent zur Basis zwei eingeben: ";
		std::cin >> exNoOfPoints;
		std::cout << "Es wird ein 1-0-Signal der Laenge 2 hoch "<< exNoOfPoints <<" erzeugt." << std::endl;

		std::vector<std::complex<decimal>> sig(myPow(2, exNoOfPoints));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j+=2) {
			sig.at(j) = 1.0;
		}
		std::cout << "Anzahl der Threads: Exponent zur Basis zwei eingeben: ";
		std::cin >> exNoOfThreads;
		std::cout << "Es wird mit 2 hoch " << exNoOfThreads << " Threads gerechnet." << std::endl;
		auto start = std::chrono::system_clock::now();
		{
			auto start = std::chrono::system_clock::now();
			Fourier f(exNoOfPoints, Mode::SingleThreaded);
			std::vector<std::complex<decimal>> res = f.Transformation(sig);
			std::chrono::duration<decimal> dur = std::chrono::system_clock::now() - start;
			std::cout << "Zeit fuer single-threaded Algorithmus: " << dur.count() << " seconds" << std::endl;
		}
		{
			auto start = std::chrono::system_clock::now();
			Fourier f(exNoOfPoints, Mode::MultiThreaded, myPow(2, exNoOfThreads));
			std::vector<std::complex<decimal>> res = f.Transformation(sig);
			std::chrono::duration<decimal> dur = std::chrono::system_clock::now() - start;
			std::cout << "Zeit fuer multi-threaded Algorithmus: " << dur.count() << " seconds" << std::endl;
		}
	}

	{
		size_t N{ 3 };
		Fourier f(N, Mode::MultiThreaded);

		std::vector<std::complex<decimal>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); ++j) {
			sig.at(j) = 1.0;
		}
		std::vector<std::complex<decimal>> res = f.Transformation(sig);

		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
		}
		for (size_t j = 0; j < sig.size(); j += 4) {
			sig.at(1 + j) = 2.0;
			sig.at(3 + j) = 0.0;
		}
		res = f.Transformation(sig);

		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		res = f.Transformation(sig);
	}
	
	{
		size_t N{ 2 };
		Fourier f(N, Mode::MultiThreaded);

		std::vector<std::complex<decimal>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j+=2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		std::vector<std::complex<decimal>> res = f.Transformation(sig);
	}

	{
		size_t N{ 2 };
		Fourier f(N, Mode::MultiThreaded);

		std::vector<std::complex<decimal>> sig(myPow(2, N));
		auto sz = sig.size();
		for (size_t j = 0; j < sig.size(); j += 2) {
			sig.at(j) = 1.0;
			sig.at(1 + j) = 0.0;
		}
		std::vector<std::complex<decimal>> res = f.Transformation(sig);
	}

	std::cout << "Bitte geben Sie ein beliebiges Zeichen ein: ";
	char ch;
	std::cin >> ch;
}