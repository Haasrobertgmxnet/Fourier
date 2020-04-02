#pragma once

//all math stuff
#include <cmath>
#include <complex>

//std containers
#include <vector>
#include <deque>
#include <memory>

//multithreading and async tasks
#include <future>
#include <thread>

#ifdef DEBUG
//Debug output
#include <iostream>
#include <mutex>
extern std::mutex coutMutex;//for std::cout from thread
#define TESTCOMMAND(x) x
#else
//no output
#define TESTCOMMAND(x) 
#endif

template<typename T> constexpr T tPi = 3.14159265358979323846264338328L;

//float or double
using decimal = float;

//Pi
constexpr decimal Pi = tPi<decimal>;

//single-threaded or multi-threaded
enum class Mode : unsigned int { SingleThreaded, MultiThreaded };

static size_t myPow(size_t x, size_t p)
{
	if (p == 0) return 1;
	if (p == 1) return x;
	size_t tmp = myPow(x, p / 2);
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}

static size_t myLog2(size_t x)
{
	size_t res{ 0 };
	size_t tmp{ x };
	while (tmp = tmp >> 1) {
		res++;
	}
	return res;
}

class Fourier
{
	struct MultiThreading {
		MultiThreading(Fourier* _parentObj, unsigned int _numberOfThreads) :
			parentObj{ _parentObj },
			numberOfThreads{ _numberOfThreads }
		{
			numberOfCPUCores = std::thread::hardware_concurrency();
			numberOfCPUCores = (numberOfCPUCores > 0) ? numberOfCPUCores : guessOfCPUCores;
			numberOfElemsPerThread = parentObj->numberOfElems / numberOfThreads;
		}

		void setNumberOfThreads(unsigned int _numberOfThreads) {
			numberOfThreads = _numberOfThreads;
		}

		Fourier* parentObj{ nullptr };
		unsigned int guessOfCPUCores{ 4 };
		unsigned int numberOfCPUCores{ guessOfCPUCores };
		unsigned int numberOfThreads{ 2 };
		size_t numberOfElemsPerThread{ 2 };

		struct WorkingPackData {
			WorkingPackData() = default;
			WorkingPackData(std::complex<decimal> _currentUnitRoot,
				std::complex<decimal> _sigValue1,
				std::complex<decimal> _sigValue2,
				size_t _indexSigValue1) :
				currentUnitRoot{ _currentUnitRoot },
				sigValue1{ _sigValue1 },
				sigValue2{ _sigValue2 },
				indexSigValue1{ _indexSigValue1 }
			{}
			std::complex<decimal> currentUnitRoot;
			std::complex<decimal> sigValue1;
			std::complex<decimal> sigValue2;
			size_t indexSigValue1;
		};

		struct WorkingPack {
			WorkingPack(std::vector<WorkingPackData> _workingData)
			{
				workingData.reset(new std::vector<WorkingPackData>(_workingData));
			}
			bool operator()() {
				std::for_each(workingData.get()->begin(),
					workingData.get()->end(), [this](WorkingPackData& w)
					{
						std::complex<decimal> u{ w.sigValue1 };
						std::complex<decimal> v{ w.sigValue2 * w.currentUnitRoot };
						w.sigValue1 = u + v;
						w.sigValue2 = u - v;
					});
				TESTCOMMAND(std::lock_guard<std::mutex> guard(coutMutex); std::cout << "Thread ID: " << std::this_thread::get_id() << std::endl;)
					return true;
			}
			std::shared_ptr<std::vector<WorkingPackData>> workingData;
		};
	};

public:
	Fourier(size_t _binExp, Mode _calcMode, unsigned int _numberOfThreads = 2) :binExp{ _binExp },
		numberOfElems{ myPow(2,binExp) },
		calcMode{ _calcMode },
		sigma{ std::vector<size_t>(numberOfElems) },
		rootOfUnity{ std::vector<std::complex<decimal>>(numberOfElems) }
	{
		if (calcMode == Mode::SingleThreaded) {
			Transformation = [this](std::vector<std::complex<decimal>> v) {return Transform<Mode::SingleThreaded>(v); };
		}
		if (calcMode == Mode::MultiThreaded) {
			multiThreadedFourier.reset(new MultiThreading(this, _numberOfThreads));
			Transformation = [this](std::vector<std::complex<decimal>> v) {return Transform<Mode::MultiThreaded>(v); };
		}
		BitReversal();
		AssignRootOfUnity();
	}

	template<Mode mod>
	std::vector<std::complex<decimal>> Transform(std::vector<std::complex<decimal>>);

	std::function<std::vector<std::complex<decimal>>(std::vector<std::complex<decimal>>)> Transformation;

private:
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
		std::complex<decimal> zeta{ std::exp(std::complex<decimal>(0, -2 * Pi / numberOfElems)) };
		std::complex<decimal> zpow{ zeta };
		rootOfUnity.at(0) = std::complex<decimal>(1.0, 0.0);
		for (size_t k = 1; k < numberOfElems; ++k) {
			rootOfUnity.at(sigma.at(k)) = zpow;
			zpow = zpow * zeta;
		}
	}

	size_t binExp{ 2 };
	size_t numberOfElems{ 4 };

	Mode calcMode;
	std::vector<size_t> sigma;
	std::vector<std::complex<decimal>> rootOfUnity;
	bool bitReversalDone{ false };
	std::shared_ptr<Fourier::MultiThreading> multiThreadedFourier{ nullptr };

	template<>
	std::vector<std::complex<decimal>> Transform< Mode::MultiThreaded >(std::vector<std::complex<decimal>> _signal) {
		if (multiThreadedFourier == nullptr) {
			return std::vector<std::complex<decimal>>(0);
		}

		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<decimal>>(0);
		}

		//some abbreviations
		auto numberOfThreads{ multiThreadedFourier->numberOfThreads };
		auto numberOfElemsPerThread{ multiThreadedFourier->numberOfElemsPerThread };

		//necessary adjustments: 
		//if there would be more threads
		//than working packages
		if (2 * numberOfThreads > numberOfElems) {
			numberOfThreads = numberOfElems / 2;
		}

		//necessary adjustments: 
		//if there would be too few threads
		if (numberOfElems>myPow(2,22)) {
			//I try this:
			//numberOfThreads = myPow(2, myLog2(numberOfElems) / 2);
			//numberOfThreads = myPow(2, myLog2(numberOfElems)- 20);
			numberOfThreads = 8;
		}

		auto _tmp{ std::vector<std::complex<decimal>>(numberOfElems) };
		size_t M{ numberOfElems / 2 };
		size_t K{ 2 };
		std::vector<MultiThreading::WorkingPackData> workingPackData;
		workingPackData.reserve(numberOfElemsPerThread);

		for (size_t gen = 1; gen <= binExp; ++gen) {
			size_t l{ 0 };
			size_t m{ 1 };

			std::vector<MultiThreading::WorkingPack> workingPacks;
			workingPacks.reserve(numberOfThreads);
			std::deque<std::packaged_task<bool()>> workingTasks;
			std::vector<std::future<bool>> workingResults;
			workingResults.reserve(numberOfThreads);

			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; ++k) {
					MultiThreading::WorkingPackData currentData(rootOfUnity.at(j), _signal.at(k), _signal.at(M + k), k);
					workingPackData.emplace_back(std::move(currentData));

					if (2 * m == numberOfElemsPerThread) {
						MultiThreading::WorkingPack currentPack(workingPackData);
						workingPacks.emplace_back(currentPack);
						//std::packaged_task<bool()> currentTask(std::move(currentPack));
						std::packaged_task<bool()> currentTask(currentPack);
						workingTasks.emplace_back(std::move(currentTask));
						workingResults.emplace_back(workingTasks.back().get_future());
						workingPackData.clear();
						m = 0;
					}
					m++;
				}
				l = l + 2 * M;
				if (l >= numberOfElems) {
					break;
				}
			}

			while (not workingTasks.empty())
			{
				std::packaged_task<bool()> currentTask = std::move(workingTasks.front());
				workingTasks.pop_front();
				std::thread fftThread(std::move(currentTask));
				fftThread.detach();
			}
			for (auto& currentResult : workingResults) {
				bool bl{ currentResult.get() };
			}

			for (auto currentPack : workingPacks) {
				for (size_t j = 0; j < currentPack.workingData.get()->size(); ++j) {
					auto currentData{ currentPack.workingData.get()->at(j) };
					size_t index1{ currentData.indexSigValue1 };
					_signal.at(index1) = currentData.sigValue1;
					_signal.at(index1 + M) = currentData.sigValue2;
				}
			}

			M = M / 2;
			K = 2 * K;
		}

		size_t j{ 0 };
		std::for_each(_tmp.begin(), _tmp.end(),
			[this, _signal, &j](std::complex<decimal>& w) {
				w = _signal.at(sigma[j++]);
			});
		//for (size_t j = 0; j < numberOfElems; j++) {
		//	_tmp.at(j) = _signal.at(sigma.at(j));
		//}

		return _tmp;
	}

	template<>
	std::vector<std::complex<decimal>> Transform< Mode::SingleThreaded >(std::vector<std::complex<decimal>> _signal) {
		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<decimal>>(0);
		}
		auto _tmp{ std::vector<std::complex<decimal>>(numberOfElems) };
		size_t M{ numberOfElems / 2 };
		size_t K{ 2 };

		for (size_t gen = 1; gen <= binExp; gen++) {
			size_t l = 0;
			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; k++) {
					std::complex<decimal> u = _signal.at(k);
					TESTCOMMAND(std::cout << "(gen,j,k,l, M): " << "("<< gen <<", " << j << ", " << k << ", " << l << ", " << M << ")" << std::endl;)
					std::complex<decimal> v = _signal.at(M + k) * rootOfUnity.at(j);
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

