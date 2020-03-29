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
#define TESTANWEISUNG(x) x
#else
//no output
#define TESTANWEISUNG(x) 
#endif

//pi
const double Pi = 3.141592653589793238462643383279502884197169399375;

enum class Mode : unsigned int { SingleThreaded, MultiThreaded, MT };

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
public:

	template<Mode mod>
	std::vector<std::complex<double>> Transform(std::vector<std::complex<double>>);
	std::function<std::vector<std::complex<double>>(std::vector<std::complex<double>>)> Transformation;

	Fourier(size_t _binExp, Mode _calcMode) :binExp{ _binExp },
		numberOfElems{ myPow(2,binExp) },
		calcMode{ _calcMode },
		sigma{ std::vector<size_t>(numberOfElems) },
		rootOfUnity{ std::vector<std::complex<double>>(numberOfElems) } {
		if (calcMode == Mode::MT) {
			Transformation = [this](std::vector<std::complex<double>> v) {return TransformMT(v); };
		}
		if (calcMode == Mode::SingleThreaded) {
			Transformation = [this](std::vector<std::complex<double>> v) {return Transform<Mode::SingleThreaded>(v); };
		}
		if (calcMode == Mode::MultiThreaded) {

			Transformation = [this](std::vector<std::complex<double>> v) {return Transform<Mode::MultiThreaded>(v); };
		}		
		BitReversal();
		AssignRootOfUnity();
	}

private:
	size_t binExp{ 2 };
	size_t numberOfElems{ 4 };
	Mode calcMode;
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
		std::complex<double> zeta{ std::exp(std::complex<double>(0, -2 * Pi / numberOfElems)) };
		std::complex<double> zpow{ zeta };
		rootOfUnity.at(0) = std::complex<double>(1.0, 0.0);
		for (size_t k = 1; k < numberOfElems; ++k) {
			rootOfUnity.at(sigma.at(k)) = zpow;
			zpow = zpow * zeta;
		}
	}

	struct WorkingPackData {
		std::complex<double> currentUnitRoot;
		std::complex<double> sigValue1;
		std::complex<double> sigValue2;
		size_t indexSigValue1;
		WorkingPackData() = default;
		WorkingPackData(std::complex<double> _currentUnitRoot,
			std::complex<double> _sigValue1,
			std::complex<double> _sigValue2,
			size_t _indexSigValue1) :
			currentUnitRoot{ _currentUnitRoot },
			sigValue1{ _sigValue1 },
			sigValue2{ _sigValue2 },
			indexSigValue1{ _indexSigValue1 }
		{}
	};

	struct WorkingPack {
		std::shared_ptr<std::vector<WorkingPackData>> workingData;
		WorkingPack(std::vector<WorkingPackData> _workingData)
		{
			workingData.reset(new std::vector<WorkingPackData>(_workingData));
		}
		bool operator()() {
			auto vec{ workingData.get() };
			for (size_t j = 0; j < vec->size(); ++j) {
				auto currentData{ vec->at(j) };
				std::complex<double> u = workingData.get()->at(j).sigValue1;
				std::complex<double> v = vec->at(j).sigValue2 * currentData.currentUnitRoot;
				workingData.get()->at(j).sigValue1 = u + v;
				vec->at(j).sigValue2 = u - v;
			}
			TESTANWEISUNG(std::lock_guard<std::mutex> guard(coutMutex); std::cout << "Thread ID: " << std::this_thread::get_id() << std::endl;)
				return true;
		}
	};

	bool _WorkingPackage(std::vector<std::complex<double>>& _signal,
		size_t K, size_t M, size_t l) {
		
	}

	std::vector<std::complex<double>> TransformMT(std::vector<std::complex<double>> _signal) {
		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<double>>(0);
		}
		unsigned int guessOfCPUCores{ 4 };
		unsigned int numberOfCPUCores{ std::thread::hardware_concurrency() };
		numberOfCPUCores = (numberOfCPUCores > 0) ? numberOfCPUCores : guessOfCPUCores;
		size_t numberOfElemsPerCPUCore{ numberOfElems / numberOfCPUCores };
		auto _tmp{ std::vector<std::complex<double>>(numberOfElems) };
		size_t M{ numberOfElems / 2 };
		size_t K{ 2 };
		std::vector<WorkingPackData> workingPackData;
		workingPackData.reserve(numberOfElemsPerCPUCore);

		for (size_t gen = 1; gen <= binExp; ++gen) {
			size_t l{ 0 };
			size_t m{ 1 };

			std::vector<WorkingPack> workingPacks;
			workingPacks.reserve(numberOfCPUCores);
			std::deque<std::packaged_task<bool()>> workingTasks;
			std::vector<std::future<bool>> workingResults;
			workingResults.reserve(numberOfCPUCores);

			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; ++k) {
					WorkingPackData currentData(rootOfUnity[j], _signal[k], _signal[M + k], k);
					workingPackData.emplace_back(std::move(currentData));

					if (2 * m == numberOfElemsPerCPUCore) {
						WorkingPack currentPack(workingPackData);
						workingPacks.emplace_back(currentPack);
						workingPackData.clear();
						//workingPackData.shrink_to_fit();
						//workingPackData.reserve(numberOfElemsPerCPUCore);
						//std::packaged_task<bool()> currentTask(std::move(currentPack));
						std::packaged_task<bool()> currentTask(workingPacks.back());
						workingTasks.emplace_back(std::move(currentTask));
						workingResults.emplace_back(workingTasks.back().get_future());
						
						
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

			size_t idx{ 0 };
			for (auto currentPack : workingPacks) {
				//auto vec{ currentPack.workingData.get() };
				for (size_t j = 0; j < currentPack.workingData.get()->size(); ++j) {
					auto currentData{ (*currentPack.workingData.get())[j] };
					size_t index1{ currentData.indexSigValue1 };
					//size_t index2{ currentData.indexSigValue2 };
					TESTANWEISUNG(std::lock_guard<std::mutex> guard(coutMutex); std::cout << "(M, index1, index2): (" << M << ", "  << index1 << ", " << index1+M << ")" << std::endl;)
						idx += 1;
					//TESTANWEISUNG(std::cout << "2. (index1, index2): (" << sigma.at(j) << ", " << index2 << ")" << std::endl;)
					_signal[index1] = currentData.sigValue1;
					_signal[index1 + M] = currentData.sigValue2;
				}
			}

			M = M / 2;
			K = 2 * K;
		}

		for (size_t j = 0; j < numberOfElems; j++) {
			_tmp[j] = _signal.at(sigma[j]);
		}

		return _tmp;
	}

	template<>
	std::vector<std::complex<double>> Transform< Mode::MultiThreaded >(std::vector<std::complex<double>> _signal) {
		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<double>>(0);
		}
		unsigned int guessOfCPUCores{ 4 };
		unsigned int numberOfCPUCores{ std::thread::hardware_concurrency() };
		numberOfCPUCores = (numberOfCPUCores > 0) ? numberOfCPUCores : guessOfCPUCores;
		size_t numberOfElemsPerCPUCore{ numberOfElems / numberOfCPUCores };
		auto _tmp{ std::vector<std::complex<double>>(numberOfElems) };
		size_t M{ numberOfElems / 2 };
		size_t K{ 2 };
		std::vector<WorkingPackData> workingPackData;
		workingPackData.reserve(numberOfElemsPerCPUCore);

		for (size_t gen = 1; gen <= binExp; ++gen) {
			size_t l{ 0 };
			size_t m{ 1 };

			std::vector<WorkingPack> workingPacks;
			workingPacks.reserve(numberOfCPUCores);
			std::deque<std::packaged_task<bool()>> workingTasks;
			std::vector<std::future<bool>> workingResults;
			workingResults.reserve(numberOfCPUCores);

			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; ++k) {
					WorkingPackData currentData(rootOfUnity.at(j), _signal.at(k), _signal.at(M + k), k);
					workingPackData.emplace_back(std::move(currentData));

					if (2*m == numberOfElemsPerCPUCore) {
						WorkingPack currentPack(workingPackData);
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

		for (size_t j = 0; j < numberOfElems; j++) {
			_tmp.at(j) = _signal.at(sigma.at(j));
		}

		return _tmp;
	}

	template<>
	std::vector<std::complex<double>> Transform< Mode::SingleThreaded >(std::vector<std::complex<double>> _signal) {
		if (_signal.size() != numberOfElems) {
			return std::vector<std::complex<double>>(0);
		}
		auto _tmp{ std::vector<std::complex<double>>(numberOfElems) };
		size_t M{ numberOfElems / 2 };
		size_t K{ 2 };

		for (size_t gen = 1; gen <= binExp; gen++) {
			size_t l = 0;
			for (size_t j = 0; j < K; j += 2) {
				for (size_t k = l; k < l + M; k++) {
					std::complex<double> u = _signal.at(k);
					//TESTANWEISUNG(std::cout << "(gen,j,k,l, M): " << "("<< gen <<", " << j << ", " << k << ", " << l << ", " << M << ")" << std::endl;)
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

