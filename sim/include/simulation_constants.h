#ifndef SIMULATION_CONSTANTS_H
#define SIMULATION_CONSTANTS_H

// All constants of a run of the simulation, the convention is that they are all prefixed by an underscore and cannot be changed

#include <array>
#include <thread>
#include <string>

/*
    The number of species in the community
*/ 
constexpr int _N = 200;
/* 
    The number of sub communities, in the limit _B to infty, 
    we set this equal to 1 and change the interaction matrix instead.
*/
constexpr int _B = 2;

/*
    These are some examples of finite _B communities in which the number of species in each group is varied.
    They produce communities which are virtually the same as _B = \infty communities with the same statistics
    and were only used to test independence from n^a, which is proven analytically in the supplementary material.
*/
// constexpr std::array<int, _B> _bucket_sizes{{41, 40, 39, 38, 37, 37, 36, 35, 34, 33, 33, 32, 31, 30, 30, 29, 28, 28, 27, 26, 26, 25, 25, 24, 23, 23, 23, 23, 22, 22}};
// constexpr std::array<int, _B> _bucket_sizes{{120, 0, 0, 38, 37, 37, 36, 35, 34, 33, 33, 32, 31, 30, 30, 29, 28, 28, 27, 26, 26, 25, 25, 24, 23, 23, 23, 23, 22, 22}};

// use this to find how many threads can be supported on your system
// const int _noOfThreads = std::thread::hardware_concurrency();
/*
    each thread is run concurrently, this is just to speed up simulations.
    You can set _noOfThreads = 1 and nothing bad will happen if you don't want to do multithreading.
*/
constexpr int _noOfThreads = 4;
constexpr int _iterationLoops = 100;
constexpr int _noOfIterations = _iterationLoops * _noOfThreads;

// concurrent processes are helled in a threadList
template <class T>
using ThreadList = std::array<T, _noOfThreads>;

/* 
    Integration quantities, the step size _dt is actually a minimum step size,
    the integrator used will take larger steps if the error if doing so is small enough.
*/
constexpr double _dt = 0.1;
constexpr double _tMin = 0.0;
constexpr double _tMax = 1500.0;
constexpr int _timeStepsKept = 50;

/*
    Some measures of stability for producing Fig. S1
    The diverging abundance cutoff (i.e. the M \to \infty limit) classes a run as diverging if any abundnace
    exceeds _maxStableAbundance.
    If a single run takes longer than _timeOutTimeInMs milliseconds then the sun is classed as linearly unstable.
    This time must be made large enough so as to only distinguish points which are deep in the linearly unstable regime.
*/
constexpr double _maxStableAbundance = 10000.0;
constexpr double _timeOutTimeInMs = 60000;


// put the full (not relative) directory you want simulation data to go in here
const std::string _dataDirectory("./simData/");

#endif