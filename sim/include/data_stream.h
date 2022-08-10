#pragma once

#include <fstream>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <sstream>

#include "simulation_constants.h"

struct DataStream
{
    std::fstream file;

    // by default will append, creates files but not directories
    // if you just want to put something in a file
    DataStream(const std::string &fileName) : file(fileName, std::fstream::app)
    {
        file.precision(5);
        if (!file.is_open())
            throw std::runtime_error("Couldn't open file");
    }

    void Add_Header(int noOfIterations = _noOfIterations, int N = _N,
                          double tMin = _tMin, double tMax = _tMax, double dt = _dt, double maxStableAbundance = _maxStableAbundance,
                          int timeStepsKept = _timeStepsKept)
    {
        file << "# Integrated from " << tMin << " to " << tMax << ". dt = " << dt << ".\n"
             << "# Iterations = " << noOfIterations << ".\n"
             << "# Infinite abundance cutoff = " << maxStableAbundance << ".\n"
             << "# Number of timesteps kept = " << timeStepsKept << "\n"
             << "# N = " << N << "\n#\n";


        file << "mu,nu,sigma,rho,gamma,convergent,fixed,unique,timedOut,";

        for (int i = 0; i != _N; i++)
        {
            file << "x_" << i << (i != _N - 1 ? ',' : '\n');
        }
    }
};

template <class T>
DataStream &operator<<(DataStream &ds, const T &x)
{
    ds.file << x;
    return ds;
}
