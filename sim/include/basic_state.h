#ifndef BASIC_STATE_VECTOR_H
#define BASIC_STATE_VECTOR_H

#include "community_matrix.h"
#include "simulation_constants.h"
#include "helper_functions.h"

#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <cmath>
using std::pow;

#include <boost/multi_array.hpp>
#include <boost/numeric/odeint.hpp>

// should be faster than a State as it only keeps track of the average value and variance of the last 100 or so timesteps
// it is unstable if any of the states x_i go over a certain value or pass a certain variance, specified in simulation constants.h
//! this is possibly all you need in general when doing this sort of thing
struct BasicState
{
    using timeSlice = std::array<double, _N>;

    bool convergent, fixedPoint, timedOut;
    timeSlice means;

    // a state to be integrated
    BasicState() : convergent(true), fixedPoint(true), timedOut(false)
    {
        means.fill(0.0);
    }

    void Infinite_Lotka_Volterra_Fill(timeSlice &IC, CommunityMatrix &A,
                                      const double &tMin = _tMin, const double &tMax = _tMax, const double &dt = _dt,
                                      const double &maxStableAbundance = _maxStableAbundance, int timeStepsKept = _timeStepsKept)
    {
        auto Lotka_Volterra = [&A, this, maxStableAbundance, dt, tMax](const timeSlice &x, timeSlice &dxdt, const double t)
        {
            for (int i = 0; i != _N; i++)
            {
                // if any of the states is blowing up (i.e going unstable) then the system is divergent
                if (x[i] > maxStableAbundance)
                {
                    convergent = false;
                    fixedPoint = false;
                }
                // if it is converging but wobbling upwards in the last 1% of the run then its not fixed
                // it can wobble down, this catches exponentially dying species who are otherwise stable
                if ((t > 0.99*tMax) && (dxdt[i] > (0.0001 / dt)) && convergent)
                {
                    fixedPoint = false;
                }
            }
            // update dynamics if stable,
            if (convergent && fixedPoint && !timedOut)
            {
                for (auto i = 0; i != _N; i++)
                {
                    dxdt[i] = x[i] * std::inner_product(x.begin(), x.end(), A.data[i].begin(), 1.0);
                }
            }
            // this is a work around not being able to stop the integration early, does the same thing
            else
            {
                for (auto i = 0; i != _N; i++)
                    dxdt[i] = 0;
            }
        };

        BasicTimer<float, std::milli> timeOutTimer;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        int averagingCount = 0;
        boost::numeric::odeint::integrate(
            Lotka_Volterra,
            IC, tMin, tMax, dt,
            [this, tMax, dt, &averagingCount, &timeOutTimer](const timeSlice &X, const double &T)
            {
                if(timeOutTimer.Read() > _timeOutTimeInMs)
                {
                    timedOut = true;
                }
                // store the mean and variance of the last (timeStepsKept) abundances
                if (T > 0.99 * tMax)
                {
                    for (int i = 0; i != _N; i++)
                    {
                        means[i] = Average(means[i], X[i], averagingCount);
                    }
                    averagingCount++;
                }
            });
    }

    double& operator[](int i)
    {
        return means[i];
    }

    void Reset() noexcept
    {
        convergent = true;
        fixedPoint = true;
        timedOut = false;
        means.fill(0.0);
    }

    bool Is_Fixed_Point() const noexcept
    {
        return convergent && fixedPoint;
    }
};

std::ostream &operator<<(std::ostream &os, const BasicState &state)
{
    for (int j = 0; j != _N; j++)
    {
        os << (state.means[j] < 1e-10 ? 0.0 : state.means[j]) << (j < _N - 1 ? "," : "");
    }
    return os;
}

#endif