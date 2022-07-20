#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

#include "simulation_constants.h"

using std::pow;

// Can most likely remove the need for both Lerp and Vector_Lerp with std::type_traits
// but i cant be bothered to do it at the moment

// returns the linear interpolation of the points (x1,y1),(x2,y2) at x
template <class T>
T Lerp(const T &x, const T &x1, const T &y1, const T &x2, const T &y2)
{
    return (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1);
}


//*/
// updates the average avg with newValue, given that avg is already an average of r values
template <class T>
T Average(const T &avg, const T &newValue, int r = 1)
{
    return (r * avg + newValue) / (r + 1);
}
template <class Vec>
Vec Vector_Average(const Vec &avg, const Vec &newValue, int r = 1)
{
    if (!(avg.size() == newValue.size()))
    {
        throw std::domain_error("Arguments to Vector_Average must have the same size.");
    }

    Vec result(avg);

    for (int i = 0; i != avg.size(); i++)
    {
        result[i] = Average(avg[i], newValue[i], r);
    }
    return result;
}

template<class T>
T Variance(T& var, T& mean, const T& newValue, int r = 1)
{
    return Average(var, pow(newValue - mean, 2), r);
}

// timers, for timing stuff
//Basic Timer for stuff, starts on construction, call read to get the time since construction
template <class storageT = float, class durationT = std::milli>
class BasicTimer
{
   public:
    //A tick is a millisecond
    using tick  = std::chrono::duration<storageT, durationT>;
    using clock = std::chrono::steady_clock;

    BasicTimer() : _startTime(clock::now()) {}
    ~BasicTimer() = default;

    void ReStart() noexcept
    {
        _startTime = clock::now();
    }

    storageT Read() const
    {
        return std::chrono::duration_cast<tick>(clock::now() - _startTime).count();
    }

    //eg will cast milliseconds to 0.001
    static constexpr storageT resolution() noexcept
    {
        return static_cast<storageT>(durationT::num) / durationT::den;
    }

   protected:
    clock::time_point _startTime;
};

//Allows you to pause the timer and not time some amount of code, then restart again
template <class storageT = float, class durationT = std::milli>
class PausableTimer : public BasicTimer<storageT, durationT>
{
   public:
    using base = BasicTimer<storageT, durationT>;

    PausableTimer() : base(), _running(true) {}
    ~PausableTimer() = default;

    void Pause()
    {
        if (_running == false)
            throw std::logic_error("Double paused a timer.");

        _auxTime = base::clock::now();
        _running = !_running;
    }

    void Play()
    {
        if (_running == true)
            throw std::logic_error("Double played a timer.");

        base::_startTime += base::clock::now() - _auxTime;
        _running = !_running;
    }

   private:
    typename base::clock::time_point _auxTime;
    bool _running;
};


#endif
