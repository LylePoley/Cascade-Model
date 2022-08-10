#pragma once

#include <array>
#include <chrono>

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

