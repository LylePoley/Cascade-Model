#include <boost/multi_array.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric> // for inner_product
#include <random>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <mutex>

#include "include/community_matrix.h"
#include "include/data_stream.h"
#include "include/basic_state.h"
#include "include/helper_functions.h"

#include "include/simulation_constants.h"

/*
    This is an example of a simple run of the dynamics with \nu varied, from it you could
    get data like that used for Fig. 1 in the main text.
*/
int main()
{
    /*
        Timers for run time reporting, just for fun.
        They time the total simulation runtime, as well as the time spent constructing
        interaction matrices and integrating the dynamics. Unsurprisingly, integrating
        the dynamics takes up the majority of simulation time.
    */
    BasicTimer<float, std::ratio<1, 1>> totalTimer;
    BasicTimer<float, std::milli> loopTimer;
    PausableTimer<float, std::ratio<1, 1>> integrationTimer, randomFillingTimer;
    integrationTimer.Pause();
    randomFillingTimer.Pause();

    // where the data goes
    DataStream file("./simData/finaltest.csv");
    // write out the data somewhere
    file.Add_Header("munu");

    std::mt19937 generator{std::random_device{}()};

    // std::uniform_real_distribution<> nuDist(-250.0, 250.0);
    // std::uniform_real_distribution<> rhoDist(0.0, std::log(1000));

    ThreadList<CommunityMatrix> alpha;
    ThreadList<BasicState::timeSlice> IC; // initial condition
    ThreadList<BasicState> x;
    ThreadList<std::thread> threads;

    // average of the randomised states
    BasicState xAvg;

    // set the model parameters here
    double mu = -1.0, nu = 0.0, sigma = 1.0, rho = 1.0, gamma = -1.0;
    /*
        outer loop can be used to change model parameters,
        inner loop generates interaciton matrices, runs the dynamics and averages the results
    */
    for (int count = 0; count != 1; count += 1)
    {
        // iterate parameters to be changed, eg
        nu += 0.1;
        // change parameters so they can be fed into the functions below
        double mu_L = mu + nu, mu_U = mu - nu,
               sigma_L = sigma * rho, sigma_U = sigma / rho;

        // read out parameters if they change
        std::cout << "count = " << count << ", mu = " << mu << ", nu = " << nu << ", sigma = " << sigma << ", rho = " << rho << ", gamma = " << gamma << '\n';

        for (int it = 0; it != _iterationLoops; it++)
        {
            // Feed this same function into multiple threads
            // mutex is so the random number generator is only accessed by one thread at a time
            std::mutex mut;

            // generates interaction matrices
            auto Generate_Quenched_Disorder = [&generator, &mut, mu_L, mu_U, sigma_U, sigma_L, gamma](CommunityMatrix &alpha, BasicState::timeSlice &initCondition)
            {
                // lock this function to a single thread so the generator isnt accessed simultaneously
                std::lock_guard<std::mutex> lg(mut);

                // generate the interaction matrix
                alpha.Infinite_Cascade_Fill(generator, mu_L, mu_U, sigma_L, sigma_U, gamma);
                // randomly fill the initial condition with values in (0,1)
                std::uniform_real_distribution<> unif(0.0, 1.0);
                for (auto &&ic : initCondition)
                {
                    ic = unif(generator);
                }
            };

            auto Integrate_Dynamics = [](CommunityMatrix &alpha, BasicState::timeSlice &initCondition, BasicState &state)
            {
                state.Infinite_Lotka_Volterra_Fill(initCondition, alpha);
            };

            randomFillingTimer.Play();
            for (int i = 0; i != _noOfThreads; i++)
            {
                threads[i] = std::thread(Generate_Quenched_Disorder, std::ref(alpha[i]), std::ref(IC[i]));
            }
            for (auto &&t : threads)
            {
                t.join();
            }
            randomFillingTimer.Pause();

            integrationTimer.Play();
            for (int i = 0; i != _noOfThreads; i++)
            {
                threads[i] = std::thread(Integrate_Dynamics, std::ref(alpha[i]), std::ref(IC[i]), std::ref(x[i]));
            }
            for (auto &&t : threads)
            {
                t.join();
            }
            integrationTimer.Pause();


            for (auto &&state : x)
            {
                file << mu << ',' << nu << ',' << sigma << ',' << rho << ',' << gamma << ','
                     << (int)state.convergent << ',' << (int)state.fixedPoint << ',' << 1.0 << ',' << (int)state.timedOut << ',';

                // if everything works then read the data out
                if (state.convergent)
                {
                    file << state << '\n';
                }
                else
                {
                    for (int i = 0; i != _N; i++)
                    {
                        file << 0.0 << (i != _N - 1 ? ',' : '\n');
                    }
                }
                state.Reset();
            }

            // how many iterations done so far out of the total
            std::cout << "Iteration " << _noOfThreads * it << " of "
                      << std::abs(_noOfIterations) << " finished after " << loopTimer.Read() << "ms" << std::endl;
            loopTimer.ReStart();
        }
    }

    std::cout << "Took " << totalTimer.Read() << "s" << std::endl;
    std::cout << "Integration took " << integrationTimer.Read() << "s" << '\n';
    std::cout << "Random filling took " << randomFillingTimer.Read() << "s" << '\n';

    return 0;
}
