/*
    This is another example, code similar to that used to get Fig. S1 in the appendix
    copy and paste it into main.cpp to run
*/

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

    // where the data goes, by default it adds data to the end of an existing file
    DataStream file("/home/lyle/Dropbox/PhD/Code/Cascade-Model/sim/data/finaltest.csv");
    // comment this out if you want to add data to an existing file
    file.Add_Header();

    std::mt19937 generator{std::random_device{}()};

    // these generate pairs of nu and rho, we choose them randomly so that data can be added at any time
    std::uniform_real_distribution<> nuDist(-5.0, 5.0);
    std::uniform_real_distribution<> rhoDist(std::log(0.1), std::log(10));

    ThreadList<CommunityMatrix> alpha;
    ThreadList<BasicState::timeSlice> IC; // initial condition
    ThreadList<BasicState> x;
    ThreadList<std::thread> threads;

    // average of the randomised states
    BasicState xAvg;

    // set the model parameters here, nu and rho are randomised
    double mu = 0.5, nu = nuDist(generator), 
            sigma = 0.8, rho = std::exp(rhoDist(generator)), gamma = 0.8;
    /*
        outer loop can be used to change model parameters,
        inner loop generates interaciton matrices, runs the dynamics and averages the results
    */
    for (int count = 0; count != 100; count += 1)
    {
        // iterate parameters to be changed, eg
        nu = nuDist(generator);
        rho = std::exp(rhoDist(generator));


        // change parameters so they can be fed into the functions below
        double  mu_L = mu + nu, 
                mu_U = mu - nu,
                sigma_L = sigma * rho, 
                sigma_U = sigma / rho;

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
                threads[i] = std::thread(Generate_Quenched_Disorder, std::ref(alpha[0]), std::ref(IC[i]));
            }
            for (auto &&t : threads)
            {
                t.join();
            }
            randomFillingTimer.Pause();

            integrationTimer.Play();
            for (int i = 0; i != _noOfThreads; i++)
            {
                threads[i] = std::thread(Integrate_Dynamics, std::ref(alpha[0]), std::ref(IC[i]), std::ref(x[i]));
            }
            for (auto &&t : threads)
            {
                t.join();
            }
            integrationTimer.Pause();


            // if all three are true then the data gets processed
            int unique = 0;
            int convergent = 0;
            int fixedPoint = 0;
            int timedOut = 0;

            // first check they all converge to a fixed point
            // non-convergent and non fixed point solutions dont get integrated fully
            // so these solutions have no valid data
            for (auto &&state : x)
            {
                convergent += (int)state.convergent;
                fixedPoint += (int)state.fixedPoint;
                timedOut += (int)state.timedOut;
            }

            // if theres data to look at, check that each abundance (x[thread][i]) ends up within 0.001 of the other three
            for (int i = 0; i != _N; i++)
            {
                for (int thread = 0; thread != _noOfThreads - 1; thread++)
                {
                    unique += (int)(std::abs(x[thread][i] - x[thread + 1][i]) < std::abs(1e-10 + x[thread][i]*0.001));
                }
            }

            // read out parameter and stability information
            file << mu << ',' << nu << ',' << sigma << ',' << rho << ',' << gamma << ','
                 << double(convergent)/_noOfThreads << ',' << double(fixedPoint)/_noOfThreads << ','
                 << double(unique)/((_noOfThreads - 1)*_N) << ',' << double(timedOut)/_noOfThreads << ',';

            // read out the abundances of each species
            for (int i = 0; i != _N; i++)
            {
                x[0][i] = (x[0][i] + x[1][i] + x[2][i] + x[3][i]) / 4.0;
            }
            file << x[0] << '\n';

            for (auto &&state : x)
            {
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
