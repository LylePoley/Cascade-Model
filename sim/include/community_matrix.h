#pragma once

#include "simulation_constants.h"

#include <boost/multi_array.hpp>

#include <random>
#include <array>
#include <numeric>
#include <cmath>

using std::pow;
using std::sqrt;

struct CommunityMatrix
{
    boost::multi_array<double, 2> data; 

    void Infinite_Cascade_Fill(std::mt19937 &generator, const double &mu_L, const double& mu_U, 
                                const double &sig_L, const double& sig_U, 
                                const double &gam)
    {
        // set the size of the matrix
        data.resize(boost::extents[_N][_N]);

        // rand() generates a mean zero variance 1 gaussian variable
        std::normal_distribution<> gaussian{0, 1};
        auto rand = [&gaussian, &generator]() { return gaussian(generator); };

        for (int a = 0; a != _N; a++)
        {
            data[a][a] = -1;
            // run through all the b's bigger than a's i.e. the upper part of the matrix
            for (int b = a + 1; b != _N; b++)
            {
                // generate two mean 0, variance 1, variables with covariance gamma (cholesky decomposition)
                double x = rand();
                double y = gam * x + sqrt(1 - gam*gam) * rand();

                // the upper and lower halves of the interaction matrix have different mean, variance and are correlated (correlation coefficient = sig_L sig_U gam)
                data[a][b] = mu_U/_N + sig_U/sqrt(_N) * x;
                data[b][a] = mu_L/_N + sig_L/sqrt(_N) * y;
            }
        }
    }

    void Infinite_Cascade_Fill_Different_Bucket_sizes(std::mt19937 &generator, const double &mu_L, const double& mu_U, 
                                const double &sig_L, const double& sig_U, 
                                const double &gam, const std::array<int, _B>& bucket_sizes)
    {
        // set the size of the matrix
        data.resize(boost::extents[_N][_N]);

        // rand() generates a mean zero variance 1 gaussian variable
        std::normal_distribution<> gaussian{0, 1};
        auto rand = [&gaussian, &generator]() { return gaussian(generator); };

        for (int i = 0; i != _N; i++)
        {
            data[i][i] = -1;
            // run through all the b's bigger than a's i.e. the upper part of the matrix
            for (int j = i + 1; j != _N; j++)
            {
                // generate two mean 0, variance 1, variables with covariance gamma (cholesky decomposition)
                double x = rand();
                double y = gam * x + sqrt(1 - gam*gam) * rand();

                // the upper and lower halves of the interaction matrix have different mean, variance and are correlated (correlation coefficient = sig_L sig_U gam)
                data[i][j] = mu_U/_N + sig_U/sqrt(_N) * x;
                data[j][i] = mu_L/_N + sig_L/sqrt(_N) * y;
            }
        }

        // now set alpha^{aa}_{ij} = 0, wasteful but this isn't the bottleneck in speed by a long way
        int pSum = 0;
        for (int a = 0; a != _B; a++)
        {
            for (int i = 0; i != bucket_sizes[a]; i++)
            {
                // don't change the diagonal
                for (int j = i + 1; j != bucket_sizes[a]; j++)
                {
                    data[pSum + i][pSum + j] = 0;
                    data[pSum + j][pSum + i] = 0;
                }
            }

            // move to the start of the next bucket
            pSum += bucket_sizes[a]; 
        }
    }

};
