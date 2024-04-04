/*
 *  utils.hpp
 *  Header file for util functions
 */

#pragma once
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <cmath>
using namespace std;
static const unsigned int global_seed = 35;
static const double sqrt1_2 = sqrt(0.5);

inline double norm255(double v)
{
    double normV = max(0.0, v);
    normV = min(normV, 255.0);
    return normV;
}

inline double normScope(double v, double *vscope)
{
    double normV = max(vscope[0], v);
    normV = min(normV, vscope[1]);
    return normV;
}

inline double normHue(const double &h)
{
    double res = fmod(h, 360.0);
    if (res < 0)
        res += 360.0;
    return res;
}

static unsigned int seed = global_seed;
inline void resetSeed()
{
    seed = global_seed;
}

inline double randomOne()
{
    double res = sin(seed++) * 10000;
    return res - floor(res);
}

inline int randomInt(int min, int max)
{
    return floor(randomOne() * (max - min + 0.999)) + min;
}

inline double randomDouble(double min, double max)
{
    return randomOne() * (max - min) + min;
}

inline double randomDoubleSym(double r)
{
    return randomDouble(-r, r);
}

struct ColorHash
{
    std::size_t operator()(const array<int, 3> &key) const
    {
        std::size_t hash = 0;
        for (int i : key)
        {
            hash ^= std::hash<int>()(i);
        }
        return hash;
    }
};

struct ColorEqual
{
    bool operator()(const array<int, 3> &lhs, const array<int, 3> &rhs) const
    {
        return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
    }
};

struct PosHash
{
    std::size_t operator()(const array<int, 2> &key) const
    {
        std::size_t hash = 0;
        for (int i : key)
        {
            hash ^= std::hash<int>()(i);
        }
        return hash;
    }
};

struct PosEqual
{
    bool operator()(const array<int, 2> &lhs, const array<int, 2> &rhs) const
    {
        return lhs[0] == rhs[0] && lhs[1] == rhs[1];
    }
};

template <typename T>
std::vector<int> argsort(const std::vector<T> &array)
{
    const int array_len(array.size());
    std::vector<int> array_index(array_len, 0);
    for (int i = 0; i < array_len; ++i)
        array_index[i] = i;

    std::sort(array_index.begin(), array_index.end(),
              [&array](int pos1, int pos2)
              { return (array[pos1] < array[pos2]); });

    return array_index;
}

inline double floor2(double x)
{
    // floor to 2 decimal places
    return floor(x * 100.0) / 100.0;
}