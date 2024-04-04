/*
 *  matrix.hpp
 *  Header file for matrix save and load
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void writeBinaryMatrix(const std::vector<std::vector<float>> &matrix, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary | std::ios::out);
    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }
    auto N = matrix.size();
    std::vector<float> flatMatrix(N * N);
    for (std::size_t i = 0; i < N; ++i)
    {
        std::memcpy(&flatMatrix[i * N], matrix[i].data(), N * sizeof(float));
    }
    file.write(reinterpret_cast<char *>(flatMatrix.data()), N * N * sizeof(float));
    file.close();
}

void loadBinaryMatrix(std::vector<std::vector<float>> &matrix, const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary | std::ios::in);
    if (!file.is_open())
    {
        std::cerr << "Error opening file for reading." << std::endl;
        return;
    }

    auto N = matrix.size();

    std::vector<float> flatMatrix(N * N);
    file.read(reinterpret_cast<char *>(flatMatrix.data()), flatMatrix.size() * sizeof(float));
    file.close();

    for (std::size_t i = 0; i < N; ++i)
    {
        std::memcpy(matrix[i].data(), &flatMatrix[i * N], N * sizeof(float));
    }
}

void writeBinarySymmetricMatrix(const std::vector<std::vector<double>> &matrix, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary | std::ios::out);
    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    std::size_t N = matrix.size();

    std::vector<float> flatMatrix(N * (N + 1) / 2);
    std::size_t index = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            flatMatrix[index++] = matrix[i][j];
        }
    }

    file.write(reinterpret_cast<const char *>(flatMatrix.data()), flatMatrix.size() * sizeof(float));
    file.close();
}

void loadBinarySymmetricMatrix(std::vector<std::vector<double>> &matrix, const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary | std::ios::in);
    if (!file.is_open())
    {
        std::cerr << "Error opening file for reading." << std::endl;
        return;
    }

    std::size_t N = matrix.size();
    std::vector<float> flatMatrix(N * (N + 1) / 2);
    file.read(reinterpret_cast<char *>(flatMatrix.data()), flatMatrix.size() * sizeof(float));
    file.close();

    std::size_t index = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            matrix[i][j] = flatMatrix[index++];
            matrix[j][i] = matrix[i][j];
        }
    }
    return;
}
