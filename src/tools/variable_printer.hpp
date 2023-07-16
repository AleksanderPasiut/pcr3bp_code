///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <fstream>

namespace Carina
{

template<typename MapT>
class VariablePrinter
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static void print(std::ofstream& fs, std::string description, const ScalarType& value);
    static void print(std::ofstream& fs, std::string description, const VectorType& value);
    static void print(std::ofstream& fs, std::string description, const MatrixType& value);

    static void print(std::string filename, std::string description, const ScalarType& value)
    {
        std::ofstream fs(filename);
        print(fs, description, value);
        fs.close();
    }

    static void print(std::string filename, std::string description, const VectorType& value)
    {
        std::ofstream fs(filename);
        print(fs, description, value);
        fs.close();
    }

    static void print(std::string filename, std::string description, const MatrixType& value)
    {
        std::ofstream fs(filename);
        print(fs, description, value);
        fs.close();
    }
};



}
