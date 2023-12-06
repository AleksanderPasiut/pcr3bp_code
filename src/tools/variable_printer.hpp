///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <fstream>

namespace CapdUtils
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief A utility to print information about specific variable to file or file stream
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
