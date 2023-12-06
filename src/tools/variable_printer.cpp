///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "variable_printer.hpp"
#include "floating_info.hpp"

#include <fstream>

#include "types.hpp"

#include <capd_utils/type_cast.hpp>
#include <capd_utils/readable_scalar.hpp>

namespace CapdUtils
{

static void print_data(std::ofstream& fs, double value)
{
    fs << "Approx. value     " << value << '\n';
    fs << "Exact hex. value  " << std::string( ReadableScalar<double>(value) ) << '\n';
    fs << "Exact frac. value " << FloatingInfo<double>(value) << '\n';
}

static void print_data(std::ofstream& fs, Interval value)
{
    fs << "Approx. left bound      " << value.leftBound() << '\n';
    fs << "Approx. right bound     " << value.rightBound() << '\n';
    fs << "Approx. span            " << span(value) << '\n';
    fs << "Exact hex. left bound   " << std::string( ReadableScalar<double>(value.leftBound()) ) << '\n';
    fs << "Exact hex. right bound  " << std::string( ReadableScalar<double>(value.rightBound()) ) << '\n';
    fs << "Exact frac. left bound  " << FloatingInfo<double>(value.leftBound()) << '\n';
    fs << "Exact frac. right bound " << FloatingInfo<double>(value.rightBound()) << '\n';
}

template<>
void VariablePrinter<RMap>::print(std::ofstream& fs, std::string description, const ScalarType& value)
{
    if (fs)
    {
        fs << "Real value\n";
        fs << "Description: " << description << '\n';
        print_data(fs, value);
    }
    else
    {
        throw std::logic_error("Failed to access file stream!");
    }
}

template<>
void VariablePrinter<IMap>::print(std::ofstream& fs, std::string description, const ScalarType& value)
{
    if (fs)
    {
        fs << "Interval value\n";
        fs << "Description: " << description << '\n';
        print_data(fs, value);
    }
    else
    {
        throw std::logic_error("Failed to access file stream!");
    }
}

template<>
void VariablePrinter<RMap>::print(std::ofstream& fs, std::string description, const VectorType& value)
{
    if (fs)
    {
        fs << "Real vector value\n";
        fs << "Description: " << description << '\n';

        for (int i = 1; i <= value.dimension(); ++i)
        {
            fs << "Component " << i << '\n';
            print_data( fs, value(i) );
        }
    }
    else
    {
        throw std::logic_error("Failed to access file stream!");
    }
}

template<>
void VariablePrinter<IMap>::print(std::ofstream& fs, std::string description, const VectorType& value)
{
    if (fs)
    {
        fs << "Interval vector value\n";
        fs << "Description: " << description << '\n';

        for (int i = 1; i <= value.dimension(); ++i)
        {
            fs << "Component " << i << '\n';
            print_data( fs, value(i) );
        }
    }
    else
    {
        throw std::logic_error("Failed to access file stream!");
    }
}

}
