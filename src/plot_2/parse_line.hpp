///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <array>
#include <string>

template<size_t N>
inline bool parse_line(std::array<double, N>& out, std::istream& istr)
{
    auto it = out.begin();
    for (size_t n = 1; n < N; ++n)
    {
        std::string line;
        std::getline(istr, line, ' ');

        if (istr.good() == false)
        {
            return false;
        }

        *(it++) = std::stod(line);
    }

    std::string line;
    std::getline(istr, line);

    if (istr.good() == false)
    {
        return false;
    }

    *(it++) = std::stod(line);

    return istr.good();
}
