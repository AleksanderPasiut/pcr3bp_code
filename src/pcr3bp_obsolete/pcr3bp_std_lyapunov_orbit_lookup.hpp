///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/misc/set_utilities.hpp>

#include "tools/test_tools.hpp"
#include "pcr3bp_std_lyapunov_orbit.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Lookup entry for pairs (x0, py0) for which px1 == 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Pcr3bpStdLyapunovOrbitLookupTable
{
public:
    struct Entry
    {
        Real x;
        Real p;

        bool operator< (const Entry& rhs) const
        {
            return this->x < rhs.x;
        }
    };

    static Pcr3bpStdLyapunovOrbitLookupTable& get()
    {
        static Pcr3bpStdLyapunovOrbitLookupTable instance;
        return instance;
    }

    Entry get_base(Real x0) const noexcept
    {
        auto key_func = [](const Entry& entry) -> double { return entry.x; };
        return Leo::get_by_double_key(m_lookup, x0, key_func);
    }

    std::set<Entry>::const_iterator begin() const
    {
        return m_lookup.begin();
    }

    std::set<Entry>::const_iterator end() const
    {
        return m_lookup.end();
    }

private:
    Pcr3bpStdLyapunovOrbitLookupTable()
        : m_lookup({
            { 0.003, 12.2758909089832 },
            { 0.004, 11.82492344575792 },
            { 0.005, 11.41970775251516 },
            { 0.006, 11.05298153176244 },
            { 0.007, 10.71900071325091 },
            { 0.008, 10.41315428537666 },
            { 0.009, 10.13169220178198 },
            { 0.010, 9.871529114789606 },
            { 0.015, 8.813287613838666 },
            { 0.020, 8.03088013154764 },
            { 0.025, 7.42183671125734 },
            { 0.030, 6.930138133968224 },
            { 0.035, 6.522231439997716 },
            { 0.040, 6.176630240445058 },
            { 0.050, 5.618732786402781 },
            { 0.100, 4.078777768218560 },
            { 0.150, 3.325150108090716 },
            { 0.200, 2.854259797779463 },
            { 0.250, 2.52300823057559 },
            { 0.300, 2.272973130460007 },
            { 0.350, 2.07524547115377 },
            { 0.400, 1.913681992537075 },
            { 0.450, 1.778498619029175 },
            { 0.500, 1.66342489408133 },
            { 0.550, 1.56430416551522 },
            { 0.600, 1.478350576300372 },
            { 0.650, 1.40369352634577 },
            { 0.700, 1.33882233724878 },
            { 0.750, 1.280329511811914 },
            { 0.800, 1.200888715248629 },
            { 0.840, 0.9112785784840014 },
            { 0.900, 0.564707906888744 },
            { 0.950, 0.2631719166397496 },
            { 0.960, 0.14779407600376 },
            { 0.970, -0.04107948287076307 },
            { 0.980, -0.4664863942570889 },
            { 0.985, -1.081852537974346 },
            { 0.986, -1.347797313114224 }
        })
    {}

    Pcr3bpStdLyapunovOrbitLookupTable(const Pcr3bpStdLyapunovOrbitLookupTable&) = delete;
    Pcr3bpStdLyapunovOrbitLookupTable& operator= (const Pcr3bpStdLyapunovOrbitLookupTable&) = delete;

    const std::set<Entry> m_lookup;
};

}
