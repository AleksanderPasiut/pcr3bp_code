///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdint.h>
#include <iostream>

namespace Leo
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Floating info base class
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Scalar>
class FloatingInfoBase
{
protected:
    struct Split
    {};

    union ScalarWithSplit
    {
        static_assert(sizeof(Scalar) == sizeof(Split), "Scalar and Split mismatch!");

        Scalar value;
        Split split;
    };
};

template<>
struct FloatingInfoBase<float>::Split
{
    uint32_t mantissa : 23;
    uint16_t exponent : 8;
    uint16_t sign : 1;

    static const uint16_t EXPONENT_OFFSET = 0x7f + 23;
    static const uint32_t MANTISSA_IMPLICIT_ONE = (uint32_t(1) << 23);
};

template<>
struct FloatingInfoBase<double>::Split
{
    uint64_t mantissa : 52;
    uint16_t exponent : 11;
    uint16_t sign : 1;

    static const uint16_t EXPONENT_OFFSET = 0x3ff + 52;
    static const uint64_t MANTISSA_IMPLICIT_ONE = (uint64_t(1) << 52);
};

template<>
struct FloatingInfoBase<long double>::Split
{
    uint64_t mantissa : 64;
    uint16_t exponent : 15;
    uint16_t sign : 1;

    static const uint16_t EXPONENT_OFFSET = 0x3fff + 63;
    static const uint64_t MANTISSA_IMPLICIT_ONE = 0;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Floating info
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Scalar>
class FloatingInfo : private FloatingInfoBase<Scalar>
{
public:
    using MantissaType = decltype(FloatingInfoBase<Scalar>::Split::mantissa);

    FloatingInfo(Scalar arg) : m_split()
    {
        m_split.value = arg;
    }

    bool get_sign() const noexcept
    {
        return m_split.split.sign != 0;
    }

    int16_t get_exponent() const noexcept
    {
        if (m_split.split.exponent > 0)
        {
            return int16_t(m_split.split.exponent) - FloatingInfoBase<Scalar>::Split::EXPONENT_OFFSET;
        }
        else
        {
            return 1 - int16_t(FloatingInfoBase<Scalar>::Split::EXPONENT_OFFSET);
        }
    }

    MantissaType get_mantissa() const noexcept
    {
        if (m_split.split.exponent > 0)
        {
            return m_split.split.mantissa | FloatingInfoBase<Scalar>::Split::MANTISSA_IMPLICIT_ONE;
        }
        else
        {
            return m_split.split.mantissa;
        }
    }

private:
    typename FloatingInfoBase<Scalar>::ScalarWithSplit m_split;
};

template<typename MantissaType, typename ExponentType>
inline void reduce(MantissaType& mantissa, ExponentType& exponent)
{
    if (mantissa > 0)
    {
        while ((mantissa & 1u) == 0)
        {
            mantissa >>= 1;
            ++exponent;
        }
    }
}

template<typename T>
inline std::ostream& operator<< (std::ostream& ostr, const FloatingInfo<T>& arg)
{
    auto mantissa = arg.get_mantissa();
    auto exponent = arg.get_exponent();
    reduce(mantissa, exponent);

    if (mantissa > 0)
    {
        ostr << (arg.get_sign() ? "-" : "+");
        ostr << mantissa;
        ostr << "\\cdot 2^{" << exponent << "}";
    }
    else
    {
        ostr << "0";
    }
    return ostr;
}

}
