///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/manifold/manifold4.hpp>

#include <functional>
#include <array>

namespace Pcr3bpProof
{

class Manifold4_Transformation
{
public:
    using Func = std::function<Lyra::Point4d(std::array<double, 4>)>;

    Manifold4_Transformation(Func func) : m_func(func), m_params_counter()
    {}

    const Func& func() const noexcept
    {
        return m_func;
    }

    void inc_params_counter() noexcept
    {
        ++m_params_counter;
    }

    size_t get_params_counter() const noexcept
    {
        return m_params_counter;
    }

private:
    Func m_func;
    size_t m_params_counter;
};

class Manifold4_Transformation_Ref
{
public:
    Manifold4_Transformation_Ref(Manifold4_Transformation const & ref) : m_ref(ref)
    {}

    bool update() noexcept
    {
        const size_t params_counter = m_ref.get_params_counter();
        if (m_params_counter != params_counter)
        {
            m_params_counter = params_counter;
            return true;
        }
        
        return false;
    }

    Manifold4_Transformation const & ref() const noexcept
    {
        return m_ref;
    }

private:
    Manifold4_Transformation const & m_ref;
    size_t m_params_counter {};
};

}
