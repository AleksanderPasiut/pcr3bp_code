///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include <lyra/core3d.hpp>

#include <capd_utils/capd/basic_types.hpp>
#include <capd_utils/capd/solution_curve.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief convert array<float> to VectorT
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename VectorT, size_t dimension>
VectorT convert(const std::array<double, dimension>& input)
{
    VectorT ret( dimension );
    for (size_t i = 0; i < dimension; ++i)
    {
        ret[i] = input[i];
    }
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief convert RVector to array<float>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t dimension>
std::array<float, dimension> convert(const CapdUtils::RVector& vector)
{
    std::array<float, dimension> ret;
    for (size_t i = 0; i < dimension; ++i)
    {
        ret[i] = vector[i];
    }
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief convert IVector to array<float>
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t dimension>
std::array<float, dimension> convert(const CapdUtils::IVector& vector)
{
    std::array<float, dimension> ret;
    for (size_t i = 0; i < dimension; ++i)
    {
        ret[i] = vector[i].mid().leftBound();
    }
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief convert IVector to RulerSet
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t dimension>
Leo::RulerSet<dimension> convert_to_ruler_set(
    const CapdUtils::IVector& vector,
    size_t count,
    size_t subcount)
{
    Leo::RulerSet<dimension> ret;
    for (size_t i = 0; i < dimension; ++i)
    {
        ret[i] = Leo::Ruler<double>(
            vector[i].leftBound(),
            vector[i].rightBound(),
            count,
            subcount);
    }
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render solution curve object with continuous line
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t image_dimension>
class CapdSolutionCurveRenderable : public Lyra::ManifoldInterface<1, image_dimension>
{
public:
    CapdSolutionCurveRenderable(
        const CapdUtils::SolutionCurve<MapT>& curve,
        double minimum,
        double maximum,
        size_t points,
        Leo::Color color)
            : Lyra::ManifoldInterface<1, image_dimension>({
                Leo::RulerSet<1>({
                    Leo::Ruler<>(minimum, maximum, points, 1)
                }),
                color
                })
            , m_curve(curve)
    {}

    void fill(float thickness)
    {
        auto func = [this](const std::array<double, 1>& in) -> std::array<float, image_dimension>
        {
            return convert<image_dimension>( m_curve( in[0] ) );
        };

        Lyra::ManifoldInterface<1, image_dimension>::fill(func, thickness);
    }

private:
    const CapdUtils::SolutionCurve<MapT>& m_curve;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render solution curve object with points
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t image_dimension>
class CapdSolutionCurvePointRenderable : public Lyra::ManifoldInterface<0, image_dimension>
{
public:
    CapdSolutionCurvePointRenderable(
        Lyra::CoreObjects<image_dimension>& objects,
        const CapdUtils::SolutionCurve<MapT>& curve,
        double minimum,
        double maximum,
        size_t points,
        Leo::Color color)
            : Lyra::ManifoldInterface<0, image_dimension>(objects, { points, color })
            , m_curve(curve)
            , m_ruler(minimum, maximum, points, 1)
    {}

    void fill(float thickness)
    {
        auto func = [this](const size_t& in) -> std::array<float, image_dimension>
        {
            const double t = m_ruler.minimum() + in * m_ruler.step();
            return convert<image_dimension>( m_curve( t ) );
        };

        Lyra::ManifoldInterface<0, image_dimension>::fill(func, thickness);
    }

private:
    const CapdUtils::SolutionCurve<MapT>& m_curve;

    Leo::Ruler<> m_ruler;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render solution curve object with continuous line and with
//! coordinate change
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t image_dimension>
class CapdSolutionCurveWithCchngRenderable : public Lyra::ManifoldInterface<1, image_dimension>
{
public:
    CapdSolutionCurveWithCchngRenderable(
        const CapdUtils::SolutionCurve<MapT>& curve,
        MapT cchng,
        double minimum,
        double maximum,
        size_t points,
        Leo::Color color)
            : Lyra::ManifoldInterface<1, image_dimension>({
                Leo::RulerSet<1>({
                    Leo::Ruler<>(minimum, maximum, points, 1)
                }),
                color
                })
            , m_cchng(cchng)
            , m_curve(curve)
    {}

    void fill(float thickness)
    {
        auto func = [this](const std::array<double, 1>& in) -> std::array<float, image_dimension>
        {
            return convert<image_dimension>( m_cchng( m_curve( in[0] ) ) );
        };

        Lyra::ManifoldInterface<1, image_dimension>::fill(func, thickness);
    }

private:
    MapT m_cchng;
    const CapdUtils::SolutionCurve<MapT>& m_curve;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render N->M map with wireframe
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t domain_dimension, size_t image_dimension>
class CapdMapRenderable : public Lyra::ManifoldInterface<domain_dimension, image_dimension>
{
public:
    CapdMapRenderable(
        Lyra::CoreObjects<image_dimension>& objects,
        MapT& map,
        const Leo::RulerSet<domain_dimension>& ruler_set,
        Leo::Color color)
            : Lyra::ManifoldInterface<domain_dimension, image_dimension>(objects, { ruler_set, color })
            , m_map(map)
    {}

    void fill(float thickness)
    {
        auto func = [this](const std::array<double, domain_dimension>& in) -> std::array<float, image_dimension>
        {
            using VectorType = typename MapT::VectorType;
            const VectorType vector = convert<VectorType, domain_dimension>(in);
            return convert<image_dimension>( m_map( vector ) );
        };

        Lyra::ManifoldInterface<domain_dimension, image_dimension>::fill(func, thickness);
    }

private:
    MapT& m_map;
};

#ifdef DISABLED

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render N->4 map with wireframe
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t domain_dimension>
class CapdMapRenderable<MapT, domain_dimension, 4> : public Lyra::Manifold4<domain_dimension>
{
public:
    static constexpr size_t image_dimension = 4;
    
    CapdMapRenderable(
        MapT& map,
        const Leo::RulerSet<domain_dimension>& ruler_set,
        Leo::Matrix4f const& rotation)
            : Lyra::Manifold4<domain_dimension>( ruler_set, rotation )
            , m_map(map)
    {}

    void fill(float thickness)
    {
        auto func = [this](const std::array<double, domain_dimension>& in) -> std::array<float, image_dimension>
        {
            using VectorType = typename MapT::VectorType;
            const VectorType vector = convert<VectorType, domain_dimension>(in);
            return convert<image_dimension>( m_map( vector ) );
        };

        Lyra::Manifold4<domain_dimension>::fill(func, thickness);
    }

private:
    MapT& m_map;
};

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render map with discrete domain as points (in 2. or 3. dim)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t image_dimension>
class CapdMapPointRenderable : public Lyra::ManifoldInterface<0, image_dimension>
{
public:
    CapdMapPointRenderable(
        MapT& map,
        size_t points,
        Leo::Color color)
            : Lyra::ManifoldInterface<0, image_dimension>({ points, color })
            , m_map(map)
    {}

    void fill(float thickness)
    {
        auto func = [this](const size_t& in) -> std::array<float, image_dimension>
        {
            return convert<image_dimension>( m_map( in ) );
        };

        Lyra::ManifoldInterface<0, image_dimension>::fill(func, thickness);
    }

private:
    MapT& m_map;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render single vector as point (in 2. or 3. dim)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename VectorType, size_t image_dimension>
class CapdVectorRenderable : public Lyra::ManifoldInterface<0, image_dimension>
{
public:
    CapdVectorRenderable(Lyra::CoreObjects<image_dimension>& objects, VectorType vector, Leo::Color color)
        : Lyra::ManifoldInterface<0, image_dimension>(objects, { 1, color })
        , m_vector(vector)
    {}

    void fill(float thickness)
    {
        auto func = [this](const size_t& in) -> std::array<float, image_dimension>
        {
            return convert<image_dimension>( m_vector );
        };

        Lyra::ManifoldInterface<0, image_dimension>::fill(func, thickness);
    };

    const VectorType & get_vector() const noexcept
    {
        return m_vector;
    }

private:
    const VectorType m_vector;
};

#ifdef DISABLED

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render single vector as point (4. dim)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename VectorType>
class CapdVectorRenderable4 : public Lyra::Manifold4<0>
{
public:
    CapdVectorRenderable4(VectorType vector, Leo::Matrix4f const& rotation)
        : Lyra::Manifold4<0>(Leo::RulerSet<0>( 1 ), rotation)
        , m_vector(vector)
    {}

    void fill(float thickness)
    {
        auto func = [this](const size_t& in) -> std::array<float, 4>
        {
            return convert<4>( m_vector );
        };

        Lyra::Manifold4<0>::fill(func, thickness);
    };

    const VectorType & get_vector() const noexcept
    {
        return m_vector;
    }

private:
    const VectorType m_vector;
};

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief render single interval vector with map (in 2. or 3. dim)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t domain_dimension, size_t image_dimension>
class CapdMapIntervalRenderable : public Lyra::ManifoldInterface<domain_dimension, image_dimension>
{
public:
    using VectorType = typename MapT::VectorType;

    CapdMapIntervalRenderable(
        MapT& map,
        const CapdUtils::IVector& ivector,
        size_t count,
        size_t subcount,
        Leo::Color color)
            : Lyra::ManifoldInterface<domain_dimension, image_dimension>({
                convert_to_ruler_set<domain_dimension>(ivector, count, subcount), color })
            , m_map(map)
    {}

    void fill(float thickness)
    {
        auto func = [this](const std::array<double, domain_dimension>& in) -> std::array<float, image_dimension>
        {
            using VectorType = typename MapT::VectorType;
            const VectorType vector = convert<VectorType, domain_dimension>(in);
            return convert<image_dimension>( m_map( vector ) );
        };

        Lyra::ManifoldInterface<domain_dimension, image_dimension>::fill(func, thickness);
    }

private:
    MapT& m_map;
};

}