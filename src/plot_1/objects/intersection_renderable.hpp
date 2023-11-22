///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "tools/test_tools.hpp"
#include "tools/list_reader.hpp"
#include "capd_renderable.hpp"

#include <carina/projection_map.hpp>

namespace Ursa
{

class IntersectionRenderable
{
public:
    IntersectionRenderable(Lyra::Core3d& core_ref, const std::list<RVector>& intersections)
        : m_core_ref(core_ref)
        , m_reader(intersections)
        , m_renderable(m_reader, intersections.size(), Leo::Color(0.0, 1.0, 1.0))
    {
        m_renderable.fill(1.2e-2f);
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~IntersectionRenderable() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

private:
    Lyra::Core3d& m_core_ref;

    using ListReaderBase = Ursa::ListReader<std::list<RVector>>;
    class ListReader : public ListReaderBase
    {
    public:
        ListReader(const std::list<RVector>& intersections)
            : ListReaderBase(intersections)
            , m_projection(Carina::ProjectionMap<RMap>::create(5, {2, 3, 0}))
        {}

        RVector operator() (size_t index)
        {
            return m_projection( (*this)[index] );
        }

    private:
        RMap m_projection;
    };

    ListReader m_reader;
    CapdMapPointRenderable<ListReader, 3> m_renderable;
};

}
