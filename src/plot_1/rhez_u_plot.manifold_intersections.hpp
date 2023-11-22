///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <lyra/core3d.hpp>

#include "objects/interp.hpp"
#include "objects/unstable_ppv_manifold.hpp"
#include "objects/intersection_renderable.hpp"

namespace Ursa
{

using UnstableManifoldPos = UnstableManifold<Pcr3bpRegPoincarePositiveU<RMap>>;
using UnstableManifoldNeg = UnstableManifold<Pcr3bpRegPoincareNegativeU<RMap>>;

class ManifoldIntersections
{
public:
    ManifoldIntersections(
        Lyra::Core3d& core_ref,
        const UnstableManifoldPos& manifold_pos,
        const UnstableManifoldNeg& manifold_neg)
    {
        m_intersections = std::make_unique<IntersectionRenderable>(
            std::ref(core_ref),
            this->get_intersections(manifold_pos, manifold_neg));
    }

private:
    static std::list<RVector> get_intersections(
        const UnstableManifoldPos& manifold_pos,
        const UnstableManifoldNeg& manifold_neg)
    {
        std::list<RVector> ret;

        for (auto& interp_pos : manifold_pos.m_interp)
        {
            for (auto& interp_neg : manifold_neg.m_interp)
            {
                for (auto& ppos : interp_pos.m_ppv)
                {
                    for (auto& npos : interp_neg.m_ppv)
                    {
                        auto intersections = ppos.find_intersections(npos);

                        if (intersections.size() > 0)
                        {
                            for (auto & v : intersections)
                            {
                                auto it = std::find_if(ret.begin(), ret.end(), [&v](const RVector& u)
                                {
                                    return (u-v).euclNorm() < 1e-3;
                                });

                                if (it == ret.end())
                                {
                                    ret.push_back(v);
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout.precision(20);
        for (const RVector& v : ret)
        {
            std::cout << v << '\n';
        }

        return ret;
    }

    std::unique_ptr<IntersectionRenderable> m_intersections;
};

}
