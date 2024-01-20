///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include "tools/auxiliary_functions.hpp"

#include "covering_relations_test_base.hpp"
#include "covering_relation_checker.hpp"
#include "parallelogram_covering_checker.hpp"

#include "scaled_local_poincare4_map.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class CoveringRelationsTest_ParallelogramCoveringDerivativeCheck : public CoveringRelationsTestBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    CoveringRelationsTest_ParallelogramCoveringDerivativeCheck(const CoveringRelationsSetup& setup)
        : CoveringRelationsTestBase<MapT>(setup)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check parallelogram coverings around fixed point
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void parallelogram_covering_derivative_check()
    {
        const ScalarType L = this->m_basic_objects.m_parallelogram_coverings_parameters.L;

        MapT eta = AuxiliaryFunctions<MapT>::eta( L );
        MapT eta_inverse = AuxiliaryFunctions<MapT>::eta( -L );

        std::list<MatrixType> der_list {};
        for (int i = 0; i < 4; ++i)
        {
            const int first = i;
            const int second = (i+1) % 4;
            
            ScaledLocalPoincare4_Map<MapT> poincare
            {
                this->m_basic_objects.m_vf_reg_pos2,
                this->m_basic_objects.m_hamiltonian_reg2,
                this->m_basic_objects.m_order,
                this->m_periodic_orbit_coordsys.at(first),
                this->m_periodic_orbit_coordsys.at(second),
                this->m_gain_factor
            };

            CapdUtils::CompositeMap<MapT, MapT&, decltype(poincare)&, MapT&> aligned_poincare
            {
                std::ref(eta),
                std::ref(poincare),
                std::ref(eta_inverse)
            };

            MatrixType der(2,2);
            aligned_poincare(N, der);
            print_var(der);

            der_list.emplace_back(der);
        }

        auto it = der_list.begin();

        MatrixType der_union = *it;
        for (++it; it != der_list.end(); ++it)
        {
            der_union = capd::vectalg::intervalHull(*it, der_union);
        }

        print_var(der_union);

        ParallelogramCoveringChecker<MapT> parallelogram_covering_checker( der_union );
    }
};

}
