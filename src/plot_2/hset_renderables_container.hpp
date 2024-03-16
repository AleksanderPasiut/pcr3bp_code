///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "plot_2/objects/hset_renderable.hpp"
#include "plot_2/objects/renderable4d_with_params.hpp"

#include "load_hset_parameters_list.hpp"
#include "hset_parameters_to_coordsys_converter.hpp"

namespace Pcr3bpProof
{

class HsetRenderablesContainer
{
public:
    struct Params
    {
        Pcr3bp::RegBasicObjects<RMap>& basic_objects;
        bool show_arg_h_sets;
        bool show_img_h_sets;
        bool show_limg_h_sets;
        bool show_rimg_h_sets;
        float reg_evo_thickness;
    };

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    HsetRenderablesContainer(
        Lyra::Core3d& core_ref,
        const Manifold4_Transformation & transformation_ref,
        std::vector<Coordsys> const & periodic_orbit_coordsys_vector,
        std::vector<Coordsys> const & homoclinic_orbit_coordsys_vector
        );

    void update(Params const & params);

    void refresh();

private:

    std::list<Renderable4d_WithParams<HsetRenderable>> m_h_sets {};

    std::list<CapdUtils::HsetParameters> m_hset_parameters_list
    {
        load_hset_parameters_list()
    };

    HsetParametersToCoordsysConverter m_hset_parameter_to_coordsys_converter;
};

}
