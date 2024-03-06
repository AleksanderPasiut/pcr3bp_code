///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_2d.hpp>
#include <lyra/core2d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

#include "tools/test_tools.hpp"
#include "pcr3bp_obsolete/pcr3bp_reg_poincare_coord.hpp"

#include "tools_cr/local_coord_vector.hpp"

#include <pcr3bp_long_path/proof_set_key.hpp>
#include "pcr3bp_primary/proof/proof_local_coords.hpp"
#include "pcr3bp_primary/proof/proof_covering_vectors.hpp"

#include "capd_renderable.hpp"

namespace Pcr3bpProof
{

class Env
{
public:
    Env()
        : m_proof_local_coords()
        , m_setup()
        , m_ppv(m_setup, 20)
        , m_covering_vectors( m_ppv, m_proof_local_coords, 4e-13 )
    {}

    const ProofLocalCoords<RMap>& get_proof_local_coords() const
    {
        return m_proof_local_coords;
    }

    const ProofCoveringVectors & get_vectors() const
    {
        return m_covering_vectors;
    }

private:
    ProofLocalCoords<RMap> m_proof_local_coords;

    Pcr3bpSetupValues<IMap> m_setup;
    Pcr3bpRegPoincarePositiveU_XiEta_CE<IMap> m_ppv;

    ProofCoveringVectors m_covering_vectors;
};

class LocalCoordRMap
{
public:
    using VectorType = RVector;
    using MatrixType = RMatrix;

    LocalCoordRMap(
        const LocalCoord& coord,
        double center_x,
        double center_y,
        double magnifier)
            : m_coord(coord)
            , m_center({ center_x, center_y })
            , m_magnifier{ { pow(2, magnifier), 0.0 }, { 0.0, pow(2, magnifier) } }
    {}

    RVector operator() (const RVector& vec) const
    {
        const RVector std = m_coord.get_directions_matrix() * vec + m_coord.get_origin();
        const RVector ret = m_magnifier * (std - m_center);
        // std::cout << ret << '\n';
        return ret;
    }

    // RVector operator() (const RVector& vec, RMatrix& mat) const
    // {
    //     mat = m_coord.get_directions_matrix();
    //     return m_coord.get_directions_matrix() * vec + m_coord.get_origin();
    // }

private:
    LocalCoord m_coord;

    RVector m_center;
    RMatrix m_magnifier;
};

class LocalCoordIVectorPlot
{
public:
    LocalCoordIVectorPlot(
        Lyra::Core2d& core,
        const LocalCoordIVector& arg,
        float thickness,
        double center_x,
        double center_y,
        double magnifier,
        Leo::Color color)
            : m_core_ref(core)
            , m_map(arg.get_local_coord(), center_x, center_y, magnifier)
            , m_renderable(m_map, arg.get_ivector(), 2, 1, color)
    {
        m_renderable.fill(thickness);
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~LocalCoordIVectorPlot() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

private:
    Lyra::Core2d& m_core_ref;

    LocalCoordRMap m_map;
    CapdMapIntervalRenderable<LocalCoordRMap, 2, 2> m_renderable;
};

class LocalCoordIVectorImage_1_2_Plot
{
public:
    LocalCoordIVectorImage_1_2_Plot(
        Lyra::Core2d& core,
        const LocalCoordIVectorImage_2dim_0& arg,
        float thickness,
        double center_x,
        double center_y,
        double magnifier,
        double multiplier,
        Leo::Color color1,
        Leo::Color color2)
            : m_image(
                core,
                LocalCoordIVector(arg.get_local_coord(), arg.get_ivector()  * pow(10, multiplier)),
                thickness,
                center_x,
                center_y,
                magnifier,
                color1)
            , m_image_left(
                core,
                LocalCoordIVector(arg.get_local_coord(), arg.get_ivector_left()  * pow(10, multiplier)),
                thickness,
                center_x,
                center_y,
                magnifier,
                color2)
            , m_image_right(
                core,
                LocalCoordIVector(arg.get_local_coord(), arg.get_ivector_right()  * pow(10, multiplier)),
                thickness,
                center_x,
                center_y,
                magnifier,
                color2)
    {}

private:
    LocalCoordIVectorPlot m_image;
    LocalCoordIVectorPlot m_image_left;
    LocalCoordIVectorPlot m_image_right;
};


class X
{
public:
    using MapT = RMap;
    using ScalarType = MapT::ScalarType;
    using VectorType = MapT::VectorType;
    using MatrixType = MapT::MatrixType;

    X(Lyra::Core2d& core_ref, const std::vector<double>& param, const Env & env)
        : m_core_ref(core_ref)
    {
        const double thickness = param[0];
        const double multiplier = param[1];

        const double center_id = param[5];

        double center_x = param[2];
        double center_y = param[3];

        const double magnifier = param[4];

        std::map<ProofSetKey, bool> show {};

        for (uint16_t i = 0; i < 19; ++i)
        {
            show[ ProofSetKey('L', i) ] = ( param[6+i] > 0 );
        }

        show[ ProofSetKey('P', 0) ] = ( param[25] > 0 );
        show[ ProofSetKey('P', 1) ] = ( param[26] > 0 );

        const ProofLocalCoords<RMap>& proof_local_coords = env.get_proof_local_coords();
        const ProofCoveringVectors& vectors = env.get_vectors();

        if (center_id < 19 && center_id >= 0)
        {
            const uint16_t id = center_id;
            center_x = proof_local_coords.get(ProofSetKey('L', id)).get_origin()[0];
            center_y = proof_local_coords.get(ProofSetKey('L', id)).get_origin()[1];
        }

        if (center_id == 19)
        {
            center_x = vectors.get_argument(ProofSetKey('P', 0)).get_local_coord().get_origin()[0];
            center_y = vectors.get_argument(ProofSetKey('P', 0)).get_local_coord().get_origin()[1];
        }

        if (center_id == 20)
        {
            center_x = vectors.get_argument(ProofSetKey('P', 1)).get_local_coord().get_origin()[0];
            center_y = vectors.get_argument(ProofSetKey('P', 1)).get_local_coord().get_origin()[1];
        }

        for (const auto& s : show)
        {
            if (s.second)
            {
                ProofSetKey key = s.first;

                const LocalCoordIVector& arg = vectors.get_argument(key);

                LocalCoordIVector argument(
                    std::ref( arg.get_local_coord() ),
                    arg.get_ivector() * pow(10, multiplier) );

                m_plots.emplace_back(
                    std::ref(m_core_ref),
                    std::cref(argument),
                    thickness,
                    center_x,
                    center_y,
                    magnifier,
                    Leo::Color(0.5, 0.0, 0.0));
                
                m_plots_img.emplace_back(
                    std::ref(m_core_ref),
                    std::cref(vectors.get_image(key)),
                    thickness,
                    center_x,
                    center_y,
                    magnifier,
                    multiplier,
                    Leo::Color(0.0, 0.5, 0.0),
                    Leo::Color(0.5, 0.5, 0.0));
            }
        }
    }

private:
    std::list<LocalCoordIVectorPlot> m_plots;
    std::list<LocalCoordIVectorImage_1_2_Plot> m_plots_img;

    Lyra::Core2d& m_core_ref;
};

class CoreInterior : CoreInteriorBase
{
private:
    Lyra::Core2d& m_core_ref;

    Env m_env;

    std::unique_ptr<X> x_ptr;

public:
    CoreInterior(Lyra::Core2d& core_ref) : CoreInteriorBase(), m_core_ref(core_ref)
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

        std::vector<double> paramv;
        paramv.reserve(this->PARAMSET_CAPACITY);

        for (size_t i = 0; i < paramv.capacity(); ++i)
        {
            paramv.push_back( this->get_param(i) );
        }

        x_ptr = std::make_unique<X>( std::ref(m_core_ref), std::cref(paramv), std::cref(m_env) );
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore2d<Pcr3bpProof::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("plot_coverings"),
        argc,
        argv,
        Leo::Color(1.0, 1.0, 1.0),
        Serpent::Grid2d::Properties{
            .ruler_x = Leo::Ruler<>(-2.0, 2.0, 9, 1),
            .ruler_y = Leo::Ruler<>(-2.0, 2.0, 9, 1),
            .grid_color = Leo::Color(0.1, 0.1, 0.1),
            .grid_width = 0.005f,
            .label_x = L"w1",
            .label_y = L"w2",
            .label_properties{
                .back_color = Leo::Color(1.0, 1.0, 1.0),
                .fore_color = Leo::Color(0.1, 0.1, 0.1),
                .halfheight = 0.03f,
                .align = Serpent::Align::Center
            }
        });

    window.show();
    window.run();
    return 0;
}
