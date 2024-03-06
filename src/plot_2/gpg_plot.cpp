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

#include <aries/hrclock.hpp>
#include <aries/smartfile.hpp>

namespace Pcr3bpProof
{

class X
{
public:
    using MapT = RMap;
    using ScalarType = MapT::ScalarType;
    using VectorType = MapT::VectorType;
    using MatrixType = MapT::MatrixType;

    X(Lyra::Core2d& core_ref, const std::vector<double>& param)
        : m_core_ref(core_ref)
    {
        m_core_ref.register_surface2(&m_func);
    }

private:
    class Func final : public Lyra::Surface2
    {
    public:
        Func() : Surface2({
            Leo::RulerSet<2>({
                Leo::Ruler<>(-2.5, 2.5, 800, 1),
                Leo::Ruler<>(-2.5, 2.5, 800, 1)
            })
            })
            , m_setup()
            , m_xieta(m_setup)
        {
            if (!try_load())
            {
                generate();
            }
        }

    private:
        bool try_load()
        {
            try
            {
                Aries::SmartfileReader fs("gpg_plot_out.bin");

                auto func = [this, &fs](const std::array<double, 2>& in) -> std::array<float, 2>
                {
                    double img_x;
                    fs.read_value<double>(img_x);

                    double img_y;
                    fs.read_value<double>(img_y);

                    return { float(img_x), float(img_y) };
                };

                fill(func);

                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        void generate()
        {
            Aries::SmartfileWriter writer("gpg_plot_out.bin");

            Aries::HrClock clock;
            clock.start();

            auto func = [this, &writer](const std::array<double, 2>& in) -> std::array<float, 2>
            {
                const RVector arg = { in[0], in[1] };
                RVector img = { NAN, NAN };

                try
                {
                    if ( std::abs(arg.euclNorm() - 2.0) > 1e-6 )
                    {
                        img = m_xieta( RVector{ in[0], in[1] } );
                    }
                    else
                    {
                        img = { NAN, NAN };
                    }
                }
                catch (...)
                {
                    img = { NAN, NAN };
                }

                writer.write_value<double>(img[0]);
                writer.write_value<double>(img[1]);

                return { float(img[0]), float(img[1]) };
            };

            fill(func);

            clock.stop();
            std::cout << "Computation time: " << clock.get_s() << '\n';
        }

        Pcr3bpSetupValues<RMap> m_setup;
        Pcr3bpRegPoincarePositiveU_XiEta_CE<RMap> m_xieta;

    } m_func;

    Lyra::Core2d& m_core_ref;
};

class CoreInterior : CoreInteriorBase
{
private:
    Lyra::Core2d& m_core_ref;

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

        x_ptr = std::make_unique<X>( std::ref(m_core_ref), std::cref(paramv) );
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore2d<Pcr3bpProof::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("plot_gpg"),
        argc,
        argv,
        Leo::Color(1.0, 1.0, 1.0),
        Serpent::Grid2d::Properties{
            .ruler_x = Leo::Ruler<>(-2.5, 2.5, 11, 1),
            .ruler_y = Leo::Ruler<>(-2.5, 2.5, 11, 1),
            .grid_color = Leo::Color( 0.2, 0.2, 0.2 ),
            .grid_width = 0.003f,
            .label_x = L"ξ",
            .label_y = L"η",
            .label_properties{
                .back_color = Leo::Color(1.0, 1.0, 1.0),
                .fore_color = Leo::Color(0.2, 0.2, 0.2),
                .halfheight = 0.05f,
                .align = Serpent::Align::Center
            }
        });

    window.get_core().configure_colortray(
        Leo::Color::Profile2d::ComplexDomainColoring,
        Leo::Color(1.0, 1.0, 1.0),
        Leo::Color(0.0, 0.0, 0.0),
        Serpent::Colortray::Placement{ 
            .scale = 0.5,
			.position = { 0.35, -0.45 }
        },
        Leo::Ruler<double>(-2.5, 2.5, 11, 5),
        Leo::Ruler<double>(-2.5, 2.5, 11, 5));

    window.show();
    window.run();
    return 0;
}
