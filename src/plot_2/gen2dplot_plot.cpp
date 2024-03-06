///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_2d.hpp>
#include <lyra/core2d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

// #include "gen2dplot_plot_x.hpp"
// #include "gen2dplot_plot_x0.hpp"
#include "gen2dplot_plot_x1.hpp"

namespace Pcr3bpProof
{

class CoreInteriorGen2dPlot : CoreInteriorBase
{
private:
    Lyra::Core2d& m_core_ref;

    std::unique_ptr<X> x_ptr;

public:
    CoreInteriorGen2dPlot(Lyra::Core2d& core_ref) : CoreInteriorBase(), m_core_ref(core_ref)
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
    using Core = Taurus::DefaultCore2d<Pcr3bpProof::CoreInteriorGen2dPlot>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("gen2dplot"), argc, argv);

    window.show();
    window.run();
    return 0;
}
