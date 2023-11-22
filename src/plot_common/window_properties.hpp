///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>

namespace Ursa
{

inline Serpent::SglHostWindowBase::Properties create_window_properties(const std::string title)
{
    Serpent::GLFW::get();
    Serpent::SglHostWindowBase::Properties properties;
    properties.width = 1024;
    properties.height = 768;
    properties.title = "Ursa - " + title;
    properties.maximized = false;
    properties.event_timeout = 0.0;

    return properties;
}

}
