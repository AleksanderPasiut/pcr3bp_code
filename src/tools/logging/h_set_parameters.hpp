///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <array>
#include <list>

namespace CapdUtils
{

enum class HsetType
{
    Argument,
    Image,
    LeftImage,
    RightImage
};

class HsetParameters
{
public:
    HsetType type;
    std::array<double, 4> coordsys_origin;
    std::array<double, 4> coordinates;

    void serialize(std::ostream& ostr) const
    {
        serialize_type(ostr, type);
        serialize_arr(ostr, coordsys_origin);
        serialize_arr(ostr, coordinates);
        ostr << '\n';
    }

    void deserialize(std::istream& istr)
    {
        type = deserialize_type(istr);
        deserialize_arr(istr, coordsys_origin);
        deserialize_arr(istr, coordinates);

        if (istr.get() != '\n')
        {
            throw std::logic_error("Unexpected new line delimiter!");
        }
    }

private:
    static void serialize_arr(std::ostream& ostr, const std::array<double, 4>& arr)
    {
        for (double v : arr)
        {
            ostr << v << ';';
        }
    }

    static void deserialize_arr(std::istream& istr, std::array<double, 4>& arr)
    {
        for (double& v : arr)
        {
            istr >> v;
            if (!istr.good())
            {
                throw std::logic_error("Parsing failed!");
            }

            if (istr.get() != ';')
            {
                throw std::logic_error("Unexpected delimiter!");
            }
        }
    }

    static void serialize_type(std::ostream& ostr, HsetType type)
    {
        switch (type)
        {
            case HsetType::Argument: ostr << "Argument;"; break;
            case HsetType::Image: ostr << "Image;"; break;
            case HsetType::LeftImage: ostr << "LeftImage;"; break;
            case HsetType::RightImage: ostr << "RightImage;"; break;
            default: throw std::logic_error("Unknown h-set type!");
        }
    }

    static HsetType deserialize_type(std::istream& istr)
    {
        std::string type_str {};
        std::getline(istr, type_str, ';');
        if (type_str == "Argument")
        {
            return HsetType::Argument;
        }
        if (type_str == "Image")
        {
            return HsetType::Image;
        }
        if (type_str == "LeftImage")
        {
            return HsetType::LeftImage;
        }
        if (type_str == "RightImage")
        {
            return HsetType::RightImage;
        }

        throw std::logic_error("Unknown h-set type string!");
    }
};

inline void serialize_hset_parameters_list(std::ostream& ostr, const std::list<HsetParameters>& hset_parameters_list)
{
    for (const HsetParameters& hp : hset_parameters_list)
    {
        hp.serialize(ostr);
    }
}

inline std::list<HsetParameters> deserialize_hset_parameters_list(std::istream& istr)
{
    std::list<HsetParameters> ret {};
    for (;;)
    {
        try
        {
            HsetParameters hp {};
            hp.deserialize(istr);
            ret.push_back(hp);
        }
        catch(...)
        {
            break;
        }
    }
    return ret;
}

}
