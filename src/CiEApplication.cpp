// --- Internal Includes ---
#include "CiEApplication/CiEApplication.hpp"

// --- CiE Includes ---
#include "version.hpp"


namespace Kratos {


CiEApplication::CiEApplication()
    : KratosApplication("CiEApplication")
{}


void CiEApplication::Register()
{}


std::string CiEApplication::Info() const {
    return "CiEApplication";
}


void CiEApplication::PrintInfo(std::ostream& rStream) const {
    rStream << std::string(cie::Version::local);
}


void CiEApplication::PrintData(std::ostream& rStream) const {}


} // namespace Kratos
