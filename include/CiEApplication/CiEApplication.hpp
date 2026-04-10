#pragma once

// --- Kratos Includes ---
#include "includes/kratos_application.h"

// --- STL Includes ---
#include <ostream>


namespace Kratos {


class KRATOS_API(CiEApplication) CiEApplication : public KratosApplication {
public:
    KRATOS_CLASS_POINTER_DEFINITION(CiEApplication);

    CiEApplication();

    void Register() override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rStream) const override;

    void PrintData(std::ostream& rStream) const override;
}; // class CiEApplication


} // namespace Kratos
