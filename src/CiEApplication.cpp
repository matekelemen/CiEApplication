// --- Kratos Includes ---
#include "factories/standard_linear_solver_factory.h"
#include "includes/kratos_components.h"

// --- Internal Includes ---
#include "CiEApplication/CiEApplication.hpp"
#include "CiEApplication/SYCLCGSolver.hpp"

// --- CiE Includes ---
#include "version.hpp"


namespace Kratos {


CiEApplication::CiEApplication()
    : KratosApplication("CiEApplication")
{}


void CiEApplication::Register() {
    std::cout << "CiE version " << cie::Version::local << std::endl;

    #ifdef CIE_ENABLE_SYCL
    {
        static auto factory = StandardLinearSolverFactory<
            TUblasSparseSpace<double>,
            TUblasDenseSpace<double>,
            SYCLCGSolver<
                TUblasSparseSpace<double>,
                TUblasDenseSpace<double>>>();
        KratosComponents<LinearSolverFactory<
            TUblasSparseSpace<double>,
            TUblasDenseSpace<double>
        >>::Add("cg-sycl", factory);
    }

    {
        static auto factory = StandardLinearSolverFactory<
            TUblasSparseSpace<float>,
            TUblasDenseSpace<double>,
            SYCLCGSolver<
                TUblasSparseSpace<float>,
                TUblasDenseSpace<double>>>();
        KratosComponents<LinearSolverFactory<
            TUblasSparseSpace<float>,
            TUblasDenseSpace<double>
        >>::Add("cg-sycl", factory);
    }
    #endif
}


std::string CiEApplication::Info() const {
    return "CiEApplication";
}


void CiEApplication::PrintInfo(std::ostream& rStream) const {
    rStream << cie::Version::local;
}


void CiEApplication::PrintData(std::ostream&) const {}


} // namespace Kratos
