// --- External Includes ---
#include <pybind11/pybind11.h>

// --- Internal Includes ---
#include "CiEApplication/CiEApplication.hpp"


namespace Kratos::Python {


PYBIND11_MODULE(KratosCiEApplication, rModule) {
    pybind11::class_<
        CiEApplication,
        CiEApplication::Pointer,
        KratosApplication>(rModule, "CiEApplication")
            .def(pybind11::init<>())
            ;
} // PYBIND11_MODULE


} // namespace Kratos::Python
