project(CiEApplication CXX)

set(application_core_name "CiECore")
set(application_name "${PROJECT_NAME}")

# Collect dependencies.
find_package(linalg REQUIRED)

# Define main library.
add_library(${application_core_name} SHARED)
target_link_libraries(
    ${application_core_name}
    PUBLIC
        KratosCore
        cie::linalg)
set_target_properties(
    ${application_core_name}
    PROPERTIES
        COMPILE_DEFINITIONS
            "CiEApplication=EXPORT_API")

# Common headers.
target_include_directories(
    ${application_core_name}
    PRIVATE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${KRATOS_SOURCE_DIR}/kratos>"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${KRATOS_SOURCE_DIR}/kratos>"
        "$<INSTALL_INTERFACE:include>")

# Common sources.
file(
    GLOB
    common_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
target_sources(
    ${application_core_name}
    PRIVATE
    ${common_sources})

# Python bindings.
file(
    GLOB
    binding_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.cpp")
pybind11_add_module(
    Kratos${application_name}
    MODULE THIN_LTO
    ${binding_sources})
target_link_libraries(
    Kratos${application_name}
    PRIVATE
        ${application_core_name})
set_target_properties(
    Kratos${application_name}
    PROPERTIES
        PREFIX
        "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(Kratos${application_name} PROPERTIES SUFFIX .pyd) # Suffix on Windows
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(Kratos${application_name} PROPERTIES SUFFIX .so) # Suffix on MacOS
endif()

kratos_python_install_directory(
    ${INSTALL_PYTHON_USING_LINKS}
    "${CMAKE_CURRENT_SOURCE_DIR}/src/scripts"
    KratosMultiphysics/${application_name})

# Install targets
install(TARGETS ${application_core_name} DESTINATION libs)
install(TARGETS Kratos${application_name} DESTINATION libs)

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};${application_core_name}" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};${application_name}" PARENT_SCOPE)
