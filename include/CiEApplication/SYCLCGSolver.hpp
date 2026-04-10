#pragma once
#ifdef CIE_ENABLE_SYCL

// --- Kratos Includes ---
#include "linear_solvers/linear_solver.h"
#include "includes/kratos_parameters.h"

// --- STL Includes ---
#include <memory>


namespace Kratos {


template <class TSparse, class TDense>
class SYCLCGSolver : public LinearSolver<TSparse,TDense> {
public:
    SYCLCGSolver();

    SYCLCGSolver(Parameters Settings);

    SYCLCGSolver(SYCLCGSolver&& rRhs);

    ~SYCLCGSolver();

    void InitializeSolutionStep(
        typename TSparse::MatrixType& rLhs,
        typename TSparse::VectorType& rSolution,
        typename TSparse::VectorType& rRhs) override;

    bool PerformSolutionStep(
        typename TSparse::MatrixType& rLhs,
        typename TSparse::VectorType& rSolution,
        typename TSparse::VectorType& rRhs) override;

    void FinalizeSolutionStep(
        typename TSparse::MatrixType& rLhs,
        typename TSparse::VectorType& rSolution,
        typename TSparse::VectorType& rRhs) override;

    virtual Parameters GetDefaultParameters() const;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class SYCLCGSolver


} // namespace Kratos


#endif
