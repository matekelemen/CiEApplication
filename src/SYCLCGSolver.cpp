#ifdef CIE_ENABLE_SYCL

// --- Kratos Includes ---
#include "spaces/ublas_space.h"

// --- CiE Includes ---
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/solvers/inc/SYCLCSROperator.hpp"
#include "packages/solvers/inc/DiagonalOperator.hpp"
#include "packages/solvers/inc/ConjugateGradients.hpp"
#include "packages/utilities/inc/reorder.hpp"

// --- Internal Includes ---
#include "CiEApplication/SYCLCGSolver.hpp"

// --- STL Includes ---
#include <optional>


namespace Kratos {


template <class TS, class TD>
struct SYCLCGSolver<TS,TD>::Impl {
    using Value = typename TS::DataType;

    using Index = typename TS::IndexType;

    using Space = cie::linalg::SYCLSpace<Value>;

    using IndexSpace = cie::linalg::SYCLSpace<Index>;

    using Operator = cie::linalg::LinearOperator<Space>;

    using Solver =cie::linalg::ConjugateGradients<Space>;

    std::shared_ptr<Space> mpSpace;

    std::shared_ptr<IndexSpace> mpIndexSpace;

    std::string mPreconditionerName;

    std::size_t mMaxIterations;

    Value mAbsoluteTolerance;

    Value mRelativeTolerance;

    int mVerbosity;

    typename Space::Vector mEntries, mSolution, mRhs, mPreconditionerData;

    typename IndexSpace::Vector mRowExtents, mColumnIndices;

    std::shared_ptr<Operator> mpLhsOperator, mpPreconditioner;

    std::shared_ptr<Solver> mpSolver;

    std::string mReordering;

    std::vector<Index> mOrdering;
}; // struct SYCLCGSolver::Impl


template <class TS, class TD>
SYCLCGSolver<TS,TD>::SYCLCGSolver()
    : SYCLCGSolver(Parameters())
{}


template <class TS, class TD>
SYCLCGSolver<TS,TD>::SYCLCGSolver(Parameters Settings)
    : mpImpl(new Impl) {
        KRATOS_TRY
            Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

            // Construct a SYCL queue.
            const std::string device_name = Settings["device"].GetString();
            std::shared_ptr<sycl::queue> p_queue;
            if (device_name.empty()) {
                p_queue = std::make_shared<sycl::queue>(sycl::default_selector_v);
            } else {
                std::optional<sycl::device> maybe_device;
                for (const sycl::device& candidate : sycl::device::get_devices()) {
                    const std::string candidate_name = candidate.get_info<sycl::info::device::name>();
                    if (candidate_name.find(device_name) != std::string::npos) {
                        maybe_device = candidate;
                        break;
                    }
                }
                KRATOS_ERROR_IF_NOT(maybe_device.has_value())
                    << "cannot find SYCL device matching " << device_name;
                p_queue = std::make_shared<sycl::queue>(maybe_device.value());
            }

            // Construct linalg spaces.
            mpImpl->mpSpace = std::make_shared<typename Impl::Space>(p_queue);
            mpImpl->mpIndexSpace = std::make_shared<typename Impl::IndexSpace>(p_queue);

            // Parse settings.
            mpImpl->mMaxIterations = Settings["max_iterations"].GetInt();
            mpImpl->mAbsoluteTolerance = Settings["absolute_tolerance"].GetDouble();
            mpImpl->mRelativeTolerance = Settings["relative_tolerance"].GetDouble();
            mpImpl->mVerbosity = Settings["verbosity"].GetInt();
            mpImpl->mPreconditionerName = Settings["preconditioner"].GetString();
            mpImpl->mReordering = Settings["reordering"].GetString();
        KRATOS_CATCH("")
}


template <class TS, class TD>
SYCLCGSolver<TS,TD>::SYCLCGSolver(SYCLCGSolver&& rRhs) = default;


template <class TS, class TD>
SYCLCGSolver<TS,TD>::~SYCLCGSolver() = default;


template <class TS, class TD>
void SYCLCGSolver<TS,TD>::InitializeSolutionStep(
    typename TS::MatrixType& rLhs,
    typename TS::VectorType&,
    typename TS::VectorType&) {
        KRATOS_TRY
            const std::size_t system_size = TS::Size1(rLhs);

            // Compute and perform reordering if requested.
            if (mpImpl->mReordering == "cuthill-mckee") {
                mpImpl->mOrdering.resize(system_size);
                cie::makeReordering<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rLhs.index1_data()[0], rLhs.index1_data().size()},
                    {&rLhs.index2_data()[0], rLhs.index2_data().size()},
                    {&rLhs.value_data()[0], rLhs.value_data().size()},
                    cie::ReorderingStrategy::CuthillMcKee);
            } else if (mpImpl->mReordering == "reverse-cuthill-mckee") {
                mpImpl->mOrdering.resize(system_size);
                cie::makeReordering<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rLhs.index1_data()[0], rLhs.index1_data().size()},
                    {&rLhs.index2_data()[0], rLhs.index2_data().size()},
                    {&rLhs.value_data()[0], rLhs.value_data().size()},
                    cie::ReorderingStrategy::CuthillMcKee);
            } else if (mpImpl->mReordering != "none") {
                KRATOS_ERROR << "unknown reordering \"" << mpImpl->mReordering << "\n";
            }

            if (!mpImpl->mOrdering.empty()) {
                cie::reorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rLhs.index1_data()[0], rLhs.index1_data().size()},
                    {&rLhs.index2_data()[0], rLhs.index2_data().size()},
                    {&rLhs.value_data()[0], rLhs.value_data().size()});
            }

            // Copy the matrix to the device.
            mpImpl->mRowExtents = mpImpl->mpIndexSpace->makeVector(rLhs.index1_data().size());
            mpImpl->mColumnIndices = mpImpl->mpIndexSpace->makeVector(rLhs.index2_data().size());
            mpImpl->mEntries = mpImpl->mpSpace->makeVector(rLhs.value_data().size());

            mpImpl->mpIndexSpace->assign(
                mpImpl->mRowExtents,
                {&rLhs.index1_data()[0], rLhs.index1_data().size()});
            mpImpl->mpIndexSpace->assign(
                mpImpl->mColumnIndices,
                {&rLhs.index2_data()[0], rLhs.index2_data().size()});
            mpImpl->mpSpace->assign(
                mpImpl->mEntries,
                {&rLhs.value_data()[0], rLhs.value_data().size()});

            const cie::linalg::CSRView<
                const typename Impl::Value,
                const typename Impl::Index
            > lhs(
                TS::Size2(rLhs),
                {mpImpl->mRowExtents.get(), mpImpl->mRowExtents.size()},
                {mpImpl->mColumnIndices.get(), mpImpl->mColumnIndices.size()},
                {mpImpl->mEntries.get(), mpImpl->mEntries.size()});

            // Construct the LHS operator.
            mpImpl->mpLhsOperator = std::make_shared<cie::linalg::SYCLCSROperator<
                typename Impl::Index,
                typename Impl::Value>
            >(lhs, mpImpl->mpSpace);

            // Allocate space for the RHS and solution vectors.
            mpImpl->mSolution = mpImpl->mpSpace->makeVector(system_size);
            mpImpl->mRhs = mpImpl->mpSpace->makeVector(system_size);

            // Construct a preconditioner.
            if (mpImpl->mPreconditionerName == "diagonal") {
                mpImpl->mpPreconditioner = std::make_shared<cie::linalg::DiagonalOperator<typename Impl::Space>>(
                    cie::linalg::makeDiagonalOperator<
                        typename Impl::Value,
                        typename Impl::Index,
                        typename Impl::Value>(
                            lhs, mpImpl->mpSpace));
            } else if (mpImpl->mPreconditionerName != "none") {
                KRATOS_ERROR
                << "unknown preconditioner \""
                << mpImpl->mPreconditionerName
                << "\"";
            }

            // Construct the solver.
            typename Impl::Solver::Statistics settings {
                .iterationCount = mpImpl->mMaxIterations,
                .absoluteResidual = mpImpl->mAbsoluteTolerance,
                .relativeResidual = mpImpl->mRelativeTolerance};
            mpImpl->mpSolver = std::make_shared<typename Impl::Solver>(
                mpImpl->mpLhsOperator,
                mpImpl->mpSpace,
                mpImpl->mpPreconditioner,
                settings,
                mpImpl->mVerbosity);
        KRATOS_CATCH("")
}


template <class TS, class TD>
bool SYCLCGSolver<TS,TD>::PerformSolutionStep(
    typename TS::MatrixType& rLhs,
    typename TS::VectorType& rSolution,
    typename TS::VectorType& rRhs) {
        KRATOS_TRY
            if (!mpImpl->mOrdering.empty()) {
                cie::reorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rRhs[0], rRhs.size()});
                cie::reorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rSolution[0], rSolution.size()});
            }

            // Transfer the initial solution and rhs to the device.
            mpImpl->mpSpace->assign(
                mpImpl->mSolution,
                {&rSolution[0], rSolution.size()});
            mpImpl->mpSpace->assign(
                mpImpl->mRhs,
                {&rRhs[0], rRhs.size()});

            // Invoke the solver.
            mpImpl->mpSolver->product(0, mpImpl->mRhs, 1, mpImpl->mSolution);

            // Fetch the solution from the device.
            mpImpl->mpSpace->assign(
                {&rSolution[0], rSolution.size()},
                mpImpl->mSolution);

            if (!mpImpl->mOrdering.empty()) {
                cie::reverseReorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rLhs.index1_data()[0], rLhs.index1_data().size()},
                    {&rLhs.index2_data()[0], rLhs.index2_data().size()},
                    {&rLhs.value_data()[0], rLhs.value_data().size()});
                cie::reverseReorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rRhs[0], rRhs.size()});
                cie::reverseReorder<typename Impl::Index,typename Impl::Value>(
                    mpImpl->mOrdering,
                    {&rSolution[0], rSolution.size()});
            }

            // Check for convergence.
            const typename Impl::Solver::Statistics stats = mpImpl->mpSolver->getStats().value();
            return stats.absoluteResidual < mpImpl->mAbsoluteTolerance
                || stats.relativeResidual < mpImpl->mRelativeTolerance;
        KRATOS_CATCH("")
}


template <class TS, class TD>
void SYCLCGSolver<TS,TD>::FinalizeSolutionStep(
    typename TS::MatrixType&,
    typename TS::VectorType&,
    typename TS::VectorType&) {
        KRATOS_TRY
            std::unique_ptr<Impl> pImpl(new Impl);
            pImpl->mpSpace = mpImpl->mpSpace;
            pImpl->mpIndexSpace = mpImpl->mpIndexSpace;
            pImpl->mPreconditionerName = mpImpl->mPreconditionerName;
            pImpl->mMaxIterations = mpImpl->mMaxIterations;
            pImpl->mAbsoluteTolerance = mpImpl->mAbsoluteTolerance;
            pImpl->mRelativeTolerance = mpImpl->mRelativeTolerance;
            pImpl->mVerbosity = mpImpl->mVerbosity;
            pImpl->mReordering = mpImpl->mReordering;
            mpImpl = std::move(pImpl);
        KRATOS_CATCH("")
}


template <class TS, class TD>
Parameters SYCLCGSolver<TS,TD>::GetDefaultParameters() const {
    return Parameters(R"({
    "solver_type" : "cg-sycl",
    "device" : "",
    "max_iterations" : 1e3,
    "absolute_tolerance" : 1e-6,
    "relative_tolerance" : 1e-6,
    "preconditioner" : "none",
    "reordering" : "none",
    "verbosity" : 1
})");
}


template class SYCLCGSolver<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class SYCLCGSolver<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos


#endif
