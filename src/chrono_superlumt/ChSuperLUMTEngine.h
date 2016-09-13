// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni, Radu Serban
// =============================================================================
// Interfacing to the SuperLU_MT solver.
// =============================================================================

#ifndef CHSUPERLUMTENGINE_H
#define CHSUPERLUMTENGINE_H

#include "chrono_superlumt/ChApiSuperLUMT.h"
#include "chrono/core/ChSparseMatrix.h"
#include "chrono/core/ChMatrixDynamic.h"
#include <slu_mt_ddefs.h>


namespace chrono {

/// @addtogroup superlumt_module
/// @{
	enum class phase_t {
		COMPLETE = 13,
		ANALYSIS_NUMFACTORIZATION = 12,
		SOLVE = 33,
	};
/// Interface class to SuperLU_MT solver.
/// This class wraps the C interface of the solver in order to fit Chrono data structures.
/// This class can still be called by the end-user in order to solve linear systems.
/// See demo_SUPERLU_Engine for the related demo.
class ChApiSuperLUMT ChSuperLUMTEngine {
  public:
    ChSuperLUMTEngine();
    ~ChSuperLUMTEngine();

    /// Set problem dimension.
    void SetProblemSize(int pb_size) { m_n = pb_size; }

    /// Set the problem matrix.
    /// This will also update the problem dimension as well as the matrix symmetry type.
    void SetMatrix(ChSparseMatrix& Z);

    /// Set directly the CSR matrix arrays.
    /// Note that it is implied that the matrix symmetry type is GENERAL.
    void SetMatrix(int pb_size, double* values, int* rowIndex, int* colIndex);

    /// Set the solution vector.
    /// Note that it is the caller's responsibility to provide an array of appropriate size.
    void SetSolutionVector(ChMatrix<>& x);
    void SetSolutionVector(double* x);

    /// Set the right-hand side vector.
    /// Note that it is the caller's responsibility to ensure that the size is appropriate.
    void SetRhsVector(ChMatrix<>& b);
    void SetRhsVector(double* b, int nrhs = 1);

    /// Set the matrix, as well as the right-hand side and solution arrays.
    void SetProblem(ChSparseMatrix& Z, ChMatrix<>& b, ChMatrix<>& x);

    /// Solver routine.
	int SuperLUMTCall(phase_t phase, int verbose = 0);

    /// Reinitializes the solver to default values.
    void ResetSolver();

	// Auxiliary functions
	/// Returns the Options vector
	superlumt_options_t& GetOptions() { return superlumt_options; }

	/// Set the number of cores for SuperLU_MT (only factorization phase is parallelized)
	void SetNumProcs(int nprocs_in);

	/// Get the number of cores for SuperLU_MT
	int GetNumProcs() const { return superlumt_options.nprocs; }

	/// Compute reverse of conditioning number of the input matrix
	void SetRCONDevaluation(bool on_off) { rcond_evaluation = on_off; }

	/// Returns the value of the conditioning number of the last evaluation
	double GetRCOND() const { return m_rcond; }

	/// Force iterative refinement of the computed solution
	void SetIterativeRefinements(bool on_off) { iterative_refinement = on_off; }


    // Output functions
    /// Calculate and return the problem residual res=b-Ax.
    /// Note that it is the caller's responsibility to provide an array of appropriate size.
    void GetResidual(ChMatrix<>& res) const;
    void GetResidual(double* res) const;

    /// Calculate and return the L2-norm of the problem residual, ||b-Ax||.
    double GetResidualNorm() const;

  private:

	// Problem properties
	int m_n = 0;     ///< (square) matrix size
	int m_nrhs = 0;  ///< number of rhs vectors

	/* SuperLU_MT data */
	int& ldx = m_n; ///< leading-dimension size of the arrays

	SuperMatrix    m_mat_Super, m_rhs_Super, m_sol_Super;
	std::vector<int> perm_c; ///< column permutation vector
	std::vector<int> perm_r; ///< row permutations from partial pivoting
	std::vector<int> etree; ///< SuperLU internal used

	int            lwork = 0; ///< allocate space internally by system malloc (don't use 'work' variable)

	std::vector<double> R, C;
	std::vector<double> ferr, berr;
	double         rpg = 0, m_rcond = 0;

	// internally used and never directly modified by user
	SuperMatrix    L, U;
	void*          work = nullptr; // don't used to allocate space
	int            info = 0;

	// SuperLU_MT datas (different from serial SuperLU)
	superlu_memusage_t    superlu_memusage;
	superlumt_options_t superlumt_options;
	equed_t     equed = NOEQUIL;
	int relax, panel_size;
	std::vector<int> colcnt_h; ///< SuperLU_MT internal used
	std::vector<int> part_super_h; ///< SuperLU_MT internal used
	Gstat_t   Gstat;

	bool iterative_refinement = false;
	bool rcond_evaluation = false;

protected:
	void
		pdgssvx_mod(int_t nprocs, superlumt_options_t *superlumt_options, SuperMatrix *A,
			int_t *perm_c, int_t *perm_r, equed_t *equed, double *R, double *C,
			SuperMatrix *L, SuperMatrix *U,
			SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth,
			double *rcond, double *ferr, double *berr,
			superlu_memusage_t *superlu_memusage, int_t *info);

	

};




	/// @} superlumt_module

}  // end of namespace chrono

#endif