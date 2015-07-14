/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <VelocityCorrectionEquationSystem.h>

#include <AssembleUcorrElemSolverAlgorithm.h>
#include <DirichletBC.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <NaluEnv.h>
#include <Realm.h>
#include <Realms.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <SolverAlgorithmDriver.h>

// user functions
#include <user_functions/SteadyThermalContactAuxFunction.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/CPUTime.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// VelocityCorrectionEquationSystem - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
VelocityCorrectionEquationSystem::VelocityCorrectionEquationSystem(
 EquationSystems& eqSystems,
 const bool managesSolve)
  : EquationSystem(eqSystems, "VelCorrEQS"),
    managesSolve_(managesSolve),
    velocity_(NULL),
    qTmp_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity_correction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_UCORR);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, "VelCorrEQS", solver);

  // push back EQ to manager
  realm_.equationSystems_.push_back(this);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
VelocityCorrectionEquationSystem::~VelocityCorrectionEquationSystem()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();
 
  // register dof (same as Momentum)
  velocity_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity", numStates));
  stk::mesh::put_field(*velocity_, *part, nDim);
 
  // delta solution for linear solver
  qTmp_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "ucTmp"));
  stk::mesh::put_field(*qTmp_, *part, nDim);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // algorithm type
  const AlgorithmType algType = INTERIOR;

  // solver
  std::map<AlgorithmType, SolverAlgorithm *>::iterator its
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleUcorrElemSolverAlgorithm *theAlg
      = new AssembleUcorrElemSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    its->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{
  // algorithm type
  const AlgorithmType algType = WALL;

  // find out if this is a wall function approach
  WallUserData userData = wallBCData.userData_;
  const bool wallFunctionApproach = userData.wallFunctionApproach_;
  
  if ( !wallFunctionApproach  ) {
    stk::mesh::MetaData &meta_data = realm_.meta_data();
    const unsigned nDim = meta_data.spatial_dimension();
    
    // extract fields
    VectorFieldType *velocityBC = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_bc");
    
    // deal with state
    VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
      
    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &velocityNp1, velocityBC, 0, nDim);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }  
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &/*inflowBCData*/)
{
  // algorithm type
  const AlgorithmType algType = INFLOW;
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *velocityBC = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_bc");

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
    
  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &velocityNp1, velocityBC, 0, nDim);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::reinitialize_linear_system()
{
  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_UCORR;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity_correction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_UCORR);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, "VelCorrEQS", solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::solve_and_update()
{
  if ( managesSolve_ )
    solve_and_update_external();
}

//--------------------------------------------------------------------------
//-------- solve_and_update_external ------------------------------------------------
//--------------------------------------------------------------------------
void
VelocityCorrectionEquationSystem::solve_and_update_external()
{
  maxIterations_ = 1;
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;
    
    // projected nodal gradient, load_complete and solve
    assemble_and_solve(qTmp_);
    
    // update
    double timeA = stk::cpu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *qTmp_,
      1.0, velocity_->field_of_state(stk::mesh::StateNP1), 
      realm_.get_activate_aura());
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);   
  }
}

} // namespace nalu
} // namespace Sierra
