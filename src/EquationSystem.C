/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <EquationSystem.h>
#include <AuxFunctionAlgorithm.h>
#include <SolverAlgorithmDriver.h>
#include <InitialConditions.h>
#include <Realm.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <LinearSystem.h>
#include <ConstantAuxFunction.h>
#include <Enums.h>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

EquationSystem::EquationSystem(
  EquationSystems& eqSystems,
  const std::string name) 
  : equationSystems_(eqSystems),
    realm_(eqSystems.realm_),
    name_(name),
    solverAlgDriver_(new SolverAlgorithmDriver(realm_)),
    timerAssemble_(0.0),
    timerLoadComplete_(0.0),
    timerSolve_(0.0),
    timerMisc_(0.0),
    avgLinearIterations_(0.0),
    maxLinearIterations_(0.0),
    minLinearIterations_(1.0e10),
    nonLinearIterationCount_(0),
    reportLinearIterations_(false),
    edgeNodalGradient_(realm_.realmUsesEdges_),
    linsys_(NULL)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- set_nodal_gradient ----------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::set_nodal_gradient(
 const std::string &dofName)
{
  // determine nodal gradient
  std::map<std::string, std::string >::iterator ii =
    realm_.solutionOptions_->nodalGradMap_.find(dofName);
  if ( ii != realm_.solutionOptions_->nodalGradMap_.end() ) {
    if ( ii->second == "edge" ) 
      edgeNodalGradient_ = true;
    else if ( ii->second == "element" )
      edgeNodalGradient_ = false;
    else 
      throw std::runtime_error("only edge or element nodal gradient supported");
  }
}

//--------------------------------------------------------------------------
//-------- system_is_converged ---------------------------------------------
//--------------------------------------------------------------------------
bool
EquationSystem::system_is_converged()
{
  bool isConverged = true;
  if ( NULL != linsys_ ) {
    isConverged = (linsys_->scaledNonLinearResidual() <  convergenceTolerance_ );
  }
  return isConverged;
}

//--------------------------------------------------------------------------
//-------- provide_scaled_norm ---------------------------------------------
//--------------------------------------------------------------------------
double
EquationSystem::provide_scaled_norm()
{
  return ( (NULL != linsys_) ? linsys_->scaledNonLinearResidual() : 0.0 );
}

//--------------------------------------------------------------------------
//-------- dump_eq_time ----------------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::dump_eq_time()
{

  double l_timer[4] = {timerAssemble_, timerLoadComplete_, timerSolve_, timerMisc_};
  double g_min[4] = {};
  double g_max[4] = {};
  double g_sum[4] = {};

  int nprocs = sierra::Env::parallel_size();

  Env::outputP0() << "Timing for Eq: " << name_ << std::endl;

  // get max, min, and sum over processes
  stk::all_reduce_sum(sierra::Env::parallel_comm(), &l_timer[0], &g_sum[0], 4);
  stk::all_reduce_min(sierra::Env::parallel_comm(), &l_timer[0], &g_min[0], 4);
  stk::all_reduce_max(sierra::Env::parallel_comm(), &l_timer[0], &g_max[0], 4);

  // output
  Env::outputP0() << "         assemble --  " << " \tavg: " << g_sum[0]/double(nprocs)
                  << " \tmin: " << g_min[0] << " \tmax: " << g_max[0] << std::endl;
  Env::outputP0() << "    load_complete --  " << " \tavg: " << g_sum[1]/double(nprocs)
                  << " \tmin: " << g_min[1] << " \tmax: " << g_max[1] << std::endl;
  Env::outputP0() << "            solve --  " << " \tavg: " << g_sum[2]/double(nprocs)
                  << " \tmin: " << g_min[2] << " \tmax: " << g_max[2] << std::endl;
  Env::outputP0() << "             misc --  " << " \tavg: " << g_sum[3]/double(nprocs)
                  << " \tmin: " << g_min[3] << " \tmax: " << g_max[3] << std::endl;

  if (reportLinearIterations_)
    Env::outputP0() << "linear iterations -- " << " \tavg: " << avgLinearIterations_
                    << " \tmin: " << minLinearIterations_ << " \tmax: "
                    << maxLinearIterations_ << std::endl;

  // reset anytime these are called
  timerAssemble_ = 0.0;
  timerLoadComplete_ = 0.0;
  timerMisc_ = 0.0;
  timerSolve_ = 0.0;
  avgLinearIterations_ = 0.0;
  minLinearIterations_ = 1.0e10;
  maxLinearIterations_ = 0.0;
  nonLinearIterationCount_ = 0;

}

void
EquationSystem::update_iteration_statistics(
  const int & iters)
{
  const double iterations = (double)iters;
  avgLinearIterations_ = (nonLinearIterationCount_*avgLinearIterations_
                          + iterations)/(nonLinearIterationCount_+1);
  maxLinearIterations_ = std::max(maxLinearIterations_,iterations);
  minLinearIterations_ = std::min(minLinearIterations_,iterations);
  nonLinearIterationCount_ += 1;
  reportLinearIterations_ = true;
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::initial_work()
{
  std::vector<Algorithm *>::iterator ii;
  for( ii=copyStateAlg_.begin(); ii!=copyStateAlg_.end(); ++ii )
    (*ii)->execute();
}

//--------------------------------------------------------------------------
//-------- assemble_and_solve ----------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::assemble_and_solve(
  stk::mesh::FieldBase *deltaSolution)
{

  int error = 0;
  
  // zero the system
  double timeA = stk::cpu_time();
  linsys_->zeroSystem();
  double timeB = stk::cpu_time();
  timerAssemble_ += (timeB-timeA);

  // apply all flux and dirichlet algs
  timeA = stk::cpu_time();
  solverAlgDriver_->execute();
  timeB = stk::cpu_time();
  timerAssemble_ += (timeB-timeA);

  // load complete
  timeA = stk::cpu_time();
  linsys_->loadComplete();
  timeB = stk::cpu_time();
  timerLoadComplete_ += (timeB-timeA);

  // solve the system; extract delta
  timeA = stk::cpu_time();
  error = linsys_->solve(deltaSolution);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_delta_solution_update(deltaSolution, linsys_->numDof());
  }

  timeB = stk::cpu_time();
  timerSolve_ += (timeB-timeA);

  // handle statistics
  update_iteration_statistics(
    linsys_->linearSolveIterations());
  
  if ( error > 0 )
    Env::outputP0() << "Error in " << name_ << "::solve_and_update()  " << std::endl;
  
}


//--------------------------------------------------------------------------
//-------- bc_data_specified ----------------------------------------------------
//--------------------------------------------------------------------------
bool
EquationSystem::bc_data_specified(
  const UserData &userData, std::string &name)
{
  bool isSpecified = false;
  std::map<std::string, bool>::const_iterator iter
    = userData.bcDataSpecifiedMap_.find(name);
  if ( iter != userData.bcDataSpecifiedMap_.end() ) {
    isSpecified = true;
  }
  return isSpecified;
}

//--------------------------------------------------------------------------
//-------- get_bc_data_type ------------------------------------------------
//--------------------------------------------------------------------------
UserDataType
EquationSystem::get_bc_data_type(
  const UserData &userData, std::string &name)
{
  UserDataType dataType = CONSTANT_UD;
  std::map<std::string, UserDataType>::const_iterator iter
    = userData.bcDataTypeMap_.find(name);
  if ( iter != userData.bcDataTypeMap_.end() ) {
    dataType = (*iter).second;
  }
  return dataType;
}

//--------------------------------------------------------------------------
//-------- get_bc_function_name --------------------------------------------
//--------------------------------------------------------------------------
std::string
EquationSystem::get_bc_function_name(
  const UserData &userData, std::string &name)
{
  std::string fcnName = "no_name";
  std::map<std::string, std::string>::const_iterator iter
    = userData.userFunctionMap_.find(name);
  if ( iter != userData.userFunctionMap_.end() ) {
    fcnName = (*iter).second;
  }
  else {
    throw std::runtime_error("wall function name not found");
  }
  
  return fcnName;
}

//--------------------------------------------------------------------------
//-------- get_bc_function_params ------------------------------------------
//--------------------------------------------------------------------------
std::vector<double>
EquationSystem::get_bc_function_params(
  const UserData &userData, std::string &name)
{
  std::vector<double> theParams;
  std::map<std::string, std::vector<double> >::const_iterator iter
    = userData.functionParams_.find(name);
  if ( iter != userData.functionParams_.end() ) {
    return (*iter).second;
  }
  else {
    return theParams;
  }
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EquationSystem::~EquationSystem()
{
  delete solverAlgDriver_;
  
  if ( NULL != linsys_ )
    delete linsys_;

  // initial conditions and bc prop
  std::vector<AuxFunctionAlgorithm *>::iterator ii;
  for( ii=bcDataAlg_.begin(); ii!=bcDataAlg_.end(); ++ii )
    delete *ii;

  std::vector<Algorithm *>::iterator iim;
  for( iim=bcDataMapAlg_.begin(); iim!=bcDataMapAlg_.end(); ++iim )
    delete *iim;

  for( iim=copyStateAlg_.begin(); iim!=copyStateAlg_.end(); ++iim )
    delete *iim;

}


} // namespace nalu
} // namespace Sierra