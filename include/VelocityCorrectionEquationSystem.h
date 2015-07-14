/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VelocityCorrectionEquationSystem_h
#define VelocityCorrectionEquationSystem_h

#include <Enums.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradAlgorithmDriver;
class AlgorithmDriver;
class EquationSystems;

class VelocityCorrectionEquationSystem : public EquationSystem {

public:

  VelocityCorrectionEquationSystem(
    EquationSystems& equationSystems,
    const bool managesSolve = false);
  virtual ~VelocityCorrectionEquationSystem();

  void register_nodal_fields(
    stk::mesh::Part *part);

  void register_interior_algorithm(
    stk::mesh::Part *part);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);

  // internal solve and update from EquationSystems
  void solve_and_update();

  // external intended to be called by another EqSystem (used when someone manages PNGEqs)
  void solve_and_update_external();

  void initialize();
  void reinitialize_linear_system();

  // who manages the solve? Often times, this is created by another EqSys
  const bool managesSolve_;

  // internal fields
  VectorFieldType *velocity_;
  VectorFieldType *qTmp_;  
};

} // namespace nalu
} // namespace Sierra

#endif
