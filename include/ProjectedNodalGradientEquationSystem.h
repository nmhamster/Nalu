/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ProjectedNodalGradientEquationSystem_h
#define ProjectedNodalGradientEquationSystem_h

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

class ProjectedNodalGradientEquationSystem : public EquationSystem {

public:

  ProjectedNodalGradientEquationSystem(
    EquationSystems& equationSystems);
  virtual ~ProjectedNodalGradientEquationSystem();

  void register_nodal_fields(
    stk::mesh::Part *part);

  void register_interior_algorithm(
    stk::mesh::Part *part);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  void solve_and_update();

  void initialize();
  void reinitialize_linear_system();

  VectorFieldType *dqdx_;
  VectorFieldType *qTmp_;  
};

} // namespace nalu
} // namespace Sierra

#endif
