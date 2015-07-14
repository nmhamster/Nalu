/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleUcorrElemSolverAlgorithm_h
#define AssembleUcorrElemSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleUcorrElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleUcorrElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleUcorrElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  VectorFieldType *velocity_;
  VectorFieldType *provisionalVelocity_;
  ScalarFieldType *density_;
  VectorFieldType *GjpK_;
  VectorFieldType *GjpNew_;
  GenericFieldType *scVolume_;
};

} // namespace nalu
} // namespace Sierra

#endif
