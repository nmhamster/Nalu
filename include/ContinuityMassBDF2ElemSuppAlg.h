/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassBDF2ElemSuppAlg_h
#define ContinuityMassBDF2ElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class ContinuityMassBDF2ElemSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassBDF2ElemSuppAlg(
    Realm &realm);

  virtual ~ContinuityMassBDF2ElemSuppAlg() {}

  virtual void setup();

  virtual void elem_resize(
    MasterElement *meSCS,
    MasterElement *meSCV);

  virtual void elem_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    MasterElement *meSCS,
    MasterElement *meSCV);
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  GenericFieldType *scVolume_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const bool useShifted_;

  // scratch space
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_qm1_;
  std::vector<double> ws_qN_;
  std::vector<double> ws_qNp1_;
  std::vector<double> ws_rhoNm1_;
  std::vector<double> ws_rhoN_;
  std::vector<double> ws_rhoNp1_;
};

} // namespace nalu
} // namespace Sierra

#endif
