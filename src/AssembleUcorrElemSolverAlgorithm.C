/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleUcorrElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleUcorrElemSolverAlgorithm - add LHS/RHS for uvw momentum pressureP
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleUcorrElemSolverAlgorithm::AssembleUcorrElemSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem)
  : SolverAlgorithm(realm, part, eqSystem),
    velocity_(NULL),
    provisionalVelocity_(NULL),
    density_(NULL),
    GjpK_(NULL),
    GjpNew_(NULL),
    scVolume_(NULL)
{
  // save off data; old pressure gradient stored in uTmp (LowMachEquationSystem::project_nodal_gradient)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  provisionalVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "provisional_velocity");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  GjpK_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "uTmp");
  GjpNew_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  scVolume_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "sc_volume");
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleUcorrElemSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleUcorrElemSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // space for LHS/RHS; nodesPerElem*nDim*nodesPerElem*nDim and nodesPerElem*nDim
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<stk::mesh::Entity> connected_nodes;

  // nodal fields to gather
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_provisionalVelocity;
  std::vector<double> ws_densityNp1;
  std::vector<double> ws_GjpK;
  std::vector<double> ws_GjpNew;

  // geometry related to populate
  std::vector<double> ws_shape_function_scv;

  // fixed size with pointers
  std::vector<double> uNp1Scv(nDim,0.0);
  std::vector<double> uProvScv(nDim,0.0);
  std::vector<double> GjpKScv(nDim,0.0);
  std::vector<double> GjpNewScv(nDim,0.0);
  double *p_uNp1Scv = &uNp1Scv[0];
  double *p_uProvScv = &uProvScv[0];
  double *p_GjpKScv = &GjpKScv[0];
  double *p_GjpNewScv = &GjpNewScv[0];
  
  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = realm_.get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;

    // mappings for this element, SCV
    const int *ipNodeMap = meSCV->ipNodeMap();

    // resize some things; matrix related
    const int lhsSize = nodesPerElement*nDim*nodesPerElement*nDim;
    const int rhsSize = nodesPerElement*nDim;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    connected_nodes.resize(nodesPerElement);

    // algorithm related
    ws_velocityNp1.resize(nodesPerElement*nDim);
    ws_provisionalVelocity.resize(nodesPerElement*nDim);
    ws_densityNp1.resize(nodesPerElement);
    ws_GjpK.resize(nodesPerElement*nDim);
    ws_GjpNew.resize(nodesPerElement*nDim);
    ws_shape_function_scv.resize(numScvIp*nodesPerElement);

    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_provisionalVelocity = &ws_provisionalVelocity[0];
    double *p_densityNp1 = &ws_densityNp1[0];
    double *p_GjpK = &ws_GjpK[0];
    double *p_GjpNew = &ws_GjpNew[0];
    double *p_shape_function_scv = &ws_shape_function_scv[0];

    // extract shape function
    meSCV->shape_fcn(&p_shape_function_scv[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get elem
      stk::mesh::Entity elem = b[k];

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // ip data for this element; scs and scv
      const double *scVolume = stk::mesh::field_data(*scVolume_, b, k );

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // set connected nodes
        connected_nodes[ni] = node;

        // pointers to real data; vector and scalar
        const double * uNp1  =  stk::mesh::field_data(velocityNp1, node);
        const double * uProv  =  stk::mesh::field_data(*provisionalVelocity_, node);
        const double * GjpK =  stk::mesh::field_data(*GjpK_, node);
        const double * GjpNew =  stk::mesh::field_data(*GjpNew_, node);
        const double  rhoNp1  =  *stk::mesh::field_data(densityNp1, node);
      
        // gather scalars
        p_densityNp1[ni] = rhoNp1;
  
        // gather vectors
        const int niNdim = ni*nDim;

        for ( int i=0; i < nDim; ++i ) {
          p_velocityNp1[niNdim+i] = uNp1[i];
          p_provisionalVelocity[niNdim+i] = uProv[i];
          p_GjpK[niNdim+i] = GjpK[i];
          p_GjpNew[niNdim+i] = GjpNew[i];
        }
      }

      // handle scv RHS/LHS
      for ( int ip = 0; ip < numScvIp; ++ip ) {

        // nearest node to ip
        const int nearestNode = ipNodeMap[ip];

        // save off some offsets and sc_volume at this ip
        const int nnNdim = nearestNode*nDim;
        const int offSetSF = ip*nodesPerElement;
        const double scV = scVolume[ip];

        // zero out scv
        for ( int j = 0; j < nDim; ++j ) {
          p_uNp1Scv[j] = 0.0;
          p_uProvScv[j] = 0.0;
          p_GjpKScv[j] = 0.0;
          p_GjpNewScv[j] = 0.0;
        }
        double rhoScv = 0.0;

        // interpolate
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function_scv[offSetSF+ic];
          rhoScv += r*p_densityNp1[ic];
          for ( int j = 0; j < nDim; ++j ) {
            p_uNp1Scv[j] += r*p_velocityNp1[ic*nDim+j];
            p_uProvScv[j] += r*p_provisionalVelocity[ic*nDim+j];
            p_GjpKScv[j] += r*p_GjpK[ic*nDim+j];
            p_GjpNewScv[j] += r*p_GjpNew[ic*nDim+j];
          }
        }
      
        // assemble rhs
        for ( int i = 0; i < nDim; ++i ) {
          const double residual = -((p_uNp1Scv[i] - p_uProvScv[i]) + projTimeScale/rhoScv*(p_GjpNewScv[i] - p_GjpKScv[i])); 
          p_rhs[nnNdim+i] -= residual*scV;
        }

        // manage LHS
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          
          const int icNdim = ic*nDim;
      
          // save off shape function
          const double r = p_shape_function_scv[offSetSF+ic];
      
          const double lhsfac = r*scV;
      
          for ( int i = 0; i < nDim; ++i ) {
            const int indexNN = nnNdim + i;
            const int rowNN = indexNN*nodesPerElement*nDim;
            const int rNNiC_i = rowNN+icNdim+i;
            lhs[rNNiC_i] += lhsfac;
          }
        }
      }

      apply_coeff(connected_nodes, rhs, lhs, __FILE__);

    }
  }
}

} // namespace nalu
} // namespace Sierra
