/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MasterElement_h
#define MasterElement_h

#include <vector>
#include <cstdlib>
#include <stdexcept>

namespace stk {
struct topology;
}

namespace sierra{
namespace nalu{

namespace Jacobian{

enum Index
{
  DXDS = 0,
  DXDT = 1,
  DYDS = 2,
  DYDT = 3
};

enum Direction
{
  DIR_S = 0,
  DIR_T = 1
};
}

class MasterElement
{
public:

  MasterElement();
  virtual ~MasterElement();

  virtual void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) {
    throw std::runtime_error("determinant not implemented");}

  virtual void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void nodal_grad_op(
    const int nelem,
    double *deriv,
    double * error ) {
    throw std::runtime_error("nodal_grad_op not implemented");}

  virtual void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("face_grad_op not implemented; avoid this element type at open bcs, walls and symms");}
  
  virtual const int * adjacentNodes() {
    throw std::runtime_error("adjacentNodes not implementedunknown bc");
    return NULL;}

  virtual const int * ipNodeMap() {
      throw std::runtime_error("ipNodeMap not implemented");
      return NULL;}

  virtual void shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shape_fcn not implemented"); }

  virtual void shifted_shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shifted_shape_fcn not implemented"); }

  virtual int opposingNodes(
    const int ordinal, const int node) {
    throw std::runtime_error("adjacentNodes not implemented"); }

  virtual int opposingFace(
    const int ordinal, const int node) {
    throw std::runtime_error("opposingFace not implemented"); 
    return 0; }

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord) {
    throw std::runtime_error("isInElement not implemented"); 
    return 1.0e6; }

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result) {
    throw std::runtime_error("interpolatePoint not implemented"); }
  
  virtual void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc) {
    throw std::runtime_error("general_shape_fcn not implement"); }

  virtual void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("general_face_grad_op not implemented");}

  virtual void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords) {
    throw std::runtime_error("sidePcoords_to_elemPcoords");}

  virtual const int * faceNodeOnExtrudedElem() {
    throw std::runtime_error("faceNodeOnExtrudedElem not implement"); }

  virtual const int * opposingNodeOnExtrudedElem() {
    throw std::runtime_error("opposingNodeOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnExtrudedElem() {
    throw std::runtime_error("faceScsIpOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnFaceEdges() {
    throw std::runtime_error("faceScsIpOnFaceEdges not implement"); }

  virtual const double * edgeAlignedArea() {
    throw std::runtime_error("edgeAlignedArea not implement"); }

  double isoparametric_mapping(const double b, const double a, const double xi);

  int nDim_;
  int nodesPerElement_;
  int numIntPoints_;
  double scaleToStandardIsoFac_;

  std::vector<int> lrscv_;
  std::vector<int> ipNodeMap_;
  std::vector<int> oppNode_;
  std::vector<int> oppFace_;
  std::vector<double> intgLoc_;
  std::vector<double> intgLocShift_;
  std::vector<double> intgExpFace_;
  std::vector<double> nodeLoc_;
  // extrusion-based scheme
  std::vector<int> faceNodeOnExtrudedElem_;
  std::vector<int> opposingNodeOnExtrudedElem_;
  std::vector<int> faceScsIpOnExtrudedElem_;
  std::vector<int> faceScsIpOnFaceEdges_;
  std::vector<double> edgeAlignedArea_;
  
};

// Hex 8 subcontrol volume
class HexSCV : public MasterElement
{
public:

  HexSCV();
  virtual ~HexSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error );

};


// Hex 8 subcontrol surface
class HexSCS : public MasterElement
{
public:

  HexSCS();
  virtual ~HexSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // extrusion data structure access
  const int * faceNodeOnExtrudedElem();
  const int * opposingNodeOnExtrudedElem();
  const int * faceScsIpOnExtrudedElem();
  const int * faceScsIpOnFaceEdges();
  const double * edgeAlignedArea();
  
  // helper
  double vector_norm( const double * vect, int len );
  double parametric_distance(const std::vector<double> &x);
  bool within_tol( const double & val, const double & tol );
};


// 3D Quad 27 subcontrol volume
class Hex27SCV : public MasterElement
{
public:
  Hex27SCV();
  ~Hex27SCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  void hex27_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// 3D Hex 27 subcontrol surface
class Hex27SCS : public MasterElement
{
public:
  Hex27SCS();
  ~Hex27SCS();

  const int * adjacentNodes();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shape_fcn(
    double *shpfc);
};

// Tet 4 subcontrol volume
class TetSCV : public MasterElement
{
public:

  TetSCV();
  virtual ~TetSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Tet 4 subcontrol surface
class TetSCS : public MasterElement
{
public:

  TetSCS();
  virtual ~TetSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tet_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );
  
  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);
};

// Pyramid 5 subcontrol volume
class PyrSCV : public MasterElement
{
public:

  PyrSCV();
  virtual ~PyrSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Pyramid 5 subcontrol surface
class PyrSCS : public MasterElement
{
public:

  PyrSCS();
  virtual ~PyrSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void pyr_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);
};

// Wedge 6 subcontrol volume
class WedSCV : public MasterElement
{
public:
  WedSCV();
  virtual ~WedSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Wedge 6 subcontrol surface
class WedSCS : public MasterElement
{
public:
  WedSCS();
  virtual ~WedSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void wedge_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  const int * adjacentNodes();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);
  
  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void wed_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

};

// 2D Quad 4 subcontrol volume
class Quad2DSCV : public MasterElement
{
public:
  Quad2DSCV();
  virtual ~Quad2DSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// 2D Quad 4 subcontrol surface
class Quad2DSCS : public MasterElement
{
public:
  Quad2DSCS();
  virtual ~Quad2DSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  const int * adjacentNodes();

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // extrusion data structure access
  const int * faceNodeOnExtrudedElem();
  const int * opposingNodeOnExtrudedElem();
  const int * faceScsIpOnExtrudedElem();
  const int * faceScsIpOnFaceEdges();
  const double * edgeAlignedArea();
};

// 2D Quad 9 subcontrol volume
class Quad92DSCV : public MasterElement
{
public:
  Quad92DSCV();
  ~Quad92DSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  void quad9_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  double jacobian_determinant(
    const double *elemNodalCoords,
    const double *shapeDerivs);

  void point_shape_deriv(
     const double *isoParCoords,
     double *deriv);

private:
  std::vector<double> ipWeight_;

};

// 2D Quad 9 subcontrol surface
class Quad92DSCS : public MasterElement
{
public:
  Quad92DSCS();
  ~Quad92DSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  const int * adjacentNodes();

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  void quad9_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  void area_vector(
    const double *elemNodalCoords,
    const double *shapeDeriv,
    double *normalVec);

  void point_shape_deriv(
    const double *isoParCoords,
    const Jacobian::Direction direction,
    double *deriv);

private:

  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

  std::vector<ContourData>ipInfo_;

};

// 2D Tri 3 subcontrol volume
class Tri2DSCV : public MasterElement
{
public:
  Tri2DSCV();
  virtual ~Tri2DSCV();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

};

// 2D Tri 3 subcontrol surface
class Tri2DSCS : public MasterElement
{
public:
  Tri2DSCS();
  virtual ~Tri2DSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tri_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);
  
  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  double tri_parametric_distance(
    const std::vector<double> &x);
  
  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

};

// 3D Quad 4
class Quad3DSCS : public MasterElement
{
public:

  Quad3DSCS();
  virtual ~Quad3DSCS();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  bool within_tol( const double & val, const double & tol );
  
  double vector_norm2( const double * vect, int len );

  void non_unit_face_normal(
    const double * par_coord,
    const double * elem_nodal_coor,
    double * normal_vector );
  
  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;
};

// 3D Quad 9
class Quad93DSCS : public MasterElement
{
public:
  Quad93DSCS();
  virtual ~Quad93DSCS();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void area_vector(
    const double *coords,
    const double s,
    double *normalVec);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad9_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// 3D Tri 3
class Tri3DSCS : public MasterElement
{
public:

  Tri3DSCS();
  virtual ~Tri3DSCS();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
     double *shpfc);

   void shifted_shape_fcn(
     double *shpfc);

   void tri_shape_fcn(
     const int &npts,
     const double *par_coord,
     double* shape_fcn);

   double isInElement(
     const double *elemNodalCoord,
     const double *pointCoord,
     double *isoParCoord);

  double parametric_distance(
    const std::vector<double> &x);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

};

// edge 2d
class Edge2DSCS : public MasterElement
{
public:
  Edge2DSCS();
  virtual ~Edge2DSCS();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;  
};

// edge 2d
class Edge32DSCS : public MasterElement
{
public:
  Edge32DSCS();
  virtual ~Edge32DSCS();

  const int * ipNodeMap();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void area_vector(
    const double *coords,
    const double s,
    double *normalVec);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

private:
  std::vector<double> ipWeight_;
};

} // namespace nalu
} // namespace Sierra

#endif
