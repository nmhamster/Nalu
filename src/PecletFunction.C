/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <PecletFunction.h>

// basic c++
#include <algorithm>
#include <cmath>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// PecletFunction - base class
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PecletFunction::PecletFunction()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PecletFunction::~PecletFunction()
{
  // nothing to do
}

//==========================================================================
// Class Definition
//==========================================================================
// ClassicPecletFunction - classic 
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ClassicPecletFunction::ClassicPecletFunction( const double A, const double hf)
  : A_(A),
    hf_(hf)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ClassicPecletFunction::~ClassicPecletFunction()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double 
ClassicPecletFunction::execute(const double pecletNumber)
{
  const double modPeclet = hf_*pecletNumber;
  return modPeclet*modPeclet/(5.0 + modPeclet*modPeclet);
}

//==========================================================================
// Class Definition
//==========================================================================
// TanhPecletFunction - classic 
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TanhPecletFunction::TanhPecletFunction(double c1, double c2)
  : c1_(c1),
    c2_(c2),
    shift_(0.0),
    delta_(1.0)
{
  // make sure c2_ is greater than something small
  c2_ = std::max(c2_, 1.0e-16);
  const double pecMin = execute(0.0);
  const double pecMax = execute(1.0e16);
  shift_ = pecMin;
  delta_ = pecMax - pecMin;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TanhPecletFunction::~TanhPecletFunction()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double 
TanhPecletFunction::execute(const double pecletNumber)
{
  return (0.50*(1.0+std::tanh((pecletNumber-c1_)/c2_))-shift_)/delta_;
}

} // namespace nalu
} // namespace Sierra
