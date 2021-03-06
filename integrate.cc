#include "integrate.h"
#include <apfMesh.h>
#include <cmath>

namespace pe {

Integrate::Integrate(int integr_ord, apf::Field* f, std::function<double(apf::Vector3 const&)> rhs_fun) :
    apf::Integrator(integr_ord),
    u(f),
    rhs(rhs_fun),
    ndims(apf::getMesh(f)->getDimension())
{
}

void Integrate::inElement(apf::MeshElement* me)
{
  e = apf::createElement(u,me);
  ndofs = apf::countNodes(e);
  fe.setSize(ndofs);
  ke.setSize(ndofs,ndofs);
  for (int a=0; a < ndofs; ++a)
  {
    fe(a) = 0.0;
    for (int b=0; b < ndofs; ++b)
      ke(a,b) = 0.0;
  }
}

void Integrate::outElement()
{
  apf::destroyElement(e);
}

void Integrate::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<double> BF;
  apf::getShapeValues(e,p,BF);

  apf::NewArray<apf::Vector3> gradBF;
  apf::getShapeGrads(e,p,gradBF);

  apf::Vector3 x;
  apf::MeshElement* me = apf::getMeshElement(e);
  apf::mapLocalToGlobal(me,p,x);

  for (int a=0; a < ndofs; ++a)
  {
    fe(a) += rhs(x) * BF[a] * w * dv;
    for (int b=0; b < ndofs; ++b)
    for (int i=0; i < ndims; ++i)
      ke(a,b) += 0.1 * gradBF[a][i] * gradBF[b][i] * w * dv +
                 1./std::sqrt(2) * gradBF[b][i] * BF[a] * w * dv; 
  }
}

//-------------------------
IntegrateNeuBC::IntegrateNeuBC(int integr_ord, apf::Field* f, std::function<double(apf::Vector3 const&)> g_neu) : 
    apf::Integrator(integr_ord),
    f(f),
    g_neu(g_neu),
    n_dims(apf::getMesh(f)->getDimension()-1)
{
}

void IntegrateNeuBC::inElement(apf::MeshElement* me)
{
  e = apf::createElement(f,me);
  n_dofs = apf::countNodes(e);
  fe.setSize(n_dofs);
  for (auto&& fe_i : fe)
    fe_i = 0.0;
}

void IntegrateNeuBC::outElement()
{
  apf::destroyElement(e);
}

void IntegrateNeuBC::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<double> BF;
  apf::getShapeValues(e,p,BF);
  
  apf::Vector3 x;
  apf::MeshElement* me = apf::getMeshElement(e);
  apf::mapLocalToGlobal(me,p,x);
  
  for (int a=0; a<n_dofs; ++a)
    fe(a) += g_neu(x) * BF[a] * w * dv;
}
}
