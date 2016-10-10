#ifndef PE_INTEGRATE_H
#define PE_INTEGRATE_H

#include <apf.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <functional>

namespace pe {

class Integrate : public apf::Integrator
{
  public:
    Integrate(int integr_ord, apf::Field* f, std::function<double(apf::Vector3 const&)> rhs_fun);
    void inElement(apf::MeshElement*) override;
    void outElement() override;
    void atPoint(apf::Vector3 const& p, double w, double dv) override;
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
  private:
    int ndofs;
    int ndims;
    apf::Field* u;
    apf::Element* e;
    std::function<double(apf::Vector3 const&)> rhs;
};

//----------------------
class IntegrateNeuBC : public apf::Integrator
{
public:
    IntegrateNeuBC(int integr_ord, apf::Field* f, std::function<double(apf::Vector3 const&)> g_neu);
    void inElement(apf::MeshElement*) override;
    void outElement() override;
    void atPoint(apf::Vector3 const& p, double w, double dv) override;
    apf::DynamicVector fe;
private:
    int n_dofs;
    int n_dims;
    apf::Field* f;
    apf::Element* e;
    std::function<double(apf::Vector3 const&)> g_neu;

};


}

#endif
