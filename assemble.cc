#include "app.h"
#include "utils.h"
#include "linsys.h"
#include "integrate.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <gmi.h>
#include <PCU.h>
#include <functional>

namespace pe {

static void addToSystem(
    apf::DynamicVector& fe,
    apf::DynamicMatrix& ke,
    apf::MeshEntity* e,
    apf::GlobalNumbering* n,
    LinSys* ls)
{
  apf::NewArray<long> numbers;
  int sz = apf::getElementNumbers(n, e, numbers);
  ls->addToVector(sz, &numbers[0], &fe[0]);
  ls->addToMatrix(sz, &numbers[0], &ke(0,0));
}

static void assembleSystem(
    int o,
    apf::Mesh* m,
    apf::Field* f,
    std::function<double(apf::Vector3 const&)> rhs,
    apf::GlobalNumbering* n,
    LinSys* ls)
{
  IntegrateInput in = { o, f, rhs};
  Integrate integrate(in);
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(m->getDimension());
  while ((elem = m->iterate(elems)))
  {
    apf::MeshElement* me = apf::createMeshElement(m, elem);
    integrate.process(me);
    addToSystem(integrate.fe, integrate.ke, elem, n, ls);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  ls->synchronize();
}

static void applyBCToSystem(
    apf::Mesh* m,
    apf::GlobalNumbering* gn,
    apf::Field* f,
    LinSys* ls)
{
  auto g_diri = [](apf::Vector3 const &p)->double { return ((p[0]<1.e-12)||(p[0]>1-1.e-12))?1.:0. ;};
  gmi_model* model = m->getModel();
  gmi_ent* boundary;
  gmi_iter* boundaries = gmi_begin(model, m->getDimension()-1);
  while ((boundary = gmi_next(model, boundaries)))
  {
    apf::DynamicArray<apf::Node> nodes;
    apf::ModelEntity* b = reinterpret_cast<apf::ModelEntity*>(boundary);
    apf::getNodesOnClosure(m, b, nodes, apf::getShape(f));
    size_t nnodes = nodes.getSize();
    std::vector<long>   v_rows; v_rows.reserve(nnodes);
    std::vector<double> v_vals; v_vals.reserve(nnodes);
    apf::Vector3 p;
    for (auto&& nd : nodes) {
            m->getPoint(nd.entity, nd.node, p);
            v_vals.push_back(g_diri(p));  
            v_rows.push_back(apf::getNumber(gn, nd));      
    }
    ls->diagMatRow(nnodes, &v_rows[0]);
    ls->setToVector(nnodes, &v_rows[0], &v_vals[0]);
  }
  gmi_end(model, boundaries);
  ls->synchronize();
}

void App::assemble()
{
  double t0 = PCU_Time();
  assembleSystem(polynomialOrder, mesh, sol, rhs, shared, linsys);
  applyBCToSystem(mesh, shared, sol, linsys);
  double t1 = PCU_Time();
  print("assembled in %f seconds", t1-t0);
}

}
