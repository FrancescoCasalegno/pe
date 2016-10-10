#include "app.h"
#include "utils.h"
#include "linsys.h"
#include "integrate.h"
#include "bd_type.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <gmi.h>
#include <apfShape.h>
#include <PCU.h>
#include <functional>
#include <mpi.h>
#include <string>
#include <cassert>
namespace pe {

static void addToRHS(
    apf::DynamicVector& fe,
    apf::MeshEntity* e,
    apf::GlobalNumbering* n,
    LinSys* ls)
{
  apf::NewArray<long> numbers;
  int sz = apf::getElementNumbers(n, e, numbers);
  ls->addToVector(sz, &numbers[0], &fe[0]);
}

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
  Integrate integrate(o, f, rhs);
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


////-----------------------------------------------------------------------------
typedef std::set<apf::MeshEntity*> EntitySet;

static void getClosureEntitiesWithNodes(
    apf::Mesh* m,
    apf::MeshEntity* e,
    EntitySet& out,
    apf::FieldShape* s)
{
  int D = getDimension(m, e);
  for (int d=0; d <= D; ++d)
    if (s->hasNodesIn(d))
    {
      apf::Downward de;
      int nde = m->getDownward(e,d,de);
      for (int i=0; i < nde; ++i)
      out.insert(de[i]);
    }
}

static void synchronizeEntitySet(
    apf::Mesh* m,
    EntitySet& set)
{
  PCU_Comm_Begin();
  APF_ITERATE(EntitySet,set,it)
    if (m->isShared(*it))
    {
      apf::Copies remotes;
      m->getRemotes(*it,remotes);
      APF_ITERATE(apf::Copies,remotes,rit)
        PCU_COMM_PACK(rit->first,rit->second);
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
  {
    apf::MeshEntity* e;
    PCU_COMM_UNPACK(e);
    set.insert(e);
  }
}

static void getNodesOnEntitySet(
    apf::Mesh* m,
    EntitySet& s,
    apf::DynamicArray<apf::Node>& n,
    apf::FieldShape* sh)
{
  size_t size = 0;
  APF_ITERATE(EntitySet,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    size += nen;
  }
  n.setSize(size);
  size_t i=0;
  APF_ITERATE(EntitySet,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    for (int j=0; j < nen; ++j)
      n[i++] = apf::Node(*it,j);
  }
  assert(i==size);
}
////-----------------------------------------------------------------------------

static void applyDirBC(
    apf::Mesh* m,
    apf::Field* f,
    apf::GlobalNumbering* gn,
    std::function<BoundaryType(apf::Vector3 const&)> bd_condition,
    std::function<double(apf::Vector3 const&)> g_dir,
    LinSys* ls)
{
    apf::DynamicArray<apf::Node> nodes_dir;
    apf::FieldShape* f_sh = apf::getShape(f);
    gmi_model* mdl = m->getModel();
    gmi_ent* bdr_it;
    gmi_iter* bdr_its = gmi_begin(mdl, m->getDimension()-1);
    EntitySet ent_set;
    while (bdr_it = gmi_next(mdl, bdr_its)) { // for each model bdr...
        apf::ModelEntity* bdr = reinterpret_cast<apf::ModelEntity*>(bdr_it);
        apf::MeshIterator* it = m->begin(m->getModelType(bdr)); 
        apf::MeshEntity* mesh_ent;
        while (mesh_ent = m->iterate(it)) { // for each entity in the model bdr...
            if (m->toModel(mesh_ent)==bdr) {
                if (bd_condition(apf::getLinearCentroid(m, mesh_ent))==DIRICHLET) { 
                    getClosureEntitiesWithNodes(m, mesh_ent, ent_set, f_sh);
                }
            }
        }
        m->end(it);
    }
    gmi_end(mdl, bdr_its);
    synchronizeEntitySet(m, ent_set);
    getNodesOnEntitySet(m, ent_set, nodes_dir, f_sh);
    size_t n_nodes = nodes_dir.getSize();
    std::vector<long>   v_rows; v_rows.reserve(n_nodes);
    std::vector<double> v_vals; v_vals.reserve(n_nodes);
    apf::Vector3 p;
    for (auto&& nd : nodes_dir) {
        m->getPoint(nd.entity, nd.node, p);
        v_vals.push_back(g_dir(p)); 
        v_rows.push_back(apf::getNumber(gn, nd));  
    }
    ls->diagMatRow(n_nodes, &v_rows[0]);
    ls->setToVector(n_nodes, &v_rows[0], &v_vals[0]);
    ls->synchronize();
}


static void applyNeuBC(
    int integr_ord,
    apf::Mesh* m,
    apf::Field* f,
    apf::GlobalNumbering* gn,
    std::function<BoundaryType(apf::Vector3 const&)> bd_condition,
    std::function<double(apf::Vector3 const&)> g_neu,
    LinSys* ls)
{
    IntegrateNeuBC integrate_neu_bc(integr_ord, f, g_neu);
    apf::FieldShape* f_sh = apf::getShape(f);
    gmi_model* mdl = m->getModel();
    gmi_ent* bdr_it;
    gmi_iter* bdr_its = gmi_begin(mdl, m->getDimension()-1);
    EntitySet ent_set;
    while (bdr_it = gmi_next(mdl, bdr_its)) { // for each model bdr...
        apf::ModelEntity* bdr = reinterpret_cast<apf::ModelEntity*>(bdr_it);
        apf::MeshIterator* it = m->begin(m->getModelType(bdr)); 
        apf::MeshEntity* mesh_ent;
        while (mesh_ent = m->iterate(it)) { // for each entity in the model bdr...
            if (m->toModel(mesh_ent)==bdr) {
                if (bd_condition(apf::getLinearCentroid(m, mesh_ent))==NEUMANN) { 
                    int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
                    apf::MeshElement* mesh_el = apf::createMeshElement(m, mesh_ent);
                    integrate_neu_bc.process(mesh_el);
                    addToRHS(integrate_neu_bc.fe, mesh_ent, gn, ls);                     
                    apf::destroyMeshElement(mesh_el);
                }
            }
        }
        m->end(it);
    }
    gmi_end(mdl, bdr_its);
    ls->synchronize();
}

static void applyBCToSystem(
    apf::Mesh* m,
    apf::GlobalNumbering* gn,
    apf::Field* f,
    std::function<BoundaryType(apf::Vector3 const&)> bd_condition,
    std::function<double(apf::Vector3 const&)> g_neu,
    std::function<double(apf::Vector3 const&)> g_dir,
    LinSys* ls)
{
    int integr_ord = 1;
    applyNeuBC(integr_ord, m, f, gn, bd_condition, g_neu, ls);
    applyDirBC(m,f,gn, bd_condition, g_dir, ls);
}

void App::assemble()
{
  double t0 = PCU_Time();
  assembleSystem(polynomialOrder, mesh, sol, rhs, shared, linsys);
  applyBCToSystem(mesh, shared, sol, bd_condition, g_neu, g_dir, linsys);
  double t1 = PCU_Time();
  print("assembled in %f seconds", t1-t0);
}

}
