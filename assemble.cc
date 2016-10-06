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
    apf::DynamicArray<apf::Node> nodes_diri;
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
    }
    gmi_end(mdl, bdr_its);
    synchronizeEntitySet(m, ent_set);
    getNodesOnEntitySet(m, ent_set, nodes_diri, f_sh);
    size_t n_nodes = nodes_diri.getSize();
    std::vector<long>   v_rows; v_rows.reserve(n_nodes);
    std::vector<double> v_vals; v_vals.reserve(n_nodes);
    apf::Vector3 p;
    for (auto&& nd : nodes_diri) {
            m->getPoint(nd.entity, nd.node, p);
            v_vals.push_back(g_dir(p)); 
            v_rows.push_back(apf::getNumber(gn, nd));      
    }
    ls->diagMatRow(n_nodes, &v_rows[0]);
    ls->setToVector(n_nodes, &v_rows[0], &v_vals[0]);
    ls->synchronize();
}

static void applyNeuBC(
    apf::Mesh* m,
    apf::MeshElement* mesh_el,
    apf::Field* u,
    std::function<double(apf::Vector3 const&)> g_neu,
    LinSys* ls)
{

}

static void applyBCToSystem(
    apf::Mesh* m,
    apf::GlobalNumbering* gn,
    apf::Field* f,
    std::function<double(apf::Vector3 const&)> g_dir,
    LinSys* ls)
{
    auto bd_condition = [](apf::Vector3 const& p)->BoundaryType{return DIRICHLET;}; 
    applyDirBC(m,f,gn, bd_condition, g_dir, ls);
  //auto bd_condition = [](apf::Vector3 const& p)->BoundaryType{ return (p[0]>1.-1.e-12)?(NEUMANN):(DIRICHLET) ; };
  //gmi_model* model = m->getModel();
  //gmi_ent* boundary;
  //gmi_iter* boundaries = gmi_begin(model, m->getDimension()-1);
  //int iter_n = 0; ///....
  //while ((boundary = gmi_next(model, boundaries)))
  //{
  //  apf::DynamicArray<apf::Node> nodes;
  //  apf::ModelEntity* b = reinterpret_cast<apf::ModelEntity*>(boundary);
  //  /// Looping over ...
  //  int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
  //  apf::MeshIterator* it = m->begin(m->getModelType(b));
  //  apf::MeshEntity* mesh_ent;
  //  while(mesh_ent = m->iterate(it))
  //  {
  //    if (m->toModel(mesh_ent)==b) {
  //      apf::MeshElement* mesh_el = apf::createMeshElement(m, mesh_ent);
  //      apf::Vector3 p = apf::getLinearCentroid(m, mesh_ent);
  //      std::string str_t; 
  //      switch (bd_condition(p)) {
  //        case NEUMANN:
  //            // applyNeuBC(...);
  //            str_t = "NEUMANN";
  //            break;
  //        case DIRICHLET:
  //            //apf::NewArray<long> numbers;
  //            //int sz = apf::getElementNumbers(gn, mesh_ent, numbers);
  //            //printf("Rank %d (iter_n=%d) has %d dofs on it: %d, %d\n", rk, iter_n, sz, numbers[0], numbers[1]);
  //            str_t = "DIRICHLET";
  //            break;
  //      }  
  //      printf("Rank %d (iter_n=%d) has zentrum: (%f, %f, %f)----->%s\n", rk, iter_n, p[0],p[1],p[2], str_t.c_str());
  //    }
  //  }
  //  /// ... end looping
  //  apf::getNodesOnClosure(m, b, nodes, apf::getShape(f));
  //  size_t nnodes = nodes.getSize();
  //  std::vector<long>   v_rows; v_rows.reserve(nnodes);
  //  std::vector<double> v_vals; v_vals.reserve(nnodes);
  //  apf::Vector3 p;
  //  for (auto&& nd : nodes) {
  //          m->getPoint(nd.entity, nd.node, p);
  //          v_vals.push_back(g_diri(p)); 
  //          v_rows.push_back(apf::getNumber(gn, nd));      
  //  }
  //  ls->diagMatRow(nnodes, &v_rows[0]);
  //  ls->setToVector(nnodes, &v_rows[0], &v_vals[0]);
  //  ++iter_n;
  //}
  //gmi_end(model, boundaries);
  //ls->synchronize();
}

void App::assemble()
{
  double t0 = PCU_Time();
  assembleSystem(polynomialOrder, mesh, sol, rhs, shared, linsys);
  applyBCToSystem(mesh, shared, sol, g_diri, linsys);
  double t1 = PCU_Time();
  print("assembled in %f seconds", t1-t0);
}

}
