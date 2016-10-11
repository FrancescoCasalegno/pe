#include "app.h"
#include "utils.h"
#include "linsys.h"
#include "integrate.h"
#include "bd_cond.h"
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

// Assemble Linear System, according to the PDE inside the domain
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


// Modify Linear System, enforcing Dirichlet boundary conditions
static void applyDirBC(
    apf::Mesh* m,
    apf::Field* f,
    apf::GlobalNumbering* gn,
    std::function<BoundaryType(apf::Vector3 const&)> bd_condition,
    std::function<double(apf::Vector3 const&)> g_dir,
    LinSys* ls)
{
    auto vec_dir_nodes = getDirNodes(m, apf::getShape(f), bd_condition);
    size_t n_nodes = vec_dir_nodes.size();
    std::vector<long>   v_rows; v_rows.reserve(n_nodes);
    std::vector<double> v_vals; v_vals.reserve(n_nodes);
    apf::Vector3 p;
    for (auto&& nd : vec_dir_nodes) {
        m->getPoint(nd.entity, nd.node, p);
        v_vals.push_back(g_dir(p)); 
        v_rows.push_back(apf::getNumber(gn, nd));  
    }
    ls->diagMatRow(n_nodes, &v_rows[0]);
    ls->setToVector(n_nodes, &v_rows[0], &v_vals[0]);
    ls->synchronize();
}


// Modify Linear System, enforcing Neumann boundary conditions
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
    auto vec_neu_ents = getNeuMeshEntities(m, bd_condition);
    for (auto&& mesh_ent : vec_neu_ents) {
        apf::MeshElement* mesh_el = apf::createMeshElement(m, mesh_ent);
        integrate_neu_bc.process(mesh_el);
        addToRHS(integrate_neu_bc.fe, mesh_ent, gn, ls);                     
        apf::destroyMeshElement(mesh_el);
    }
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
  assembleSystem(integrationOrder, mesh, sol, rhs, shared, linsys);
  applyNeuBC(integrationOrder, mesh, sol, shared, bd_condition, g_neu, linsys);
  applyDirBC(mesh, sol, shared, bd_condition, g_dir, linsys);
  double t1 = PCU_Time();
  print("assembled in %f seconds", t1-t0);
}

}
