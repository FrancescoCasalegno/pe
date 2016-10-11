#include <PCU.h>
#include "bd_cond.h"

//-----------------------------------------------------------------------------
std::vector<apf::MeshEntity*> getNeuMeshEntities(
        apf::Mesh* m,
        std::function<BoundaryType(apf::Vector3 const&)> bd_condition
        )
{
    gmi_model* mdl = m->getModel();
    gmi_ent* bdr_it;
    gmi_iter* bdr_its = gmi_begin(mdl, m->getDimension()-1);
    std::vector<apf::MeshEntity*> ent_vec;

    while (bdr_it = gmi_next(mdl, bdr_its)) { 
        apf::ModelEntity* bdr = reinterpret_cast<apf::ModelEntity*>(bdr_it);
        apf::MeshIterator* it = m->begin(m->getModelType(bdr)); 
        apf::MeshEntity* mesh_ent;
        while (mesh_ent = m->iterate(it)) { 
            if (m->toModel(mesh_ent)==bdr) {
                if (bd_condition(apf::getLinearCentroid(m, mesh_ent))==NEUMANN) { 
                    ent_vec.push_back(mesh_ent);
                }
            }
        }
        m->end(it);
    }
    gmi_end(mdl, bdr_its);

    return ent_vec;
}

//-----------------------------------------------------------------------------
static void getClosureEntitiesWithNodes(
    apf::Mesh* m,
    apf::MeshEntity* e,
    std::set<apf::MeshEntity*>& out,
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
void synchronizeEntSet(
    apf::Mesh* m,
    std::set<apf::MeshEntity*>& set)
{
  PCU_Comm_Begin();
  APF_ITERATE(std::set<apf::MeshEntity*>,set,it)
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

std::vector<apf::Node> getNodesOnEntSet(
    apf::Mesh* m,
    std::set<apf::MeshEntity*>& s,
    apf::FieldShape* sh)
{
  std::vector<apf::Node> vec_neu_nodes;
  size_t size = 0;
  APF_ITERATE(std::set<apf::MeshEntity*>,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    size += nen;
  }
  vec_neu_nodes.reserve(size);
  APF_ITERATE(std::set<apf::MeshEntity*>,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    for (int j=0; j < nen; ++j)
        vec_neu_nodes.emplace_back(*it,j);
  }
  return vec_neu_nodes;
}

//-----------------------------------------------------------------------------
std::vector<apf::Node> getDirNodes( 
        apf::Mesh* m,
        apf::FieldShape* f_sh,
        std::function<BoundaryType(apf::Vector3 const&)> bd_condition
        )
{
    gmi_model* mdl = m->getModel();
    gmi_ent* bdr_it;
    gmi_iter* bdr_its = gmi_begin(mdl, m->getDimension()-1);
    std::set<apf::MeshEntity*> ent_set;
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
    synchronizeEntSet(m, ent_set);
    return getNodesOnEntSet(m, ent_set, f_sh);
}
