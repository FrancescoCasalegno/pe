#ifndef PE_BD_COND_H
#define PE_BD_COND_H

#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi.h>
#include <functional>

enum BoundaryType{
    DIRICHLET,
    NEUMANN
};


std::vector<apf::MeshEntity*> getNeuMeshEntities(
        apf::Mesh* m,
        std::function<BoundaryType(apf::Vector3 const&)> bd_condition
        );


std::vector<apf::Node> getDirNodes( 
        apf::Mesh* m,
        apf::FieldShape* f_sh,
        std::function<BoundaryType(apf::Vector3 const&)> bd_condition
        );


#endif // PE_BD_COND_H
