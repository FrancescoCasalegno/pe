#include "app.h"
#include "utils.h"
#include "bd_type.h" 
#include <petscsys.h>
#include <apf.h>
#include <apfShape.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <cmath>
#include <functional>

namespace {

//auto bd_condition = [](apf::Vector3 const& p)->BoundaryType{ return (p[0]>1.-1.e-12)?(NEUMANN):(DIRICHLET) ; };
auto bd_condition = [](apf::Vector3 const& p)->BoundaryType{ return DIRICHLET; };

auto g_neu = [](apf::Vector3 const& p)->double{ return 1.; };

auto g_dir = [](apf::Vector3 const& p)->double{ return 0.; };
// auto g_dir = [](apf::Vector3 const& p)->double{ return  (p[0]>1.-1.e-12)?(1):(0) ; };

// auto rhs    = [](apf::Vector3 const& p)->double{ return 2.*(p[0]*(1.-p[0])+p[1]*(1.-p[1])); };
// auto rhs = [](apf::Vector3 const& p)->double{ return 1.; };
// auto rhs = [](apf::Vector3 const& p)->double{ return std::sin(2*M_PI*p[0])*std::sin(2*M_PI*p[1]); };

// This RHS corresponds to the exact solution u(x,y) = x^4*(1-x)*y*(1-y)
auto rhs = [](apf::Vector3 const& p)->double{ return -2*p[0]*p[0]*( p[0]*p[0]*p[0] - p[0]*p[0] + 10*p[0]* (p[1]-1) * p[1] - 6*p[1]*(p[1]-1)); };


void initialize()
{
  MPI_Init(0,0);
  PCU_Comm_Init();
  PetscInitialize(0,0,0,0);
}

void finalize()
{
  PetscFinalize();
  PCU_Comm_Free();
  MPI_Finalize();
}

}



int main(int argc, char** argv)
{
  ASSERT(argc == 4);
  const char* geom = argv[1];
  const char* mesh = argv[2];
  const char* out = argv[3];
  const int fem_ord = 1;
  const int integr_ord = 1;  
  initialize();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(geom, mesh);
  pe::AppInput in = { m, fem_ord, integr_ord, bd_condition, g_dir, rhs, out };
  pe::App app(in);
  app.run();
  m->destroyNative();
  apf::destroyMesh(m);
  finalize();
}
