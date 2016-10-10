#include "app.h"
#include "linsys.h"
#include "utils.h"
#include <PCU.h>

namespace pe {

App::App(apf::Mesh* m, 
        int pol_o, 
        int integr_o, 
        std::function<BoundaryType(apf::Vector3 const&)> bd_cond,  
        std::function<double(apf::Vector3 const&)> neu_fun,  
        std::function<double(apf::Vector3 const&)> dir_fun, 
        std::function<double(apf::Vector3 const&)> rhs_fun, 
        const char* out_name) :
  mesh(m),
  polynomialOrder(pol_o),
  integrationOrder(integr_o),
  bd_condition(bd_cond),
  g_neu(neu_fun),
  g_dir(dir_fun),
  rhs(rhs_fun),
  out(out_name)
{
  print("solvifying poisson's equation!");
}

void App::run()
{
  pre();
  assemble();
  linsys->solve();
  post();
}

}
