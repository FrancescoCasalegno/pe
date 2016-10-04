#include "app.h"
#include "linsys.h"
#include "utils.h"
#include <PCU.h>

namespace pe {

App::App(AppInput& in) :
  mesh(in.mesh),
  polynomialOrder(in.polynomialOrder),
  integrationOrder(in.integrationOrder),
  rhs(in.rhs),
  out(in.out)
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
