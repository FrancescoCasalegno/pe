#include "app.h"
#include "linsys.h"
#include "utils.h"
#include <PCU.h>

namespace pe {

App::App(AppInput& in) :
  mesh(in.mesh),
  polynomialOrder(in.polynomialOrder),
  integrationOrder(in.integrationOrder),
  bd_condition(in.bd_condition),
  g_neu(in.g_neu),
  g_dir(in.g_dir),
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
