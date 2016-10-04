#ifndef PE_APP_H
#define PE_APP_H

#include <functional>

namespace apf {
class Mesh;
class Field;
class Vector3;
template <class T> class NumberingOf;
typedef NumberingOf<long> GlobalNumbering;
}

namespace pe {

class LinSys;

struct AppInput
{
  apf::Mesh* mesh;
  int polynomialOrder;
  int integrationOrder;
  std::function<double(apf::Vector3 const&)> g_diri;
  std::function<double(apf::Vector3 const&)> rhs;
  const char* out;
};

class App
{
  public:

    App(AppInput& in);
    void run();

  private:

    void pre();
    void assemble();
    void post();

    apf::Mesh* mesh;
    apf::Field* sol;
    apf::GlobalNumbering* owned;
    apf::GlobalNumbering* shared;

    int polynomialOrder;
    int integrationOrder;

    LinSys* linsys;

    std::function<double(apf::Vector3 const&)> g_diri;
    std::function<double(apf::Vector3 const&)> rhs;

    const char* out;
};

}

#endif
