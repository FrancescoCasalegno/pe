#ifndef PE_APP_H
#define PE_APP_H

#include "bd_cond.h"
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

class App
{
  public:

    App(apf::Mesh* m, 
        int pol_o, 
        int integr_o, 
        std::function<BoundaryType(apf::Vector3 const&)> bd_cond,  
        std::function<double(apf::Vector3 const&)> neu_fun,  
        std::function<double(apf::Vector3 const&)> dir_fun, 
        std::function<double(apf::Vector3 const&)> rhs_fun, 
        const char* out_name);
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

    std::function<BoundaryType(apf::Vector3 const&)> bd_condition;
    std::function<double(apf::Vector3 const&)> g_neu;
    std::function<double(apf::Vector3 const&)> g_dir;
    std::function<double(apf::Vector3 const&)> rhs;

    const char* out;
};

}

#endif
