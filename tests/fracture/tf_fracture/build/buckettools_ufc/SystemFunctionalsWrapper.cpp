
#include "SystemFunctionalsWrapper.h"
#include "BoostTypes.h"
#include "Logger.h"
#include <dolfin.h>


namespace buckettools
{
  // A function to return a functionspace (for a coefficient) from a system given a mesh, a functionalname and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &functionalname, const std::string &uflsymbol, 
                                                               Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                               const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "Elasticity")
    {
        tf_err("No functionals in ufc_fetch_coefficientspace_from_functional", "Functional name: %s", functionalname.c_str());
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_coefficientspace_from_functional", "System name: %s", systemname.c_str());
    }
    return coefficientspace;
  }

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a coefficientname and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_constant_functional(const std::string &systemname, const std::string &coefficientname, const std::string &uflsymbol, 
                                                               Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                               const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "Elasticity")
    {
      if (coefficientname ==  "CriticalStress")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "AlphaLambda")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "G")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "Kon2G")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "u0")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "u1")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else if (coefficientname ==  "StabilizationCoefficient")
      {
        tf_err("No functionals for coefficient in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
      else
      {
        tf_err("Unknown coefficientname in ufc_fetch_coefficientspace_from_constant_functional", "Coefficient name: %s", coefficientname.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_coefficientspace_from_constant_functional", "System name: %s", systemname.c_str());
    }
    return coefficientspace;
  }

  // A function to return a functional from a system-function set given a mesh and a functionalname.
  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &functionalname, Mesh_ptr mesh)
  {
    Form_ptr functional;
    if (systemname ==  "Elasticity")
    {
        tf_err("No functionals in ufc_fetch_functional", "Functional name: %s", functionalname.c_str());
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_functional", "System name: %s", systemname.c_str());
    }
    return functional;
  }

  // A function to return a functional for a constant from a system-function set given a mesh.
  Form_ptr ufc_fetch_constant_functional(const std::string &systemname, const std::string &coefficientname, Mesh_ptr mesh)
  {
    Form_ptr functional;
    if (systemname ==  "Elasticity")
    {
      tf_err("Unknown coefficientname in ufc_fetch_constant_functional", "Coefficient name: %s", coefficientname.c_str());
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_constant_functional", "System name: %s", systemname.c_str());
    }
    return functional;
  }

}

