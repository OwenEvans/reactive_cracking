
#include "SystemSolversWrapper.h"
#include "BoostTypes.h"
#include "Logger.h"
#include <dolfin.h>
#include "ElasticitySNES.h"

namespace buckettools
{
  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).
  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname,                                             Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                            const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr functionspace;
    if (systemname ==  "Elasticity")
    {
      // All solvers within a system should return the same functionspace so just take the first one
      if (periodicmap && facetdomains)
      {
        functionspace.reset( new ElasticitySNES::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
      }
      else
      {
        functionspace.reset( new ElasticitySNES::FunctionSpace(mesh) );
      }
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_functionspace", "System name: %s", systemname.c_str());
    }
    return functionspace;
  }

  // A function to return a functionspace from a system given a mesh and a solvername.
  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname, const std::string &solvername, 
                                            Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                            const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr functionspace;
    if (systemname ==  "Elasticity")
    {
      // All solvers within a system should return the same functionspace
      if (solvername ==  "SNES")
      {
        if (periodicmap && facetdomains)
        {
          functionspace.reset(new ElasticitySNES::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
        }
        else
        {
          functionspace.reset(new ElasticitySNES::FunctionSpace(mesh));
        }
      }
      else
      {
        tf_err("Unknown solvername in ufc_fetch_functionspace", "Solver name: %s", solvername.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_functionspace", "System name: %s", systemname.c_str());
    }
    return functionspace;
  }

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_solver(const std::string &systemname, const std::string &solvername, const std::string &uflsymbol, 
                                                           Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                           const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "Elasticity")
    {
      if (solvername ==  "SNES")
      {
        if (uflsymbol ==  "sigmac")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_sigmac(mesh));
        }
        else if (uflsymbol ==  "G")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_G(mesh));
        }
        else if (uflsymbol ==  "Kon2G")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_Kon2G(mesh));
        }
        else if (uflsymbol ==  "u0")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_u0(mesh));
        }
        else if (uflsymbol ==  "u1")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_u1(mesh));
        }
        else if (uflsymbol ==  "s")
        {
          coefficientspace.reset(new ElasticitySNES::CoefficientSpace_s(mesh));
        }
        else
        {
          tf_err("Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver", "UFL symbol: %s", uflsymbol.c_str());
        }
      }
      else
      {
        tf_err("Unknown solvername in ufc_fetch_coefficientspace_from_solver", "Solver name: %s", solvername.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_coefficientspace_from_solver", "System name: %s", systemname.c_str());
    }
    return coefficientspace;
  }

  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.
  Form_ptr ufc_fetch_form(const std::string &systemname, const std::string &solvername, const std::string &solvertype, const std::string &formname, const FunctionSpace_ptr functionspace)
  {
    Form_ptr form;
    if (systemname ==  "Elasticity")
    {
      if (solvername ==  "SNES")
      {
        if (solvertype == "SNES")
        {
          if (formname == "Residual")
          {
            form.reset(new ElasticitySNES::Form_F(functionspace));
          }
          else if (formname == "Jacobian")
          {
            form.reset(new ElasticitySNES::Form_J(functionspace, functionspace));
          }
          else
          {
            tf_err("Unknown formname in ufc_fetch_form", "Form name: %s", formname.c_str());
          }
        }
        else
        {
          tf_err("Unknown solvertype in ufc_fetch_form", "Form name: %s", formname.c_str());
        }
      }
      else
      {
        tf_err("Unknown systemname in ufc_fetch_form", "System name: %s", systemname.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in ufc_fetch_form", "System name: %s", systemname.c_str());
    }
    return form;
  }

}

