
#include "SystemExpressionsWrapper.h"
#include "BoostTypes.h"
#include "Logger.h"
#include <dolfin.h>


namespace buckettools
{
  // A function to return an expression for a coefficient from a system given a systemname and a functionname (and its size, shape and private members bucket, system and time.
  Expression_ptr cpp_fetch_expression(const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname, const std::size_t &size, const std::vector<std::size_t> &shape, const Bucket *bucket, const SystemBucket *system, const double_ptr time)
  {
    Expression_ptr expression;
    if (systemname ==  "Elasticity")
    {
      if (functionname == "u")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "CriticalStress")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "AlphaLambda")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "G")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "Kon2G")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "u0")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "u1")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "StabilizationCoefficient")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_fetch_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else
      {
        tf_err("Unknown functionname in cpp_fetch_expression.", "Function name: %s", functionname.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in cpp_fetch_expression", "System name: %s", systemname.c_str());
    }
    return expression;
  }

  // A function to initialize an expression for a cpp expression given a systemname and a functionname (and a boost shared pointer to the expression to initialize.
  void cpp_init_expression(Expression_ptr expression, const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname)
  {
    if (systemname ==  "Elasticity")
    {
      if (functionname == "u")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "CriticalStress")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "AlphaLambda")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "G")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "Kon2G")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "u0")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "u1")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else if (functionname ==  "StabilizationCoefficient")
      {
        if (expressiontype == "initial_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "boundary_condition")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else if (expressiontype == "value")
        {
          tf_err("Unknown expressionname in cpp_fetch_expression.", "Expression name: %s", functionname.c_str());
        }
        else
        {
          tf_err("Unknown expressiontype in cpp_init_expression.", "Expression type: %s", expressiontype.c_str());
        }
      }
      else
      {
        tf_err("Unknown functionname in cpp_init_expression.", "Function name: %s", functionname.c_str());
      }
    }
    else
    {
      tf_err("Unknown systemname in cpp_init_expression", "System name: %s", systemname.c_str());
    }
  }

}

