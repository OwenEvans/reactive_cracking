
#include "VisualizationWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include "_VisualizationOnMeshMesh.h"

namespace buckettools
{
  // A function to return a functionspace for visualization given a mesh and a mesh name.
  FunctionSpace_ptr ufc_fetch_visualization_functionspace(const std::string &meshname, Mesh_ptr mesh)
  {
    FunctionSpace_ptr functionspace;
    if (meshname ==  "Mesh")
    {
      functionspace.reset( new _VisualizationOnMeshMesh::FunctionSpace(mesh) );
    }
    else
    {
      dolfin::error("Unknown meshname in ufc_fetch_visualization_functionspace");
    }
    return functionspace;
  }

}

