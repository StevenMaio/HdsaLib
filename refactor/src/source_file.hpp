#include "trilinos_adapter/core/HDSA_Ptr.hpp"
#include "trilinos_adapter/core/HDSA_Comm.hpp"
#include "trilinos_adapter/core/HDSA_ParameterList.hpp"
#include "trilinos_adapter/core/HDSA_Dense_Matrix.hpp"

#include "core/base/HDSA_Vector.hpp"
#include "core/base/HDSA_Linear_Operator.hpp"
#include "trilinos_adapter/adapters/HDSA_Belos_Adapter.hpp"
#include "trilinos_adapter/core/HDSA_Linear_Algebra.hpp"

#include "core/base/HDSA_Std_Vector.hpp"
#include "core/base/HDSA_MultiVector.hpp"
#include "core/base/HDSA_Stream.hpp"
#include "core/base/HDSA_Randomized_GSVD.hpp"
#include "core/base/HDSA_Hessian_Inversion.hpp"

#include "core/model_discrepancy/HDSA_MD_Interface.hpp"
#include "core/model_discrepancy/HDSA_MD_Interface_Elliptic_Prior.hpp"
#include "core/model_discrepancy/HDSA_Bayes_Posterior_Data.hpp"
#include "core/model_discrepancy/HDSA_MD_Update.hpp"

#include "../interfaces/rol/HDSA_ROL_Vector.hpp"
#include "../interfaces/rol/HDSA_ROL_MD_Interface.hpp"

#include "../interfaces/HDSA_Vector_Mrhyde.hpp"
#include "../interfaces/HDSA_Mrhyde_State_Vector.hpp"
#include "../interfaces/obj_mrhyde.hpp"
#include "../interfaces/MD_Interface_Mrhyde.hpp"



