#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <math.h>
#include <petsctime.h>
#include "def.h"


#undef __FUNCT__
#define __FUNCT__ "phi_ex"
PetscErrorCode phi_ex(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	PetscErrorCode 	ierr;
	PetscInt       	i_old;
	PetscInt 		order = min(user->curN+1, user->BDFiA);
	int             rank;
	PetscScalar 	BDF_coef[10] = {0};
    PetscScalar 	ex_coef[10] = {0};
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//PetscPrintf(PETSC_COMM_WORLD, "Rank:%d, phi_ex in\n", rank);
    
	cal_BDF_coef(BDF_coef, order);
	cal_ex_coef(ex_coef, order);
	
    for (i_old = 0; i_old < order; ++i_old){
    	ierr = VecAXPY(user->x_sum, BDF_coef[i_old], user->x_olds[i_old]); CHKERRQ(ierr);
	}
	
	for (i_old = 0; i_old < order; ++i_old){
		ierr = VecAXPY(user->x_ex, ex_coef[i_old], user->x_olds[i_old]); CHKERRQ(ierr);
	}

    //PetscPrintf(PETSC_COMM_WORLD, "Rank:%d, phi_ex out\n", rank);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BDFi"
PetscScalar BDFi(PetscInt order, PetscReal* phi_old){
	PetscReal	phi_sum;
	switch(order){
  		case 1:
  			phi_sum = phi_old[0];
			return phi_sum;
  			break;
  		case 2:
  			phi_sum = 4.0/3.0*phi_old[0]
			-1.0/3.0*phi_old[1];
			return phi_sum;
  			break;
  		case 3:
  			phi_sum = 18.0/11.0*phi_old[0]
			  -9.0/11.0*phi_old[1]
			  +2.0/11.0*phi_old[2];
			return phi_sum;
			break;
		default:
			phi_sum = 48.0/25.0*phi_old[0]
			-36.0/25.0*phi_old[1]
			+16.0/25.0*phi_old[2]
			-3.0/25.0*phi_old[3];
			return phi_sum;
	}
	return 0;	
}

#undef __FUNCT__
#define __FUNCT__ "cal_BDF_coef"
PetscErrorCode cal_BDF_coef(PetscScalar *BDF_coef, PetscInt order){
	switch(order){
  		case 1:
  			BDF_coef[0] = 1;
  			break;
  		case 2:
  			BDF_coef[0] = 4.0/3.0;
  			BDF_coef[1] = -1.0/3.0;
  			break;
  		case 3:
  			BDF_coef[0] = 18.0/11.0;
  			BDF_coef[1] = -9.0/11.0;
  			BDF_coef[2] = 2.0/11.0;
			break;
		default:
			BDF_coef[0] = 48.0/25.0;
  			BDF_coef[1] = -36.0/25.0;
  			BDF_coef[2] = 16.0/25.0;
  			BDF_coef[3] = -3.0/25.0;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "cal_ex_coef"
PetscErrorCode cal_ex_coef(PetscScalar *ex_coef, PetscInt order){
	PetscInt i, j, l;
	for (i = 0; i < order; ++i){
		ex_coef[i] = 0.0;
	}
	ex_coef[0] = 1.0;
	l = 1;
	for (i = 1; i < order; ++i){
		for (j = l; j > 0; --j){
			ex_coef[j] -= ex_coef[j-1];
		}
		ex_coef[0] += 1.0;
		++l;
	}
	/*
	for (i = 0; i < order; ++i){
		PetscPrintf(PETSC_COMM_WORLD, "ex_coef[%d] = %g\n", i, ex_coef[i]);
	}
	*/
	PetscFunctionReturn(0);
}
