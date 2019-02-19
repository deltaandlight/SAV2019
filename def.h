#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <math.h>
#include <petsctime.h>

#undef __STRUCT__
#define __STRUCT__ "tsCtx"
typedef struct {
	PetscReal      	dt, curT, endT;
  	PetscInt 		tsmax, tsstep, InterP;
  	PetscReal 		fnorm0, fnorm_ratio, norm_tol;
  	PetscReal		tmin, tmax;
  	PetscReal		alpha, Zalpha;
} tsCtx;

#undef __STRUCT__
#define __STRUCT__ "AppCtx"
typedef struct {
	tsCtx		*ts;
	KSP			ksp;
	DM			da;
	Vec         x_olds[10];
	Vec			x, B, C, X, Y, b;
	Vec			x_ex, x_sum;
	Mat			A;

    PassiveReal Length;          /* test problem parameter */
    PetscInt   	curN, order;
    PetscInt    		BDFi, BDFiA, BDFAorB;
    PetscScalar beta, eps, gamma, c0;
	PetscScalar stept[20000], Tenergy[20000];
	PetscScalar R, Rold, Rold2, Rold3;
    PetscScalar	SE1, SE, SE1_nplus1_2;/*about energy*/
    PetscScalar	SBxphi, SBxphiP;/*temp variable*/
} AppCtx;

extern PetscErrorCode FormInitial(DM,AppCtx*);

extern PetscErrorCode Update(void *ptr);
extern PetscErrorCode SetKSP(void *ptr);
extern PetscErrorCode CalB(void *ptr);
extern PetscErrorCode Getdt(void *ptr);
extern PetscErrorCode CalC(void *ptr);
extern PetscErrorCode Calnewx(void *ptr);

extern PetscErrorCode phi_ex(void *ptr);
extern PetscScalar 	  BDFi(PetscInt, PetscReal*);
extern PetscErrorCode cal_ex_coef(PetscScalar *ex_coef, PetscInt order);
extern PetscErrorCode cal_BDF_coef(PetscScalar *BDF_coef, PetscInt order);

extern PetscScalar 	  max(PetscScalar, PetscScalar);
extern PetscScalar    min(PetscScalar, PetscScalar);
extern PetscErrorCode DataSaveASCII(Vec, char*);
extern PetscErrorCode GetFilename(PetscInt, char*);
extern PetscErrorCode testout(void *ptr);
