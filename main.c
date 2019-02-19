#undef __CONST__
#define __CONST__ "help"
static char help[] = "We solve the  equation in 2D rectangular domain with\n\
the SAV method\n\n";

#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <math.h>
#include <petsctime.h>
#include "def.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{                        
    tsCtx          ts;
  	AppCtx         user;                /* user-defined work context */
  	PetscErrorCode ierr;
  	MPI_Comm       comm;
  	int 		   rank;
    
    PetscLogDouble t1,t2,t3,t4;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return(1);
    
    PetscFunctionBeginUser;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    ierr = PetscMalloc(sizeof(AppCtx), &user); CHKERRQ(ierr);
    /*
     Create DM
     */
    ierr = DMDACreate2d(comm, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX,-4,-4, PETSC_DECIDE, PETSC_DECIDE,1,2,0,0, &user.da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(user.da); CHKERRQ(ierr);
    ierr = DMSetUp(user.da); CHKERRQ(ierr);
    ierr = DMCreateMatrix(user.da, &user.A);
    /*ierr = KSPCreate(PETSC_COMM_WORLD, &user.ksp); CHKERRQ(ierr);*//*put it in SetKSP*/
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize problem parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
    user.BDFi = 4;
    user.BDFAorB = 1;
    user.beta = 0.0;
    user.eps = 0.1;
    user.Length = 2*M_PI;
    user.c0 = 0;
    
    ts.endT = 0.032;
    ts.curT = 0.0;
    ts.InterP = 10;
    ts.tmin = 0.0001;
	ts.tmax = 0.0001;
	ts.dt = 0.0001;
	ts.Zalpha = 400000;
	ts.alpha = 4000000;
    
    ierr = PetscOptionsGetInt(NULL, "-BDFi", &user.BDFi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, "-BDFAorB", &user.BDFAorB, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, "-interP", &ts.InterP, NULL); CHKERRQ(ierr);
    
    ierr = PetscOptionsGetReal(NULL, "-beta", &user.beta, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-eps", &user.eps, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-c0", &user.c0, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, "-Length", &user.Length, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-endT", &ts.endT, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-Tmin", &ts.tmin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-Tmax", &ts.tmax, NULL); CHKERRQ(ierr);
    
    ierr = PetscPrintf(comm, "discrete form = BDF%d%c\n", (int)user.BDFi, (char)(user.BDFAorB+'A')); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "beta          = %g\n", (double)user.beta); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "epsilon       = %g\n", (double)user.eps); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "L             = %g\n", (double)user.Length); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "endTime       = %g\n", (double)ts.endT); CHKERRQ(ierr);
    ierr = PetscPrintf(comm, "interP        = %d\n", ts.InterP); CHKERRQ(ierr);
    
    user.ts = &ts;
    user.BDFiA = user.BDFi - user.BDFAorB;
    user.order = 1;
    ierr = DMSetApplicationContext(user.da, &user); CHKERRQ(ierr);

    /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA; then duplicate for remaining
     vectors that are the same types
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMCreateGlobalVector(user.da, &user.x); CHKERRQ(ierr);
	
	PetscInt i_old;
	for (i_old = 0; i_old < user.BDFi; ++i_old){
		ierr = VecDuplicate(user.x, &user.x_olds[i_old]); CHKERRQ(ierr);
	}
    
    ierr = VecDuplicate(user.x, &user.x_ex); CHKERRQ(ierr);
    ierr = VecDuplicate(user.x, &user.x_sum); CHKERRQ(ierr);

    ierr = PetscOptionsGetReal(NULL,"-dt", &ts.dt, NULL); CHKERRQ(ierr);
    
	ierr = FormInitial(user.da, &user); CHKERRQ(ierr);
	for (i_old = 0; i_old < user.BDFi; ++i_old){
		ierr =  VecCopy(user.x, user.x_olds[i_old]); CHKERRQ(ierr);
	}

	ierr = VecCopy(user.x, user.x_ex); CHKERRQ(ierr);
	ierr = VecCopy(user.x, user.x_sum); CHKERRQ(ierr);

    char  filenameascii[PETSC_MAX_PATH_LEN] = "ASCIsol_0000";
    ierr = DataSaveASCII(user.x, filenameascii);
    
	user.curN = 0;
    PetscTime(&t1);
    do{
        PetscTime(&t3);
        Update(&user);
        PetscTime(&t4);
		PetscPrintf(PETSC_COMM_WORLD,"INFO:\n time of this step equals %f seconds.\n",t4-t3);
        user.stept[user.curN] = ts.dt;
        user.curN = user.curN +1;
        user.order = min(user.curN+1, user.BDFiA);
        ts.curT = ts.curT + ts.dt;
        ierr = PetscPrintf(comm,"curN = %d, deltaT = %g, curT=%g, Tenergy=%g\n\n", user.curN, (double)ts.dt,ts.curT, user.Tenergy[user.curN-1]/user.Length/user.Length); CHKERRQ(ierr);
        if (user.curN%ts.InterP==0)
        {
			ierr = GetFilename(user.curN,filenameascii);
            ierr = DataSaveASCII(user.x, filenameascii);
        }
    }while(ts.curT < ts.endT);
    PetscTime(&t2);
    
    PetscPrintf(PETSC_COMM_WORLD,"total time equals %f seconds.\n", t2-t1);
    char fileLastuASCII[PETSC_MAX_PATH_LEN] = "ASCIsol_0000";
    ierr = GetFilename(user.curN,fileLastuASCII);
    ierr = DataSaveASCII(user.x, fileLastuASCII);
    
    PetscInt i;
    FILE *fp;
    char fileTE[PETSC_MAX_PATH_LEN] = "BinaStepAndEnergy";
    
    if (!rank) {
        ierr = PetscFOpen(PETSC_COMM_SELF, fileTE, "w", &fp);
        for (i=0; i < user.curN; i++)
        {
            PetscFPrintf(PETSC_COMM_SELF, fp, "%g %g\n", user.stept[i], user.Tenergy[i]);
        }
        PetscFClose(PETSC_COMM_SELF, fp);
    }
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*ierr = KSPDestroy(&user.ksp); CHKERRQ(ierr);*/
	ierr = VecDestroy(&user.x); CHKERRQ(ierr);
	for (i_old = 0; i_old < user.BDFi; ++i_old){
		ierr =  VecDestroy(&user.x_olds[i_old]); CHKERRQ(ierr);
	}
    ierr = VecDestroy(&user.x_ex); CHKERRQ(ierr);
    ierr = VecDestroy(&user.x_sum); CHKERRQ(ierr);
    ierr = VecDestroy(&user.B); CHKERRQ(ierr);
    ierr = VecDestroy(&user.C); CHKERRQ(ierr);
    ierr = VecDestroy(&user.X); CHKERRQ(ierr);
    ierr = VecDestroy(&user.Y); CHKERRQ(ierr);
    ierr = MatDestroy(&user.A); CHKERRQ(ierr);
    ierr = DMDestroy(&user.da); CHKERRQ(ierr);
    
    ierr = PetscFinalize();
    return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitial"
/*
 FormInitialGuess - Forms initial values.
 
 Output Parameter:
 user.Xpre - vector
 */
PetscErrorCode FormInitial(DM da,AppCtx *user)
{
    PetscInt       i, j, Mx, My, xs, ys, xm, ym,use_random=0;
    PetscErrorCode ierr;
    PetscReal      L, hx, hy;
    PetscScalar    **x;
    
    ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE); CHKERRQ(ierr);
    
	if (use_random)
	{
	    PetscRandom rctx;
	    PetscReal sum, mean;
        PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
	    PetscRandomSetFromOptions(rctx);
	    PetscRandomSetInterval(rctx,-0.15,0.65);
	    VecSetRandom(user->x, rctx);
        VecSum(user->x, &sum);
	    mean=sum/(double)(Mx*My);
	    VecScale(user->x,0.25/mean);
	}
	else
	{
	    PetscReal lx,ly;
	    L = user->Length;
        hx = L/(PetscReal)Mx;
        hy = L/(PetscReal)My;

	    VecSet(user->x,0);
        ierr = DMDAVecGetArray(da, user->x, &x); CHKERRQ(ierr);
        ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
            
        for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				lx = i*hx;
				ly = j*hy;
		        x[j][i] = 0.05*sin(lx)*sin(ly);
			}
	    }
        ierr = DMDAVecRestoreArray(da, user->x, &x); CHKERRQ(ierr);
	}
    PetscFunctionReturn(0);
}
