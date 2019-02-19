#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <math.h>
#include <petsctime.h>
#include "def.h"


#undef __FUNCT__
#define __FUNCT__ "DataSaveASCII"

PetscErrorCode DataSaveASCII(Vec x, char *filename)
{
    PetscErrorCode ierr;
    PetscViewer    dataviewer;

    PetscFunctionBegin;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename, &dataviewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(dataviewer,PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);
    ierr = VecView(x,dataviewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&dataviewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataLoadASCII"
PetscErrorCode DataLoadASCII(Vec x, char *filename)
{
    PetscErrorCode ierr;
    PetscViewer    dataviewer;

    PetscFunctionBegin;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename, &dataviewer); CHKERRQ(ierr);
    ierr = VecLoad(x,dataviewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&dataviewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GetFilename"
PetscErrorCode GetFilename(PetscInt i, char *filename)
{
    PetscInt i1,i2,i3,i4;

    PetscFunctionBegin;
    i1 = i%10; i/=10;
    i2 = i%10; i/=10;
    i3 = i%10; i/=10;
    i4 = i%10;
   /* filename[9] = '0' + i4;
    filename[10] = '0' + i3;
    filename[11] = '0' + i2;
    filename[12] = '0' + i1;*/
 /* filename[21] = '0' + i4;
    filename[22] = '0' + i3;
    filename[23] = '0' + i2;
    filename[24] = '0' + i1;*/

    filename[8] = '0' + i4;
    filename[9] = '0' + i3;
    filename[10] = '0' + i2;
    filename[11] = '0' + i1;
    PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "max"
PetscScalar max(PetscScalar a, PetscScalar b)
{
    if (a>b){return a;}
    else{return b;}
}

#undef __FUNCT__
#define __FUNCT__ "min"
PetscScalar min(PetscScalar a, PetscScalar b)
{
    if (a>b){return b;}
    else{return a;}
}

#undef __FUNCT__
#define __FUNCT__ "testout"
PetscErrorCode testout(void *ptr){
	AppCtx *user = (AppCtx*)ptr;
	tsCtx  *ts = user->ts;
	PetscErrorCode 	ierr;
	char  filenamesum[PETSC_MAX_PATH_LEN] = "0000sum_0000";
	char  filenameex[PETSC_MAX_PATH_LEN]  = "00000ex_0000";
    if (user->curN%ts->InterP==0 && user->curN/ts->InterP < 5)
    {
		ierr = GetFilename(user->curN,filenamesum); CHKERRQ(ierr);
        ierr = DataSaveASCII(user->x, filenamesum); CHKERRQ(ierr);
        ierr = GetFilename(user->curN,filenameex); CHKERRQ(ierr);
        ierr = DataSaveASCII(user->x, filenameex); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}
