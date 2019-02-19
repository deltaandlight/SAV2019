#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <math.h>
#include <petsctime.h>
#include "def.h"

#undef __FUNCT__
#define __FUNCT__ "Update"
/*
 * time stepping function
 */
PetscErrorCode Update(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	PetscErrorCode 	ierr;
  	
    for (user->BDFiA = user->BDFi - user->BDFAorB; user->BDFiA <= user->BDFi; user->BDFiA++ ){
    	user->order = min(user->curN+1, user->BDFiA);
    	ierr = PetscPrintf(PETSC_COMM_WORLD, "order = %d, discrete form = BDF%d%c\n", (int)user->order, (int)user->BDFi, (char)(user->BDFAorB+'A')); CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_WORLD, "calculate BDF%d%c\n", (int)user->BDFi, (char)(user->BDFAorB+'A')); CHKERRQ(ierr);
		if (user->BDFiA == user->BDFi - user->BDFAorB){
			ierr = VecSet(user->x_sum, 0); CHKERRQ(ierr);
			ierr = VecSet(user->x_ex, 0); CHKERRQ(ierr);
    		ierr = phi_ex(user); CHKERRQ(ierr);
    	}
    	//testout(user);
    	ierr = KSPCreate(PETSC_COMM_WORLD, &user->ksp); CHKERRQ(ierr);
		/*----KSP Setting！！！----------------------------------------------*/
		//PetscPrintf(PETSC_COMM_WORLD,"SetKSP in\n");
		SetKSP(user);
		//PetscPrintf(PETSC_COMM_WORLD,"SetKSP out\n");
    	/*-------------------------------------------------------------------------------------------*/
     	/*----B, C, X, b分配空间 (x已经分配空间了)----------------------------------------------*/
		ierr = VecDuplicate(user->x, &user->B); CHKERRQ(ierr);
    	ierr = VecDuplicate(user->x, &user->C); CHKERRQ(ierr);
    	ierr = VecDuplicate(user->x, &user->X); CHKERRQ(ierr);
    	ierr = VecDuplicate(user->x, &user->Y); CHKERRQ(ierr);
    	//PetscPrintf(PETSC_COMM_WORLD,"CalB\n");
		ierr = CalB(user); CHKERRQ(ierr);
		/*----compute time step------------------------*/
 		ierr = Getdt(user); CHKERRQ(ierr);
 		/*-------------------------------------------------------------------------------------------------*/
 		//PetscPrintf(PETSC_COMM_WORLD,"CalC\n");
    	ierr = CalC(user); CHKERRQ(ierr);
		/*----------------------------------------------------------------------------------------------------------------------------------------*/
		//PetscPrintf(PETSC_COMM_WORLD,"Calnewx\n");
    	ierr = Calnewx(user); CHKERRQ(ierr);
	 	/*----------------------------------------------------------------------------------------------------------------------------------------*/
		ierr = KSPDestroy(&user->ksp); CHKERRQ(ierr);
	}
    
  	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetKSP"
/*
 * 
 */
PetscErrorCode SetKSP(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	tsCtx          	*ts = user->ts;
  	PetscErrorCode 	ierr;
	PetscInt       	i, j, xe, ye;
	PetscScalar		dt, eps, beta, h;
	PetscReal		pardt,par00,par01,par02,par11;
	DMDALocalInfo	info;
	 
	dt = ts->dt;
	eps = user->eps;
	beta = user->beta;
	
  	ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
	h = user->Length/info.mx;
	xe = info.xs+info.xm;
  	ye = info.ys+info.ym;
  	switch((int)user->order){
  		case 1:
  			pardt = 1.0;
  			break;
  		case 2:
  			pardt = 2.0/3.0;
  			break;
		case 3:
  			pardt = 6.0/11.0;
  			break;
  		default:
  			pardt = 12.0/25.0;
	}

    par00 = 1.0+(4.0*beta/(eps*eps*h*h)+20.0/(h*h*h*h))*pardt*dt;
    par01 = (-beta/(eps*eps*h*h)-8.0/(h*h*h*h))*pardt*dt;
    par11 = 2.0/(h*h*h*h)*pardt*dt;
    par02 = 1.0/(h*h*h*h)*pardt*dt;
   	/*------------------------------------------*/
	for (j=info.ys; j<ye; j++) {
		for(i=info.xs;i<xe;i++){
			MatStencil	row = {0}, col[13] = {{0}};
			PetscScalar	v[13];
			PetscInt	ncols = 0;
			row.j = j;
			row.i = i;
			/*(0,0)点*/
			col[ncols].j=j;col[ncols].i=i;v[ncols++] = par00;
			/*(0,1)点*/
			col[ncols].j = j;   	col[ncols].i = i-1; v[ncols++] = par01;
			col[ncols].j = j;   	col[ncols].i = i+1; v[ncols++] = par01;
			col[ncols].j = j-1; 	col[ncols].i = i;   v[ncols++] = par01;
			col[ncols].j = j+1; 	col[ncols].i = i;   v[ncols++] = par01;
			/*(0,2)点*/
			col[ncols].j = j;   	col[ncols].i = i-2; v[ncols++] = par02;
			col[ncols].j = j;   	col[ncols].i = i+2; v[ncols++] = par02;
			col[ncols].j = j-2; 	col[ncols].i = i;   v[ncols++] = par02;
			col[ncols].j = j+2; 	col[ncols].i = i;   v[ncols++] = par02;
			/*(1,1)点*/
			col[ncols].j = j-1;   	col[ncols].i = i-1; v[ncols++] = par11;
			col[ncols].j = j+1;   	col[ncols].i = i+1; v[ncols++] = par11;
			col[ncols].j = j-1; 	col[ncols].i = i+1;   v[ncols++] = par11;
			col[ncols].j = j+1; 	col[ncols].i = i-1;   v[ncols++] = par11;
			ierr = MatSetValuesStencil(user->A,1, &row,ncols, col,v, INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(user->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(user->A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	/*----KSP的设置，注意其中的矩阵A是定值！！！----------------------------------------------*/
	ierr = KSPSetTolerances(user->ksp,1.e-2/((info.mx+1)*(info.my+1)),1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(user->ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(user->ksp, user->A, user->A); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalB"
/*
 * 
 */
PetscErrorCode CalB(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	PetscErrorCode 	ierr;
	PetscInt       	i, j, mx, my, xs, ys, xm, ym, xe, ye;
	PetscScalar		beta, eps, h;
	PetscReal		rtemp_nplus1_2, rtemp_n, Lapphi, GE1, GE1_nplus1_2, GE, GBxphi;/*temp*/
	Vec				localx;
	//PetcsInt		Ii,Jj;
	
  	ierr = DMDAGetInfo(user->da,0, &mx, &my,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
  	ierr = DMDAGetCorners(user->da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
	
	h = user->Length/mx;
	beta = user->beta;
	eps = user->eps;
	xe = xs+xm;
  	ye = ys+ym;
  	GE = 0;
  	GE1 = 0;
  	GE1_nplus1_2 = 0;
  	
    /*----point to x,B,----------------------------------------------*/
    PetscScalar **phi, **phi_ex, **aB, **phi_sum;/*phi=local_x*/
 
    ierr = DMGetLocalVector(user->da, &localx); CHKERRQ(ierr);
    
    ierr = DMGlobalToLocalBegin(user->da, user->x, INSERT_VALUES,localx); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, user->x, INSERT_VALUES,localx); CHKERRQ(ierr);

    
    ierr = DMDAVecGetArray(user->da,localx, &phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->x_ex, &phi_ex); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->x_sum, &phi_sum); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->B, &aB);

    for (j=ys; j<ye; j++) {
        for (i=xs; i<xe; i++) {
            rtemp_nplus1_2 = (phi_ex[j][i]*phi_ex[j][i]-1.0-beta);
            aB[j][i] 	   = rtemp_nplus1_2*phi_ex[j][i]/(eps*eps);	/*the B is U*/
            GE1_nplus1_2   = GE1_nplus1_2+h*h*rtemp_nplus1_2*rtemp_nplus1_2/(eps*eps*4.0);	/*local E1n+1/2*/
            
            rtemp_n  	   = (phi[j][i]*phi[j][i]-1.0-beta);
            Lapphi 		   = (phi[j+1][i]+phi[j-1][i]+phi[j][i+1]+phi[j][i-1]-4.0*phi[j][i])/(h*h);
            GE1 		   = GE1 + h*h*rtemp_n*rtemp_n/(eps*eps*4.0);	/*local E1*/
			GE  		   = GE + h*h*rtemp_n*rtemp_n/(eps*eps*4.0)
			+ h*h*phi[j][i]*(phi[j][i]*beta/(eps*eps)-Lapphi)/2.0;
            GBxphi		   = GBxphi + phi_sum[j][i]*aB[j][i]*h*h;
        }
    }
    ierr = MPI_Allreduce(&GE1_nplus1_2, &user->SE1_nplus1_2, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
 	ierr = MPI_Allreduce(&GE1, &user->SE1, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
 	ierr = MPI_Allreduce(&GE, &user->SE, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
 	ierr = MPI_Allreduce(&GBxphi, &user->SBxphi, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
 	
 	if (user->BDFiA == user->BDFi){
 		switch((int)user->order){
  			case 1:
  				user->Rold3	= sqrt(user->SE1+user->c0);
  				user->Rold2	= sqrt(user->SE1+user->c0);
  				user->Rold	= sqrt(user->SE1+user->c0);
 				user->R		= sqrt(user->SE1+user->c0);
  				break;
  			case 2:
  				user->Rold3	= sqrt(user->SE1+user->c0);
  				user->Rold2	= sqrt(user->SE1+user->c0);
  				user->Rold	= user->R;
  				user->R		= sqrt(user->SE1+user->c0);
  				break;
			case 3:
				user->Rold3	= sqrt(user->SE1+user->c0);
  				user->Rold2	= user->Rold;
  				user->Rold	= user->R;
  				user->R		= sqrt(user->SE1+user->c0);
  				break;
  			default:
  				user->Rold3	= user->Rold2;
  				user->Rold2	= user->Rold;
  				user->Rold	= user->R;
  				user->R		= sqrt(user->SE1+user->c0);
		}
		user->Tenergy[user->curN] = user->SE;
 		ierr = PetscPrintf(PETSC_COMM_WORLD,"R = %g,E = %g\n", user->R, user->SE); CHKERRQ(ierr);
	}
 	
 	user->SBxphi=user->SBxphi/sqrt(user->SE1_nplus1_2+user->c0);
 			
 	/*----这里才得到真正的B!!!------------------------------*/	
 	ierr = DMDAVecRestoreArray(user->da,localx, &phi); CHKERRQ(ierr);
 	ierr = DMDAVecRestoreArray(user->da, user->B, &aB); CHKERRQ(ierr);
 	ierr = DMDAVecRestoreArray(user->da, user->x_ex, &phi_ex); CHKERRQ(ierr);
 	ierr = DMDAVecRestoreArray(user->da, user->x_sum, &phi_sum); CHKERRQ(ierr);
 	ierr = DMRestoreLocalVector(user->da, &localx);
 	ierr = VecScale(user->B,1./sqrt(user->SE1_nplus1_2+user->c0)); CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "Getdt"
/*
 * 
 */
PetscErrorCode Getdt(void *ptr){
	AppCtx         	*user = (AppCtx*)ptr;
  	tsCtx          	*ts = user->ts;
  	PetscReal		ta,tb,norm,dG, GV, GW,dtl;
	PetscBool		flg = PETSC_TRUE;
	PetscErrorCode 	ierr;
 	if (user->curN == 0)
    {
        ts->dt = ts->tmin;
    }
    else{
    	dtl  = user->stept[user->curN - 1];
   		if (flg){
   			GV = user->SE;
     		GW = user->Tenergy[user->curN-1];
       		dG = fabs((GV-GW)/dtl);
        	ta = ts->tmin;
        	tb = ts->tmax/sqrt(1+ts->Zalpha*dG*dG);
        	ts->dt = max(ta,tb);
    	}
    	else{
    		ierr = VecAXPY(user->x_olds[0], -1.0, user->x_olds[1]); CHKERRQ(ierr);
        	ierr = VecNorm(user->x_olds[0], NORM_2, &norm); CHKERRQ(ierr);
        	PetscPrintf(PETSC_COMM_WORLD,"norm of two step solutions = %g\n",norm);
        	ta = ts->tmin;
        	dG = norm/dtl;
        	tb = ts->tmax/sqrt(1+ts->alpha*dG*dG);
        	ts->dt = max(ta,tb);
    	}
    }
    PetscFunctionReturn(0);
}
    
#undef __FUNCT__
#define __FUNCT__ "CalC"
/*
 * 
 */
PetscErrorCode CalC(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	tsCtx          	*ts = user->ts;
  	PetscErrorCode 	ierr;
	PetscInt       	i, j, mx, my, xs, ys, xm, ym, xe, ye;
	PetscReal		h,R_old[4],R_sum,pardt;
	PetscReal		Lgamma, Ggamma;/*temp*/
	Vec				localB, LapB,localX;
	
  	ierr = DMDAGetInfo(user->da,0, &mx, &my,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
  	ierr = DMDAGetCorners(user->da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
	h = user->Length/mx;
	xe = xs+xm;
  	ye = ys+ym;
  	
  	R_old[0] = user->R;
	R_old[1] = user->Rold;
	R_old[2] = user->Rold2;
	R_old[3] = user->Rold3;
  	switch((int)user->order){
  		case 1:
  			pardt = 1.0;
  			R_sum = BDFi(1, R_old);
  			break;
  		case 2:
  			pardt = 2.0/3.0;
  			R_sum = BDFi(2, R_old);
  			break;
  		case 3:
  			pardt = 6.0/11.0;
  			R_sum = BDFi(3, R_old);
  			break;
  		default:
  			pardt = 12.0/25.0;
  			R_sum = BDFi(4, R_old);
	}
  	
  	
    /*----point to x,B,C----------------------------------------------*/
    PetscScalar **phi_sum,**alocalB,**aC,**aLapB,**alocalX;/*phi=local_x*/ 
    ierr = DMGetGlobalVector(user->da, &LapB); CHKERRQ(ierr);
    
    ierr = DMGetLocalVector(user->da, &localB); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da, user->B, INSERT_VALUES,localB); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, user->B, INSERT_VALUES,localB); CHKERRQ(ierr);
	
    ierr = DMDAVecGetArray(user->da, LapB, &aLapB); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, localB, &alocalB); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->C, &aC); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->x_sum, &phi_sum); CHKERRQ(ierr);
	 
    
    for (j=ys; j<ye; j++) {
        for (i=xs; i<xe; i++) {           
            aLapB[j][i] = (alocalB[j+1][i]+alocalB[j-1][i]+alocalB[j][i+1]+alocalB[j][i-1]
			-4.0*alocalB[j][i])/(h*h);
            aC[j][i] = phi_sum[j][i]
			+ pardt*ts->dt*aLapB[j][i]*(R_sum-1.0/2.0*user->SBxphi);
        }
    }
    ierr = DMDAVecRestoreArray(user->da, LapB, &aLapB); CHKERRQ(ierr);
    
    ierr = KSPSolve(user->ksp, LapB, user->X); CHKERRQ(ierr);	/*得到的是一个中间参数!!!A-1Gbn,用来算gamma的，暂时用一下X*/
    
	ierr = DMGetLocalVector(user->da, &localX); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->da, user->X, INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, user->X, INSERT_VALUES,localX); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->da, localX, &alocalX); CHKERRQ(ierr);
	for (j=ys; j<ye; j++) {
        for (i=xs; i<xe; i++) {
            Lgamma = alocalB[j][i]*alocalX[j][i];/*这里没加负号*/
            Ggamma = Ggamma + Lgamma*h*h;
        } 
    }	
    		/*----计算 b_n*phi_(n+1)----------------------------*/ 
    ierr = MPI_Allreduce(&Ggamma, &user->gamma, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD);
    ierr = DMDAVecRestoreArray(user->da, localX, &alocalX); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(user->da, &localX);
    ierr = VecDestroy(&LapB); CHKERRQ(ierr);
    
 	ierr = DMDAVecRestoreArray(user->da, localB, &alocalB); CHKERRQ(ierr);
 	ierr = DMDAVecRestoreArray(user->da, user->C, &aC); CHKERRQ(ierr);
 	ierr = DMDAVecRestoreArray(user->da, user->x_sum, &phi_sum); CHKERRQ(ierr);
 	ierr = DMRestoreLocalVector(user->da, &localB);
 	
	
    ierr = KSPSolve(user->ksp, user->C, user->Y); CHKERRQ(ierr);	/*得到的是一个中间参数X*/
		   
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Calnewx"
/*
 * 
 */
PetscErrorCode Calnewx(void *ptr){
  	AppCtx         	*user = (AppCtx*)ptr;
  	tsCtx          	*ts = user->ts;
  	PetscErrorCode 	ierr;
	PetscInt       	i, j, mx, my, xs, ys, xm, ym, xe, ye, i_old;
	PetscReal		h,pardt;
	PetscReal		LBxphiP, GBxphiP;/*temp*/
	Vec				localB, localY;
	
  	ierr = DMDAGetInfo(user->da, 0, &mx, &my, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
  	ierr = DMDAGetCorners(user->da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
	
	h=user->Length/mx;
	xe = xs+xm;
  	ye = ys+ym;
  	switch((int)user->order){
  		case 1:
  			pardt = 1.0;
  			break;
  		case 2:
  			pardt = 2.0/3.0;
  			break;
  		case 3:
  			pardt = 6.0/11.0;
  			break;
  		default:
  			pardt = 12.0/25.0;
	}
  	
    /*----point to x,B,C----------------------------------------------*/
    PetscScalar **alocalB,**alocalY;/*phi=local_x*/ 
     
    ierr = DMGetLocalVector(user->da, &localY); CHKERRQ(ierr);
	ierr = VecDuplicate(localY, &localB); CHKERRQ(ierr);
	
    ierr = DMGlobalToLocalBegin(user->da, user->B, INSERT_VALUES,localB); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, user->B, INSERT_VALUES,localB); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da, user->Y, INSERT_VALUES,localY); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, user->Y, INSERT_VALUES,localY); CHKERRQ(ierr);
	
    ierr = DMDAVecGetArray(user->da,localB, &alocalB); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da,localY, &alocalY); CHKERRQ(ierr);
	 
	GBxphiP = 0.0;	/*初始化 b_n*phi_(n+1) */ 
    for (j=ys; j<ye; j++) {
        for (i=xs; i<xe; i++) {
            LBxphiP = alocalB[j][i]*alocalY[j][i]/(1-ts->dt*pardt*user->gamma/2.0);
            GBxphiP = GBxphiP + LBxphiP*h*h;
        } 
    }	
    		/*----计算 b_n*phi_(n+1)----------------------------*/ 
    ierr = MPI_Allreduce(&GBxphiP, &user->SBxphiP, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD);
    		
    		/*----计算最后的方程------------------*/ 
	ierr = VecAXPY(user->Y, ts->dt*pardt*user->SBxphiP/2.0, user->X); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->da,localB, &alocalB); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da,localY, &alocalY); CHKERRQ(ierr);
 	ierr = DMRestoreLocalVector(user->da, &localB);
 	ierr = DMRestoreLocalVector(user->da, &localY);
 	
	if (user->BDFiA == user->BDFi){
		for (i_old = 0; i_old < user->BDFi-1; ++i_old){
			ierr = VecCopy(user->x_olds[i_old], user->x_olds[i_old+1]); CHKERRQ(ierr);
		}
		ierr = VecCopy(user->Y, user->x_olds[0]); CHKERRQ(ierr);
		ierr = VecCopy(user->Y, user->x); CHKERRQ(ierr);
	}
	else {
		ierr = VecCopy(user->Y, user->x_ex); CHKERRQ(ierr);
	}
    
    	 /*----------------------------------------------------------------------------------------------------------------------------------------*/
      /*----------------------------------------------------------------------------------------------------------------------------------------*/
	PetscFunctionReturn(0);
}


