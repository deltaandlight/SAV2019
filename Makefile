CFLAGS    =
FFLAGS  	= 
CPPFLAGS  =
FPPFLAGS  =
 
export PETSC_DIR=/home/yangchao/zl/petsc-3.6.4
export PETSC_ARCH=linux-gnu-c-debug
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
 
objects = main.o BDF.o file_func.o iteration.o
 
main: $(objects) chkopts
	-${CLINKER} -o main -g $(objects) ${PETSC_LIB}
	${RM} $(objects)

$(objects): def.h

.PHONY: run run_serial mvdata cleandata cpmain

n_process = 4
n_grid = 1024
BDFi = 3
BDFAorB = 0
beta = 0.0
eps = 0.1
endT = 0.032
dt = 1e-3
interP = 1
ksp_setting = -ksp_monitor -ksp_converged_reason\
							-ksp_atol 1.e-13 -ksp_rtol 1.e-6\
							-ksp_type gmres -ksp_gmres_restart 30\
							-ksp_pc_side right -pc_type asm -pc_asm_type restrict -pc_asm_overlap 2\
							-sub_ksp_type preonly -sub_pc_type lu
run:
	module add mpich/3.2.1
	mpiexec -n $(n_process) ./main $(ksp_setting)\
		-da_grid_x $(n_grid) -da_grid_y $(n_grid)\
		-BDFi $(BDFi) -BDFAorB $(BDFAorB)\
		-beta $(beta) -eps $(eps)\
		-endT $(endT) -Tmin $(dt) -Tmax $(dt) -interP $(interP)

run_serial:
	module add mpich/3.2.1
	mpiexec -n 1 ./main $(ksp_setting)\
		-da_grid_x $(n_grid) -da_grid_y $(n_grid)\
		-BDFi $(BDFi) -BDFAorB $(BDFAorB)\
		-beta $(beta) -eps $(eps)\
		-endT $(endT) -Tmin $(dt) -Tmax $(dt) -interP $(interP)

mvdata:
	-mv ASCIsol_* ~/zl/src/SAV/2019/data/BDF/$(dt)/$(n_grid) -f
	-mv BinaStepAndEnergy ~/zl/src/SAV/2019/data/BDF/$(dt)/$(n_grid) -f

cleandata:
	-rm *_0*

cpmain:
	cp main ./3A/main
	cp main ./3B/main
	cp main ./4A/main
	cp main ./4B/main
