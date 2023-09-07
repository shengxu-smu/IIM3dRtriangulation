#################################################################
#  October 2020
###################################################################
F90 = gfortran
FFLAGS   = -O3

OBJ =   iim3d.o run_input.o mesh.o allocate_field.o loadobject.o\
	allocate_surface.o  coord4intersection.o \
	dotproduct.o crossproduct.o pointinside.o raycrossing.o \
	surface_initialize.o field_initial.o velocity_reset.o \
	cfl.o rk4.o JC4uvwp.o allocate_Lagrange.o resetLagrange.o \
	dudnPJC.o interpolate.o  dudnnPJC.o cjc4d2u.o CG4p.o\
	c2_solver.o JCinterpolate.o JCinterpolate2.o IrrP_JC2.o \
	IrrP_pressureJC2.o jc_velocity.o correction_interpolate.o \
	correction_strain.o u_interpolate.o u_strain.o udu_surface.o \
	v_interpolate.o v_strain.o vdv_surface.o w_interpolate.o \
	w_strain.o wdw_surface.o jc_pressure.o correction_difference.o \
	old_save.o singular_call.o initial.o allocate_vfft.o \
	pressure.o pbc.o field_interpolate.o poisson_fft.o rhs.o \
	surface_move.o \
	deallocate_Lagrange.o ubc.o vbc.o wbc.o \
	animation.o vorticity.o vffpack.o\
	data_read.o data_write.o surface_plot.o field_plot.o\
	run_output.o \
	CartJump1st4u.o CartJump2nd4u.o CartJump2nd4p.o rhs4p.o\
	dotproduct2.o time_output.o Force4Panel.o

MOD1 = para.o
MOD2 = field.o
MOD3 = Lagrange.o
MOD4 = para.o field.o
MOD5 = para.o Lagrange.o
MOD6 = para.o field.o Lagrange.o
MOD7 = para.o field.o Lagrange.o vfft.o
MOD8 = para.o vfft.o
MOD9 = para.o field.o vfft.o
iim3d.exe : $(OBJ) $(MOD7) 
	$(F90) $(FFLAGS) -o $@ $^

para.o : para.f90
	$(F90) -c para.f90
field.o : field.f90
	$(F90) -c field.f90
Lagrange.o : Lagrange.f90
	$(F90) -c Lagrange.f90
vfft.o : vfft.f90
	$(F90) -c vfft.f90
vffpack.o : vffpack.f
	$(F90) -c vffpack.f
CartJump1st4u.o : CartJump1st4u.f90
	$(F90) -c CartJump1st4u.f90
CartJump2nd4u.o : CartJump2nd4u.f90
	$(F90) -c CartJump2nd4u.f90
CartJump2nd4p.o : CartJump2nd4p.f90
	$(F90) -c CartJump2nd4p.f90
rhs4p.o : rhs4p.f90
	$(F90) -c rhs4p.f90
run_input.o : run_input.f90 $(MOD1)
	$(F90) -c run_input.f90
allocate_field.o : allocate_field.f90 $(MOD4) 
	$(F90) -c allocate_field.f90
mesh.o : mesh.f90 $(MOD4)
	$(F90) -c mesh.f90
loadobject.o : loadobject.f90 $(MOD5)
	$(F90) -c loadobject.f90
allocate_surface.o : allocate_surface.f90 $(MOD3)
	$(F90) -c allocate_surface.f90
dotproduct2.o : dotproduct2.f90
	$(F90) -c dotproduct2.f90
dotproduct.o : dotproduct.f90
	$(F90) -c dotproduct.f90
crossproduct.o : crossproduct.f90
	$(F90) -c crossproduct.f90
pointinside.o : pointinside.f90
	$(F90) -c pointinside.f90
coord4intersection.o : coord4intersection.f90
	$(F90) -c coord4intersection.f90
raycrossing.o : raycrossing.f90 $(MOD6) 
	$(F90) -c raycrossing.f90

surface_initialize.o : surface_initialize.f90 $(MOD6)
	$(F90) -c surface_initialize.f90
field_initial.o : field_initial.f90 $(MOD4)
	$(F90) -c field_initial.f90
velocity_reset.o: velocity_reset.f90 $(MOD6) 
	$(F90) -c velocity_reset.f90
cfl.o : cfl.f90 $(MOD4) 
	$(F90) -c cfl.f90
rk4.o : rk4.f90 $(MOD6) 
	$(F90) -c rk4.f90
JC4uvwp.o : JC4uvwp.f90 $(MOD6) 
	$(F90) -c JC4uvwp.f90
allocate_Lagrange.o : allocate_Lagrange.f90 $(MOD3)
	$(F90) -c allocate_Lagrange.f90
resetLagrange.o : resetLagrange.f90 $(MOD3)
	$(F90) -c resetLagrange.f90
dudnPJC.o : dudnPJC.f90 $(MOD6) 
	$(F90) -c dudnPJC.f90
dudnnPJC.o : dudnnPJC.f90 $(MOD6) 
	$(F90) -c dudnnPJC.f90
interpolate.o : interpolate.f90
	$(F90) -c interpolate.f90
cjc4d2u.o : cjc4d2u.f90
	$(F90) -c cjc4d2u.f90
CG4p.o : CG4p.f90 $(MOD3)
	$(F90) -c CG4p.f90
c2_solver.o : c2_solver.f90
	$(F90) -c c2_solver.f90
JCinterpolate.o : JCinterpolate.f90
	$(F90) -c JCinterpolate.f90
JCinterpolate2.o : JCinterpolate2.f90
	$(F90) -c JCinterpolate2.f90
IrrP_JC.o : IrrP_JC.f90 $(MOD6) 
	$(F90) -c IrrP_JC.f90
IrrP_JC2.o : IrrP_JC2.f90 $(MOD6) 
	$(F90) -c IrrP_JC2.f90
IrrP_pressureJC.o : IrrP_pressureJC.f90 $(MOD6)
	$(F90) -c IrrP_pressureJC.f90
IrrP_pressureJC2.o : IrrP_pressureJC2.f90 $(MOD6)
	 $(F90) -c IrrP_pressureJC2.f90
jc_velocity.o : jc_velocity.f90 $(MOD5)
	$(F90) -c jc_velocity.f90
correction_interpolate.o : correction_interpolate.f90 $(MOD6)
	$(F90) -c correction_interpolate.f90
correction_strain.o : correction_strain.f90 $(MOD6) 
	$(F90) -c correction_strain.f90
u_interpolate.o : u_interpolate.f90 $(MOD6) 
	$(F90) -c u_interpolate.f90
u_strain.o : u_strain.f90 $(MOD6) 
	$(F90) -c u_strain.f90
udu_surface.o : udu_surface.f90 $(MOD6) 
	$(F90) -c udu_surface.f90
v_interpolate.o : v_interpolate.f90 $(MOD6) 
	$(F90) -c v_interpolate.f90
v_strain.o : v_strain.f90 $(MOD6)
	$(F90) -c v_strain.f90
vdv_surface.o : vdv_surface.f90 $(MOD6) 
	$(F90) -c vdv_surface.f90
w_interpolate.o : w_interpolate.f90 $(MOD6)
	$(F90) -c w_interpolate.f90
w_strain.o : w_strain.f90 $(MOD6)
	$(F90) -c w_strain.f90
wdw_surface.o : wdw_surface.f90 $(MOD6) 
	$(F90) -c wdw_surface.f90
jc_pressure.o : jc_pressure.f90 $(MOD5) 
	$(F90) -c jc_pressure.f90
correction_difference.o : correction_difference.f90 $(MOD6)
	$(F90) -c correction_difference.f90
old_save.o : old_save.f90 $(MOD6)
	$(F90) -c old_save.f90 
singular_call.o : singular_call.f90 $(MOD6)
	$(F90) -c singular_call.f90
initial.o : initial.f90 $(MOD7) 
	$(F90) -c initial.f90
allocate_vfft.o : allocate_vfft.f90 $(MOD8)
	$(F90) -c allocate_vfft.f90
pressure.o : pressure.f90 $(MOD7)
	$(F90) -c pressure.f90
pbc.o : pbc.f90 $(MOD4)
	$(F90) -c pbc.f90
field_interpolate.o : field_interpolate.f90 $(MOD6)
	$(F90) -c field_interpolate.f90
poisson_fft.o : poisson_fft.f90 $(MOD9)
	$(F90) -c poisson_fft.f90
rhs.o : rhs.f90 $(MOD6)
	$(F90) -c rhs.f90
surface_move.o : surface_move.f90 $(MOD5)
	$(F90) -c surface_move.f90
deallocate_Lagrange.o : deallocate_Lagrange.f90 $(MOD3)
	$(F90) -c deallocate_Lagrange.f90
ubc.o : ubc.f90 $(MOD4)
	$(F90) -c ubc.f90
vbc.o : vbc.f90 $(MOD4)
	$(F90) -c vbc.f90
wbc.o : wbc.f90 $(MOD4)
	$(F90) -c wbc.f90
animation.o : animation.f90 $(MOD6)
	$(F90) -c animation.f90
vorticity.o : vorticity.f90 $(MOD4)
	$(F90) -c vorticity.f90
data_read.o : data_read.f90 $(MOD6)
	$(F90) -c data_read.f90
data_write.o : data_write.f90 $(MOD6)
	$(F90) -c data_write.f90
surface_plot.o : surface_plot.f90 $(MOD5)
	$(F90) -c surface_plot.f90
field_plot.o : field_plot.f90 $(MOD4)
	$(F90) -c field_plot.f90
run_output.o : run_output.f90 $(MOD1)
	$(F90) -c run_output.f90
time_output.o : time_output.f90 $(MOD6)
	$(F90) -c time_output.f90
Force4Panel.o : Force4Panel.f90 $(MOD6)
	$(F90) -c Force4Panel.f90
iim3d.o : iim3d.f90 $(MOD7)  
	$(F90) -c iim3d.f90

all : iim3d.exe 

clean:
	rm -f *~
	rm -f *.mod
	rm -f *.o
	rm *.exe

####### End of Makefile #######


