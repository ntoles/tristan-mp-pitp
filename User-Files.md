## User file

The user file, `user_*.F90`, such as `user_shock.F90` or `user_weibel.F90`, is used for readily expanding the behavior of the code without having to alter any of the code-flow logic. The user defines routines for loading and injecting particles here, as well special boundary conditions on fields and particles (e.g. reflecting walls) and external fields.

```fortran
subroutine get_external_fields
```

This subroutine compute external magnetic fields to be added to the particle mover. These fields do not evolve via Maxwell equations. Note that it uses local coordinates corresponding to different CPU, so one should convert that to global coordinates first.

**Example:** This example sets an external `B_z` field, that is an `bz_ext0*sin(x)*sin(y)` function of global coordinate x and y.

```fortran
subroutine get_external_fields
  real,intent(inout):: bx_ext,by_ext,bz_ext, ex_ext, ey_ext, ez_ext
  real, intent(in):: x,y,z
  ex_ext=0.
  ey_ext=0.
  ez_ext=0.
  bx_ext=0.
  by_ext=0.
  bz_ext=bz_ext0*sin(x)*sin(y+modulo(rank,sizey)*(myall-5)) 
  !bz_0 is defined in the input reading part above, 
  !note the coordinate transformation rule.
end subroutine get_external_fields
```


```fortran
subroutine init_EMfields_user()
```
This subroutine sets the electromagnetic fields of any specific user purpose.

**Example:** This example sets an initial jump for in-plane magnetic field by. Also component perpendicular to interface surface `bx` is added.

```fortran
subroutine init_EMfields_user()
  integer :: i, j, k, jglob, kglob
  !initialize B field to be set by Binit and the inclination angle
  !-- used for shock problems
  Binit=sqrt((gamma0)*ppc0*.5*c2*(mi+me)*sigma) 
  do k=1,mz
    do j=1,my
      do i=1,mx
        jglob=j+modulo(rank,sizey)*(myall-5) 
        !global j,k coords are defined from local ones
        kglob=k+(rank/sizey)*(mzall-5)
        !magnetic fields components are initialized
        bx(i,j,k)=0.5*Binit 
        bz(i,j,k)=0.*Binit*sin(btheta)*sin(bphi)
        by(i,j,k)=2.*Binit*tanh((i-mx/2)/4.)
        ex(i,j,k)=0. !electric field is initialized
        ey(i,j,k)=0.
        ez(i,j,k)=-0.
      enddo
    enddo
  enddo
end subroutine init_EMfields_user
```