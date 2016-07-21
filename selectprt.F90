! initial and final time and timestep can be modified
program selectprt
  !select the particles at a given time which fulfil some requirement on energy
  
  logical exst, statusfile, twod
  character fnametst*22, fnamestep*15
  integer rank,numfiles,num
  integer(kind=4) lap,num1,proc,prevlap,lapsel,subpart,nchild,powcut
  real*8 x
  real y,z,u,v,w,gamma,bx0,by0,bz0,ex0,ey0,ez0
  integer nseli,ntoti,nsele,ntote
  real, dimension (4) :: gami_cut, game_cut
  integer, dimension (4) :: whcut
  
  twod=.false. !.true. for 2d.sig0
  subpart=50 !downsampling, for highest energies
  nchild=4 !if>1, then save more particles at higher energies (3 or 4)
  gami_cut=(/50.,100.,200.,400./)!(/40.,80.,160.,320./)
  game_cut=gami_cut         !(/40.,80.,160.,320./)      
  lapsel=100000
  
  whcut=0
  ntoti=0
  nseli=0
  ntote=0
  nsele=0
  
  print *, "gami_cut:", gami_cut
  print *, "game_cut:", game_cut
  
  rank=0
  exst=.true.
  do while (exst) 
     write(fnametst,"(a9,i3.3)") "test.prt.",rank
     inquire(file=fnametst,exist=exst)
     rank=rank+1
  enddo
  
  print *, "Detected ",rank, " files"
  
  !      goto 607 
  
  !     write the test.spl. file at the selected time
  numfiles=rank
  
  prevlap=0
  
  do num=0,numfiles-1
     exst=.false.
     write(fnametst,"(a9,i3.3)") "test.prt.",num
     
     open(unit=3, file=fnametst,form='formatted')
     print *, "opening ",fnametst
     
301  continue
     if (twod) then
        read(3,*,err=398,end=399)lap,num1,proc,x,y,u,v,gamma, &
             bz0,ex0,ey0
     else
        read(3,*,err=398,end=399)lap,num1,proc,x,y,z,u,v,w,gamma, &
             bx0,by0,bz0,ex0,ey0,ez0
     endif
     
     if(lap .eq. lapsel) then
        
        write(fnamestep,"(a9,i6.6)") "test.spl.",lap
        
        inquire(file=fnamestep,exist=exst)
        
        if(lap .ne. prevlap) then 
           print *, fnametst, lap
           if(exst) close(22)
           
           if(num .eq. 0) then
              open(unit=22,file=fnamestep,form='formatted')
              exst=.true.
           else
              inquire(file=fnamestep,exist=exst)
              if(exst) then 
                 open(unit=22,file=fnamestep,status='old' &
                      ,position='APPEND')
              else
                 open(unit=22,file=fnamestep,form='formatted')
                 exst=.true.
              endif
           endif
        endif
        
        if(exst) then
           if (twod) then 
              write(22,fmt="(i7, i12,i5,2(F10.3,' '),6(E15.5,' ') &
                   )")lap,num1,proc,x,y,u,v,gamma,bz0 &
                   ,ex0,ey0
           else
              write(22,fmt="(i7, i12,i5,3(F10.3,' '),10(E15.5, &
                   ' '))")lap,num1,proc,x,y,z,u,v,w,gamma, &
                   bx0,by0,bz0,ex0,ey0,ez0                    
           endif
        endif
        
        prevlap=lap
     endif
     
     goto 301
398  print *,"error reading the test.prt. file"
399  continue
     
  enddo !do num=0,numfiles-1
  close(22)
  
607 continue
  
  !     write the selected particles in test.sel.
  
  write(fnamestep,"(a9,i6.6)") "test.sel.",lapsel
  open(unit=23,file=fnamestep)
  close(23,status='delete') 
  open(unit=23,file=fnamestep)
  
  
  write(fnametst,"(a9,i6.6)") "test.spl.",lapsel
  open(unit=4, file=fnametst,form='formatted')
  print *, "opening ",fnametst
  
401 continue
  if (twod) then
     read(4,*,err=498,end=499)lap,num1,proc,x,y,u,v,gamma, &
          bz0,ex0,ey0
  else
     read(4,*,err=498,end=499)lap,num1,proc,x,y,z,u,v,w,gamma, &
          bx0,by0,bz0,ex0,ey0,ez0
  endif
  
  
  write(fnamestep,"(a9,i6.6)") "test.sel.",lapsel
  
  inquire(file=fnamestep,exist=exst) 
  
  statusfile=0
  if(exst) statusfile=1
  if(statusfile .eq. 0) then 
     open(unit=23,file=fnamestep,form='formatted')
  else
     open(unit=23,file=fnamestep,status='old', &
          position='APPEND')
  endif
  
  if (modulo(abs(num1)-1,2) .eq. 0) then !ions
     ntoti=ntoti+1
     where(gami_cut/gamma .gt. 1.) whcut=1
     if (whcut(1) .eq. 0) then !above gami_cut(1)
        powcut=sum(whcut)
        if (modulo(abs(num1),subpart*nchild**powcut) .eq. 1) then
           nseli=nseli+1
           write(23,fmt="(i7, i12,i5,F10.3)" &
                )lap,num1,proc,gamma
        endif
     endif
     whcut=0
  endif
  if (modulo(abs(num1),2) .eq. 0) then !lecs
     ntote=ntote+1
     where(game_cut/gamma .gt. 1.) whcut=1
     if (whcut(1) .eq. 0) then !above game_cut(1)
        powcut=sum(whcut)
        if (modulo(abs(num1),subpart*nchild**powcut) .eq. 0) then
           nsele=nsele+1
           write(23,fmt="(i7, i12,i5,F10.3)" &
                )lap,num1,proc,gamma
        endif
     endif
     whcut=0
  endif
  
  goto 401
498 print *,"error reading the test.spl. file"
499 continue         
  
  print*,"number of selected ions:",nseli
  print*,"fraction of selected ions (%):", &
       100.*real(nseli)/real(ntoti)
  print*,"number of selected lecs:",nsele
  print*,"fraction of selected lecs (%):", &
       100.*real(nsele)/real(ntote)
end program selectprt
