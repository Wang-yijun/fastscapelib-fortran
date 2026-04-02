subroutine Fastscape_Named_VTK (additional_outputs, n_additional_outputs, ids, vex, istep, foldername, k)

    ! subroutine to create a simple VTK file for plotting inside of aspect folder.
    ! To use, add this file into the fastscape src directory, and then add
    ! the line ${FASTSCAPELIB_SRC_DIR}/Named_VTK.f90 to the appropriate place
    ! in cmakelists.txt

    use FastScapeContext
    implicit none

    integer, intent(in) :: k, istep !,atime
    double precision, intent(in) :: vex
    double precision, intent(in), dimension(*) :: additional_outputs
    integer, intent(in) :: n_additional_outputs
    integer, intent(in), dimension(*) :: ids
    character(len=k), intent(in) :: foldername
    character(len=7) cstep

    integer nheader,nfooter,npart1,npart2, ftime
    character header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12
    integer i,j
    character(len=10) fwtime
    double precision dx,dy

    integer :: o,ip
    character(len=64) :: varname

    dx = xl/(nx - 1)
    dy = yl/(ny - 1)

    ! Finding the fastscape time incremented by dt. dt may be a double
    ! that is converted to an integer, but we don't care about then
    ! exact time that much.
    ! ftime = (atime + step*dt)

    !Writing the fastscape time to a string, up to 1 billion years. (add another zero if needed eventually)
   ! write (fwtime,'(i10)') ftime
   ! if (ftime.lt.10) fwtime(1:9)='000000000'
   ! if (ftime.lt.100) fwtime(1:8)='00000000'
   ! if (ftime.lt.1000) fwtime(1:7)='0000000'
   ! if (ftime.lt.10000) fwtime(1:6)='000000'
   ! if (ftime.lt.100000) fwtime(1:5)='00000'
   ! if (ftime.lt.1000000) fwtime(1:4)='0000'
   ! if (ftime.lt.10000000) fwtime(1:3)='000'
   ! if (ftime.lt.100000000) fwtime(1:2)='00'
   ! if (ftime.lt.1000000000) fwtime(1:1)='0'

    write (cstep,'(i7)') istep
    if (istep.lt.10) cstep(1:6)='000000'
    if (istep.lt.100) cstep(1:5)='00000'
    if (istep.lt.1000) cstep(1:4)='0000'
    if (istep.lt.10000) cstep(1:3)='000'
    if (istep.lt.100000) cstep(1:2)='00'
    if (istep.lt.1000000) cstep(1:1)='0'

    ! Use the output folder of ASPECT and create a VTK folder inside of it.
    call system ('mkdir -p '//trim(foldername)//'/VTK')

    write (nxc,'(i6)') nx
    write (nyc,'(i6)') ny
    write (nnc,'(i12)') nn

    header(1:1024) = ''
    header = '# vtk DataFile Version 3.0'//char(10)//'FastScape'//char(10) &
          //'BINARY'//char(10)//'DATASET STRUCTURED_GRID'//char(10) &
          //'DIMENSIONS '//nxc//' '//nyc//' 1'//char(10)//'POINTS' &
          //nnc//' float'//char(10)
    nheader = len_trim(header)

    footer(1:1024) = ''
    footer = 'POINT_DATA'//nnc//char(10)
    nfooter = len_trim(footer)

    part1(1:1024) = ''
    part1 = 'SCALARS '
    npart1 = len_trim(part1) + 1

    part2(1:1024) = ''
    part2 = ' float 1'//char(10)//'LOOKUP_TABLE default'//char(10)
    npart2 = len_trim(part2)

    ! ----------------------------------------------------------------------
    ! Topography VTK (now written as stream, so we do NOT need recl anymore)
    ! ----------------------------------------------------------------------
    open(unit=77, file=(foldername//'/VTK/Topography'//cstep//'.vtk'), status='unknown', &
        form='unformatted', access='stream')

    ! Header + grid points
    write(77) header(1:nheader)
    write(77) ((sngl(dx*(i-1)), sngl(dy*(j-1)), sngl(h(i+(j-1)*nx)*abs(vex)), i=1,nx), j=1,ny)
    write(77) footer(1:nfooter)

    ! Fixed scalar fields (same as before)
    write(77) part1(1:npart1)//'topography'//part2(1:npart2)
    write(77) sngl(h(1:nn))

    write(77) part1(1:npart1)//'basement'//part2(1:npart2)
    write(77) sngl(b(1:nn))

    write(77) part1(1:npart1)//'erosion_rate'//part2(1:npart2)
    write(77) sngl(erate(1:nn))

    write(77) part1(1:npart1)//'total_erosion'//part2(1:npart2)
    write(77) sngl(etot(1:nn))

    write(77) part1(1:npart1)//'drainage_area'//part2(1:npart2)
    write(77) sngl(a(1:nn))

    write(77) part1(1:npart1)//'catchment'//part2(1:npart2)
    write(77) sngl(catch(1:nn))

    ! ----------------------------------------------------------------------
    ! Additional outputs: written in the order provided by C++,
    ! but named by ids(o) so the meaning is explicit.
    !
    ! Layout expected from C++:
    ! additional_outputs( (o-1)*nn + p )  with o=1..n_additional_outputs, p=1..nn
    ! ----------------------------------------------------------------------
    if (n_additional_outputs .gt. 0) then
      do o = 1, n_additional_outputs

        varname = 'additional_output'

        select case (ids(o))

          case (0)
            varname = 'bedrock_river_incision_rate'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(additional_outputs((o-1)*nn + ip)), ip=1, nn )

          case (1)
            varname = 'bedrock_transport_coefficient'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(additional_outputs((o-1)*nn + ip)), ip=1, nn )

          case (2)
            varname = 'marine_sand_transport_coefficient'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(kdsea1(ip)), ip=1, nn )

          case (3)
            varname = 'marine_silt_transport_coefficient'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(kdsea2(ip)), ip=1, nn )

          case (4)
            varname = 'uplift_rate'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(additional_outputs((o-1)*nn + ip)), ip=1, nn )

          case default
            varname = 'additional_output'

            write(77) part1(1:npart1)//trim(varname)//part2(1:npart2)
            write(77) ( sngl(additional_outputs((o-1)*nn + ip)), ip=1, nn )

        end select

      end do
    end if

    close(77)

    ! ----------------------------------------------------------------------
    ! Optional extra VTK files if vex < 0 (Basement + SeaLevel)
    ! These are also written as stream now.
    ! ----------------------------------------------------------------------
    if (vex .lt. 0.d0) then

      open(unit=77, file=(foldername//'/VTK/Basement'//cstep//'.vtk'), status='unknown', &
          form='unformatted', access='stream')

      write(77) header(1:nheader)
      write(77) ((sngl(dx*(i-1)), sngl(dy*(j-1)), sngl(b(i+(j-1)*nx)*abs(vex)), i=1,nx), j=1,ny)
      write(77) footer(1:nfooter)

      write(77) part1(1:npart1)//'B'//part2(1:npart2)
      write(77) sngl(b(1:nn))

      close(77)

      open(unit=77, file=(foldername//'/VTK/SeaLevel'//cstep//'.vtk'), status='unknown', &
          form='unformatted', access='stream')

      write(77) header(1:nheader)
      write(77) ((sngl(dx*(i-1)), sngl(dy*(j-1)), sngl(sealevel*abs(vex)), i=1,nx), j=1,ny)
      write(77) footer(1:nfooter)

      write(77) part1(1:npart1)//'SL'//part2(1:npart2)
      write(77) ( sngl(sealevel), i=1, nn )

      close(77)

    end if

    return

end subroutine Fastscape_Named_VTK
