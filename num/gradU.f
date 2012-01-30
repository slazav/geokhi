C ======================================================================
      Subroutine draw_GU2(filename)
        include 'th.fh'
        character*(*) filename
        Call isolines_demo(SOL_GU2, nv,vrt, nt,tri, nb,bnd,
     &         filename, 50, '')
      End

C ======================================================================
      Subroutine calc_GU2
        include 'th.fh'

c Local variables
        Integer n,i
        Real*8 GrP0(nt), Gr(2)
        integer iErr

        Real*8 xy1(2),xy2(2),xy3(2), u1,u2,u3

C ======================================================================

        Do n = 1, nt
          do i=1,2
            xy1(i)=vrt(i,tri(1,n))
            xy2(i)=vrt(i,tri(2,n))
            xy3(i)=vrt(i,tri(3,n))
          enddo
          u1=SOL_U(tri(1,n))
          u2=SOL_U(tri(2,n))
          u3=SOL_U(tri(3,n))

          call GradUtri(xy1,xy2,xy3, u1,u2,u3, Gr)
          GrP0(n)=Gr(1)**2+Gr(2)**2
        end do

        call P02P1(nv,nt,vrt,tri, GrP0, SOL_GU2, MaxWi, iW, iErr)

        Return
      End
