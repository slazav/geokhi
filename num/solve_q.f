C ======================================================================
      Subroutine draw_Q(filename)
        include 'th.fh'
        character*(*) filename
        integer ntri(3,nt),nnt,n,i

        nnt=0
        do n=1,nt
          if (labelT(n).eq.1) then
            nnt=nnt+1
            do i=1,3
              ntri(i,nnt)=tri(i,n)
            enddo
          endif
        enddo

        Call isolines(SOL_Q, nv,vrt, nnt,ntri, nb,bnd, filename, 50)
      End

C ======================================================================
      Subroutine solve_Q
        include 'th.fh'

c Local variables
        Integer n,i,j
        Real*8 GrP0(nt), Gr(2)
        integer iErr

        Real*8 xy1(2),xy2(2),xy3(2), u1,u2,u3, S
        Real*8 dd(nv), dist

C ======================================================================

        Write(*,*) 'finding Q'

        Do n = 1, nt
          do i=1,2
            xy1(i)=vrt(i,tri(1,n))
            xy2(i)=vrt(i,tri(2,n))
            xy3(i)=vrt(i,tri(3,n))
          enddo
          u1=SOL_U(tri(1,n))
          u2=SOL_U(tri(2,n))
          u3=SOL_U(tri(3,n))

          if (labelT(n).eq.1) then 
            S = PHYS_SIGMA
          else
            S = PHYS_SIGMAd
          endif

          call GradUtri(xy1,xy2,xy3, u1,u2,u3, Gr)
          GrP0(n)=(Gr(1)**2+Gr(2)**2)*S
        end do

c      call P02P1(nv,nt,vrt,tri, GrP0, SOL_Q, MaxWi, iW, iErr)
c       recalc P0 -> P1 (don't use labelT==2 values)
        Do n = 1, nv
          dd(n)=0D0
          SOL_Q(n)=0D0
        enddo

        Do n = 1, nt
          if (labelT(n).eq.1) then
c           find tri center
            xy1(1)=0
            xy1(2)=0
            do i=1,3
              xy1(1)=xy1(1) + vrt(1,tri(i,n))
              xy1(2)=xy1(2) + vrt(2,tri(i,n))
            enddo
            xy1(1)=xy1(1)/3D0
            xy1(2)=xy1(2)/3D0

            do i=1,3
              j = tri(i,n)
              dist = dsqrt((vrt(1,j)-xy1(1))**2 +
     &                     (vrt(2,j)-xy1(2))**2)
              SOL_Q(j) = SOL_Q(j) + GrP0(n)/dist
              dd(j) = dd(j)+1/dist
            enddo
          endif
        enddo

        Do n = 1, nv
          if (dd(n).gt.0) then
            SOL_Q(n) = SOL_Q(n)/dd(n)
          else
            SOL_Q(n) = 0
          endif
        enddo

        call draw_q('ps/sol_q.ps')

        Return
      End


