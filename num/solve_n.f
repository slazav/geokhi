C ======================================================================
      Subroutine draw_n(filename)
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

        Call isolines(SOL_N, nv,vrt, nnt,ntri, nb,bnd, filename, 50)
      End

C ======================================================================
      Subroutine solve_N
        include 'th.fh'

c Local variables
        Integer n

C ======================================================================

        Do n = 1, nv
            SOL_N(n) = PHYS_N0 + SOL_T(n)*PHYS_DNDT
        enddo

        Return
      End


