      Subroutine mmul(A1,B1,C1,D1, A2,B2,C2,D2)
        real*8 A1,B1,C1,D1, A2,B2,C2,D2, A,B,C,D
        A = A1 * A2 + B1 * C2
        B = A1 * B2 + B1 * D2
        C = C1 * A2 + D1 * C2
        D = C1 * B2 + D1 * D2
        A2 = A
        B2 = B
        C2 = C
        D2 = D
      end

C ======================================================================
      Subroutine solve_ABCD(mx, nx,ny)
        include 'th.fh'

        integer i,j, nx,ny
        real*8 mx,y, N0, A0, E
        real*8 N, A2,  dy
        real*8 AA,BB,CC,DD

        open (57, FILE='dat/fit_n.dat')
        write (57,*) '## z N0 A0 ERR'

        AA=1D0
        BB=0D0
        CC=0D0
        DD=1D0

        write(*,*) '>>>',AA,BB,CC,DD
          call mmul(1D0,0D0,0D0,1D0,  AA, BB, CC, DD)
        write(*,*) '>>>',AA,BB,CC,DD

        do i= -ny, ny
          y=dabs((DIM_L - DIM_Lw) *i)/ny/2.0
          call fit_line(mx, y, nx, SOL_N,  N0, A0, E)

c         we need N(1 + A2 r^2/2)
c         we have PHYS_N0 + N0 - A0 r^2
          N = PHYS_N0 + N0
          A2 = - 2.0D0 * A0/N

          write (57,*) y*i/abs(i), N, A2, E

          dy = (DIM_L - DIM_Lw)/ny/2.0

          call mmul(1D0 + A2 * dy**2, dy/N, N * A2 * dy, 1D0,
     c              AA, BB, CC, DD)
        enddo
        close(57)



        open (57, FILE='dat/fit_abcd.dat')
        write(57,*) AA,BB,CC,DD
        call mmul(1, DIM_Ln, 0, 1,  AA, BB, CC, DD)
        write(57,*) AA,BB,CC,DD
        close(57)
        Return
      End

c      fit with  N = N0 - A0 r^2
      subroutine fit_line(mx,y,np,u,  N0, A0, maxer)
        include 'th.fh'
        real*8 mx, y, u(*)
        real*8 xy(2,np), out(np)
        integer np,info(3)
        integer i,fd

        real*8 s0,s1,s2,sf,sfr, N0, A0, er, maxer

        do i=1,np
          xy(1,i)=mx*(i-1)/(np-1)
          xy(2,i)=y
        end do
        info(1)=1
        info(2)=1
        call LINTRP(nt,tri,nv,vrt,1,u, np,xy, out, iW,rW, info)

        s0=0
        s1=0
        s2=0
        sf=0
        sfr=0
        do i=1,np
          s0 = s0 + 1D0
          s1 = s1 + xy(1,i)**2
          s2 = s2 + xy(1,i)**4
          sf = sf + out(i)
          sfr = sfr + out(i)*xy(1,i)**2
        enddo
        A0 = (s0*sfr-s1*sf)/(s1**2-s0*s2)
        N0 = (A0*s1 + sf)/s0
c       check max error
        maxer = 0
        do i=1,np
          er = dabs((N0 - A0*xy(1,i)**2 - out(i)) / out(i))
          if (er.gt.maxer) maxer=er
        enddo
      end
