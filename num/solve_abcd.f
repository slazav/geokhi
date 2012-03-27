      Subroutine mmul(A, B)
        real*8 A(2,2), B(2,2), C(2,2)
        do i=1,2
          do j=1,2
            C(i,j) = A(1,j) * B(i,1) + A(2,j) * B(i,2)
          enddo
        enddo
        do i=1,2
          do j=1,2
            B(i,j) = C(i,j)
          enddo
        enddo
      end
C ======================================================================
      subroutine m_space(A, L, N0)
        real*8 A(2,2), L, N0
        A(1,1)=1D0
        A(2,1)=L/N0
        A(1,2)=0D0
        A(2,2)=1D0
      end
C ======================================================================
      subroutine m_lens(A, F)
        real*8 A(2,2), F
        A(1,1)=1D0
        A(2,1)=0D0
        A(1,2)=-1D0/F
        A(2,2)=1D0
      end
C ======================================================================
      subroutine m_quad(A, L, N0, A0) ! N = N0(1 + A0 R^2 / 2)
        real*8 A(2,2), L, N0, A0
        A(1,1)=1D0 + A0 * L**2
        A(2,1)=L/N0
        A(1,2)=N0 * A0 * L
        A(2,2)=1D0
      end
C ======================================================================
      function do_abcd(q, A)
        complex*16 q, do_abcd
        real*8 A(2,2)
        do_abcd = (A(1,1) * q + A(2,1))/(A(1,2) * q + A(2,2))
      end
C ======================================================================



C ======================================================================
      Subroutine solve_ABCD(mx, nx,ny)
        include 'th.fh'

        integer i,j, nx,ny
        real*8 mx,y, N0, A0, E
        real*8 N, A2,  dy
        real*8 A(2,2), DA(2,2)
        complex*16 Qi, Qo
        complex*16 do_abcd

        open (57, FILE='dat/fit_n.dat')
        write (57,*) '## z N0 A0 ERR'

        call m_space(A, DIM_Li, 1D0)

        do i= -ny, ny
          dy = (DIM_L - 2*DIM_Lw)/ny/2.0
          y=dabs((DIM_L - 2*DIM_Lw) *i)/ny/2.0
          call fit_line(mx, y, nx, SOL_N,  N0, A0, E)

c         we need N(1 + A2 r^2/2)
c         we have PHYS_N0 + N0 - A0 r^2
          N = PHYS_N0 + N0
          A2 = - 2.0D0 * A0/N

          write (57,*) y*i/abs(i), N, A2, E

          call m_quad(DA, dy, N, A2)
          call mmul(DA, A)
        enddo
        close(57)

        call m_space(DA, DIM_Lo, 1D0)
        call mmul(DA, A)

C        write(*,*) '>>>',AA,BB,CC,DD
        Qi=1D0/CMPLX(0D0, RAY_L/PI/RAY_W0**2)
        Qo=do_abcd(Qi, A)

        open (57, FILE='dat/fit_abcd.dat')
        write(57,*) A(1,1),A(2,1),A(1,2),A(2,2)

        write(*,*) 'R= ', 1D0/REAL(1D0/Qo), ' m'
        write(*,*) 'W= ', DSQRT(1D0/AIMAG(1D0/Qo) * RAY_L / PI)
     &    * 1000, ' mm'


        write(57,*) A(1,1),A(2,1),A(1,2),A(2,2)
        close(57)

        call m_space(A, DIM_Li, 1D0)
        call m_space(DA, DIM_L - 2D0 * DIM_Lw, PHYS_N0)
        call mmul(DA, A)
        call m_space(DA, DIM_Lo, 1D0)
        call mmul(DA, A)

c        write(*,*) '>>>',AA,BB,CC,DD
c        write(*,*) 'Q>>>', Qi, ' ', Qo

        Qi=1D0/CMPLX(0D0, RAY_L/PI/RAY_W0**2)
        Qo=do_abcd(Qi, A)

c        write(*,*) 'Q>>>', Qi, ' ', Qo

        write(*,*) 'R0= ', 1D0/REAL(1D0/Qo), 'm'
        write(*,*) 'W0= ', DSQRT(1D0/AIMAG(1D0/Qo) * RAY_L / PI)
     &    * 1000D0, ' mm'

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
