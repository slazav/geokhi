c ======================================================================
      Program Th
c ======================================================================
      include 'th.fh'

      Integer  iLoop, nLOOPs, i
      Real*8   Quality

c ======================================================================
c number of adaptive loops
      nLOOPs = 6

      call read_cfg
      call create_mesh

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> U LOOP: ', iLoop
         call solve_u
         call solve_q
         If(iLoop.ne.nLOOPs) call adapt_mesh(Quality, SOL_Q)
      End do

      call draw_q('res/sol_q.ps')
      Call draw_u('res/sol_u.ps')
      Call draw_mesh('res/mesh_u.ps', 'final')


      call write_lines(SOL_U, 'u')
      call write_lines(SOL_Q, 'q')

c === testing the results
      If(Quality.LT.0.5) write (*,*) '!!! low mesh quality ', Quality


      Write(*,'(/,A,I2)') 'SOLVE T'
      call solve_t
      Call draw_t('res/sol_t.ps')

      call write_lines(SOL_T, 't')


      Stop
      End


      subroutine write_line(x1,x2,y1,y2,np,u, filename)
        include 'th.fh'
        character*(*) filename
        real*8 x1,x2,y1,y2,dx,dy
        real*8 xy(2,np), out(np)
        real*8 u(*)
        integer np,info(3)
        integer i

        do i=1,np
          xy(1,i)=x1 + (x2-x1)*(i-1)/(np-1)
          xy(2,i)=y1 + (y2-y1)*(i-1)/(np-1)
        end do
        info(1)=1
        info(2)=1
        call LINTRP(nt,tri,nv,vrt,1,u, np,xy, out, iW,rW, info)
        open (57, FILE=filename)
        do i=1,np
          write (57,*) xy(1,i), xy(2,i), out(i)
        enddo
        close (57)
      end


      subroutine write_lines(u,par)
        include 'th.fh'
        real*8 x(5),y1,y2
        double precision u(*)
        character*(64) filename
        character*(*) par

        integer i, np

        x(1)=0D0
        x(2)=DIM_Dd*0.25
        x(3)=DIM_Dd*0.40
        x(4)=DIM_Dd*0.48
        x(5)=DIM_Dd*0.50
        y1=0D0
        y2=0.005D0
        np=500

        do i=1,5
          write(filename,'(A,A,I1,A)') 'res/sol_', par, i,'.dat'
          call write_line(x(i),x(i),y1,y2, np, u, filename)
        enddo
      end
