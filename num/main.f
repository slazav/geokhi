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
         call calc_gu2
         If(iLoop.ne.nLOOPs) call adapt_mesh(Quality, SOL_GU2)
      End do

      call draw_gu2('res/sol_gu2.ps')
      Call draw_u('res/sol_u.ps')
      Call draw_mesh('res/mesh_u.ps', 'final')


      call write_line(0D0,0D0, 0D0,0.005D0,
     &   500,SOL_U, 'res/sol_u_1.dat')
      call write_line(Dd*0.25,Dd*0.25, 0D0,0.005D0,
     &   500,SOL_U, 'res/sol_u_2.dat')
      call write_line(Dd*0.40,Dd*0.40, 0D0,0.005D0,
     &   500,SOL_U, 'res/sol_u_3.dat')
      call write_line(Dd*0.48,Dd*0.48, 0D0,0.005D0,
     &   500,SOL_U, 'res/sol_u_4.dat')
      call write_line(Dd*0.50,Dd*0.50, 0D0,0.005D0,
     &   500,SOL_U, 'res/sol_u_5.dat')

      call write_line(0D0,0D0, 0D0,0.003D0,
     &   500,SOL_GU2, 'res/sol_gu_1.dat')
      call write_line(Dd*0.25,Dd*0.25, 0D0,0.003D0,
     &   500,SOL_GU2, 'res/sol_gu_2.dat')
      call write_line(Dd*0.40,Dd*0.40, 0D0,0.003D0,
     &   500,SOL_GU2, 'res/sol_gu_3.dat')
      call write_line(Dd*0.48,Dd*0.48, 0D0,0.003D0,
     &   500,SOL_GU2, 'res/sol_gu_4.dat')
      call write_line(Dd*0.50,Dd*0.50, 0D0,0.003D0,
     &   500,SOL_GU2, 'res/sol_gu_5.dat')

c === testing the results
      If(Quality.LT.0.5) write (*,*) '!!! low mesh quality ', Quality


      Write(*,'(/,A,I2)') 'SOLVE T'
      call solve_t
      Call draw_t('res/sol_t.ps')

      call write_line(0D0,0D0, 0D0,0.005D0,
     &   500,SOL_T, 'res/sol_t_1.dat')
      call write_line(Dd*0.25,Dd*0.25, 0D0,0.005D0,
     &   500,SOL_T, 'res/sol_t_2.dat')
      call write_line(Dd*0.40,Dd*0.40, 0D0,0.005D0,
     &   500,SOL_T, 'res/sol_t_3.dat')
      call write_line(Dd*0.48,Dd*0.48, 0D0,0.005D0,
     &   500,SOL_T, 'res/sol_t_4.dat')
      call write_line(Dd*0.50,Dd*0.50, 0D0,0.005D0,
     &   500,SOL_T, 'res/sol_t_5.dat')


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


