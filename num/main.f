c ======================================================================
      Program Th
c ======================================================================
      include 'th.fh'

      Integer  iLoop, nLOOPs
      Real*8   Quality

c ======================================================================
c number of adaptive loops
      nLOOPs = 5

      call read_cfg
      call create_mesh

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

         call solve

c === draw mesh and solution as ps-figures
         If(iLoop.EQ.1) Then
            Call isolines_demo(SOL_U, nv,vrt, nt,tri, nb,bnd,
     &           'sol_u_ini.ps', 50, '')
         Else
            Call isolines_demo(SOL_U, nv,vrt, nt,tri, nb,bnd,
     &           'sol_u_fin.ps', 50, '')
            Call draw_mesh('meshF.ps', 'final')
         End if
         If(iLoop.eq.nLOOPs) Goto 500

         call adapt_mesh(Quality)

  500    continue
      End do

c === testing the results
      If(Quality.LT.0.5) Stop 911
      Stop

      End

