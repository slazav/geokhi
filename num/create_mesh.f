C ====== make initial mesh from cell dimensions
      subroutine create_mesh()
        include 'th.fh'
        Integer   aft2dfront
        EXTERNAL  aft2dfront
        double precision vbr(2,nbmax)
        integer nbr

        integer  ipIRE,ipWork, MaxWiWork
        DATA nvfix/0/, fixedV/0/,  nbfix/0/, ntfix/0/, nc/0/

        Write (*,*) 'Creating mesh.'
        Write (*,*) 'Cell dimensions:'
        Write (*,*) '  L = ', L
        Write (*,*) '  D = ', D
        Write (*,*) '  Lw = ', Lw
        Write (*,*) '  Dw = ', Dw
        Write (*,*) '  Ld = ', Ld
        Write (*,*) '  Dd = ', Dd

        vbr(1,1) = 0
        vbr(2,1) = 0

        vbr(1,2) = 0
        vbr(2,2) = (L-Lw)/2.0

        vbr(1,3) = Dw/2.0
        vbr(2,3) = vbr(2,2)

        vbr(1,4) = vbr(1,3)
        vbr(2,4) = L/2.0

        vbr(1,5) = D/2.0
        vbr(2,5) = L/2.0

        vbr(1,6) = D/2.0
        vbr(2,6) = Ld/2.0

        vbr(1,7) = Dd/2.0
        vbr(2,7) = Ld/2.0

        vbr(1,8) = Dd/2.0
        vbr(2,8) = 0

        vbr(1,9) = 0
        vbr(2,9) = 0

        nbr = 9


C Generate a mesh
        ierr=aft2dfront(
     &           0, dummy, nbr, vbr,   ! segment data
     &           nv, vrt,              ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
        If (ierr.ne.0) stop ' error in function aft2dfront'

        Write(*,5101) Nbr, nt, nv

        Write(*,*) '   Writing mesh into mesh0.ps'
        Call graph_demo(nv,vrt,nt,tri, 'mesh0.ps', '')

        labelB(3) = 2
        labelB(8) = 3

C Refine mesh
         ipIRE = 1
         ipWork = ipIRE + 3 * nt
         MaxWiWork = MaxWi - 3 * nt


         call refine_mesh

c         Call uniformRefinement(
c     &        nv, nvmax, nb, nbmax, nt, ntmax,
c     &        vrt, bnd, labelB, tri, labelT,
c     &        ANI_CrvFunction, crv, labelC, iW(ipIRE),
c     &        rW, 1, iW(ipWork), MaxWiWork)

c        Write(*,5102) nt, nv
c        Write(*,*) '   Writing refined mesh into mesh1.ps'
c        Call graph_demo(nv,vrt,nt,tri, 'mesh1.ps', '')

c        Call Delaunay(
c     &        nv, nt, vrt,  tri, MaxWi, iW)
c        Call smoothingMesh(
c     &        nv, nt, vrt,  tri, MaxWi, iW)


        Write(*,5102) nt, nv
        Write(*,*) '   Writing smoothed mesh into mesh2.ps'
        Call graph_demo(nv,vrt,nt,tri, 'mesh2.ps', '')

        return

c =======
 5101   format(
     &  '  The initial front has', I4, ' edges.',/,
     &  '  The initial mesh has', I4, ' triangles and',
     &  I4, ' vertices.')
 5102   format(
     &  '  The refined mesh has', I4, ' triangles and',
     &  I4, ' vertices.')


      end


      subroutine refine_mesh()
        include 'th.fh'
        Integer  nEStar, iERR
        Integer  control(6)
        Real*8   Quality
        Integer  MetricFunction
        External MetricFunction

c === generate adaptive mesh
        control(1) = 500     !  MaxSkipE
        control(2) = 50000   !  MaxQItr
        control(3) = 1       !  status
        control(4) = 1       !  flagAuto
        control(5) = 1       !  iPrint:   average level of output information
        control(6) = 0       !  iErrMesgt: only critical termination allowed

        Quality = 0.8D0      !  request shape-regular triangles in metric
        nEStar  = 200       !  desired number of triangles

      Call mbaAnalytic(
c group (M)
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, CrvFunction,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &      nEStar, Quality, control, MetricFunction,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)

      end
 

C =====================================================================
      Integer Function MetricFunction(x, y, Metric)
C =====================================================================
C  This routine creates a metric at the given point (x,y). The
C  metric is a 2x2 positive definite symmetric tensor:
C                M11   M12
C      Metric =
C                M12   M22
C
C  Only the upper triangular part of 2x2 array Metric may be defined.
C
C  In this example, the metric is constant and isotropic.
C =====================================================================
      Real*8  x, y, Metric(2, 2)
c      Metric(1,1) = (D-x)*(L/2.0-y)
      Metric(1,1) = (x*y)**2
      Metric(1,2) = 0
      Metric(2,1) = 0
      Metric(2,2) = Metric(1,1)
      MetricFunction = 0
      Return
      End
