
C =====================================================================
      subroutine draw_mesh(filename, desc)
        include 'th.fh'
        character*(*) filename, desc
        Write(*,*) '   Writing ', desc,' mesh into ', filename
        Write(*,*) nt, ' triangles and ', nv, 'vertices.'
        call graph_demo(nv,vrt, nt,tri, filename, '')
        return
      end

C ====== make initial mesh from cell dimensions
      subroutine create_mesh()
        include 'th.fh'
        Integer   aft2dfront
        EXTERNAL  aft2dfront
        double precision vbr(2,nbmax)
        integer nbr, dummy

        integer  ipIRE,ipWork, MaxWiWork, iERR
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
        Write(*,5101) Nbr
        iERR=aft2dfront(
     &           0, dummy, nbr, vbr,   ! segment data
     &           nv, vrt,              ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
        If (iERR.ne.0) stop ' error in function aft2dfront'

        Call draw_mesh('mesh0.ps', 'initial')

        labelB(3) = 2
        labelB(8) = 3

C Refine mesh

        call refine_mesh
        call draw_mesh('mesh1.ps', 'refined')

        return

c =======
 5101   format(
     &  '  The initial front has ', I4, ' edges.')

      end


C =====================================================================
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
        control(5) = 0       !  iPrint:   average level of output information
        control(6) = 0       !  iErrMesgt: only critical termination allowed

        Quality = 0.99D0      !  request shape-regular triangles in metric
        nEStar  = 3000       !  desired number of triangles

      Call mbaAnalytic(
c group (M)
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, ANI_CrvFunction,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &      nEStar, Quality, control, MetricFunction,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)

c        Call smoothingMesh(
c     &        nv, nt, vrt,  tri, MaxWi, iW)
c        Call Delaunay(
c     &        nv, nt, vrt,  tri, MaxWi, iW)

      end
C =====================================================================
      Integer Function MetricFunction(x, y, M)
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
      include 'th.fh'
      Real*8  x, y, M(2, 2)

c      Metric(1,1) = (D/2.0-x)*(L/2.0-y)/D/L*4
      M(1,1) = 1
      M(1,2) = 0
      M(2,1) = 0
      M(2,2) = M(1,1)
      MetricFunction = 0
      Return
      End


C =====================================================================
      subroutine adapt_mesh(Quality)
        include 'th.fh'

        Real*8   Lp
        Real*8 Metric(3,nvmax)

        Integer  control(6), nEStar, iERR
        Real*8   Quality

        Write(*,*) 'create metric from solution'
c  ===  generate metric (from SOL) optimal for the L_p norm
c       Lp = 0             ! maximum norm
        Lp = 1             ! L_1 norm
        Call Nodal2MetricVAR(SOL_U,
     &          vrt, nv, tri, nt, bnd, nb, Metric,
     &          MaxWr, rW, MaxWi, iW)

        If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


        Write(*,*) 'generate the adaptive mesh to u'
c === generate the adaptive mesh to u
         nEStar = 8000
         control(1) = nEStar/10  ! MaxSkipE
         control(2) = nEStar*10  ! MaxQItr
         control(3) = 16+1    ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 0       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         Quality = 0.6

         Call mbaNodal(
c group (M)
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, ANI_CrvFunction,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &        nEStar, Quality, control, Metric,
c group (W)
     &        MaxWr, MaxWi, rW, iW, iERR)

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error if mbaNodal')

        return
      end
