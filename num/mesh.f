
C =====================================================================
      subroutine draw_mesh(filename, desc)
        include 'th.fh'
        character*(*) filename, desc
        Write(*,*) '   Writing ', desc,' mesh into ', filename
        Write(*,*) nt, ' triangles and ', nv, 'vertices.'
        Write(*,*) nc, ' curves ', nb, 'boundary edges'
        call graph_demo(nv,vrt, nt,tri, filename, '')
c        call graph(nv,vrt, nt,tri, filename)
        return
      end

C ====== make initial mesh from cell dimensions
      subroutine create_mesh()
        include 'th.fh'
        Integer   aft2dfront
        EXTERNAL  aft2dfront
        double precision vbr(2,nbmax), tmp(2)
        integer nbr, dummy

        integer  ipIRE,ipWork, MaxWiWork, iERR
        DATA nvfix/0/, fixedV/0/,  nbfix/0/, ntfix/0/, nc/0/

        Write (*,*) 'Creating mesh.'

        nbr=1

        vbr(1,nbr) = 0D0
        vbr(2,nbr) = 0D0

        nbr=nbr+1
        vbr(1,nbr) = 0D0
        vbr(2,nbr) = (DIM_L-DIM_Lw)/2.0D0

        nbr=nbr+1
        vbr(1,nbr) = DIM_Dw/2.0D0
        vbr(2,nbr) = vbr(2,2)

        nbr=nbr+1
        vbr(1,nbr) = vbr(1,3)
        vbr(2,nbr) = DIM_L/2.0D0

        nbr=nbr+1
        vbr(1,nbr) = DIM_D/2.0D0
        vbr(2,nbr) = DIM_L/2.0D0

        nbr=nbr+1
        vbr(1,nbr) = DIM_D/2.0D0
        vbr(2,nbr) = DIM_Ld/2.0D0


        nbr=nbr+1 
        call crv_func(0,tmp,1)
        vbr(1,nbr) = tmp(1)
        vbr(2,nbr) = tmp(2)

        nbr=nbr+1 
        call crv_func(dasin(1D0)/2D0,tmp,1)
        vbr(1,nbr) = tmp(1)
        vbr(2,nbr) = tmp(2)

        nbr=nbr+1 
        call crv_func(dasin(1D0),tmp,1)
        vbr(1,nbr) = tmp(1)
        vbr(2,nbr) = tmp(2)


        nbr=nbr+1
        vbr(1,nbr) = DIM_Dd/2.0D0
        vbr(2,nbr) = 0D0

        nbr=nbr+1
        vbr(1,nbr) = vbr(1,1)
        vbr(2,nbr) = vbr(2,1)

C Generate a mesh
        Write(*,5101) Nbr
        iERR=aft2dfront(
     &           0, dummy, nbr, vbr,   ! segment data
     &           nv, vrt,              ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
        If (iERR.ne.0) stop ' error in function aft2dfront'

        Call draw_mesh('res/mesh0.ps', 'initial')

        labelB(1) = 5     ! center line
        labelB(2) = 4     ! cell walls
        labelB(4) = 4
        labelB(5) = 4
        labelB(3) = 3     ! u=u0/2
        labelB(7) = 10
        labelB(8) = 10
        labelB(10) = 2 ! center plane (u=0)

        nc = 2
        crv(1,7) = 0D0
        crv(2,7) = dasin(1D0)/2D0
        labelC(7) = 1
        crv(1,8) = dasin(1D0)/2D0
        crv(2,8) = dasin(1D0)
        labelC(8) = 1

c       add

        nv=nv+1
        vrt(1,nv) = DIM_L/2.0D0
        vrt(2,nv) = 0
        labelV(nv) = 0

        nb=nb+1
        bnd(1,nb) = 6
        bnd(2,nb) = nv
        labelB(nb) = 11

        nb=nb+1
        bnd(1,nb) = nv
        bnd(2,nb) = 10
        labelB(nb) = 11


        nt=nt+1
        tri(1,nt) = 10
        tri(2,nt) = 9
        tri(3,nt) = 8
        labelT(nt) = 2

        nt=nt+1
        tri(1,nt) = 10
        tri(2,nt) = 8
        tri(3,nt) = 7
        labelT(nt) = 2

        nt=nt+1
        tri(1,nt) = 10
        tri(2,nt) = 7
        tri(3,nt) = 6
        labelT(nt) = 2

        nt=nt+1
        tri(1,nt) = 10
        tri(2,nt) = 6
        tri(3,nt) = nv
        labelT(nt) = 2

        Call draw_mesh('res/mesh0a.ps', 'initial')
      Call saveMani(nv, nvfix, vrt, labelV, fixedV,
     &              nb, nbfix, bnd, labelB, fixedB,
     &              nc,        crv, labelC,
     &              nt, ntfix, tri, labelT, fixedT,
     &     "res/mesh0a.ani")

c Refine mesh
        call refine_mesh
        call draw_mesh('res/mesh1.ps', 'refined')

        return

c =======
 5101   format(
     &  '  The initial front has ', I4, ' edges.')

      end

C =====================================================================
      subroutine crv_func(tc, xyc, iFunc)
        include 'th.fh'
        Real*8 tc, xyc(2)
        integer iFunc
c       tc = 0..pi/2
        xyc(1) = DIM_Dd/2D0 + DIM_Dr * (1D0 - dsin(tc))
        xyc(2) = DIM_Ld/2D0 - DIM_Dr * (1D0 - dcos(tc))
      end

C =====================================================================
      subroutine refine_mesh()
        include 'th.fh'
        Integer  nEStar, iERR
        Integer  control(6)
        Real*8   Quality
        Integer  metric_func
        External metric_func
        External crv_func

c === generate adaptive mesh
        nEStar  = 20        !  desired number of triangles (not used?)
        control(1) = nEStar/10   !  MaxSkipE
        control(2) = nEStar*10   !  MaxQItr
        control(3) = 1       !  status
        control(4) = 1       !  flagAuto
        control(5) = 0       !  iPrint:   average level of output information
        control(6) = 0       !  iErrMesgt: only critical termination allowed

        Quality = 0.60D0     !  request shape-regular triangles in metric

c      Call mbaAnalytic(
      call mbaFixShape(
c group (M)
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, crv_func,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &      nEStar, Quality, control, metric_func,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)
      write (*,*) 'quality after refining: ', Quality

c        Call smoothingMesh(
c     &        nv, nt, vrt,  tri, MaxWi, iW)
c        Call Delaunay(
c     &        nv, nt, vrt,  tri, MaxWi, iW)

      end
C =====================================================================
      Integer Function metric_func(x, y, M)
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
      metric_func = 0
      Return
      End


C =====================================================================
      subroutine adapt_mesh(Quality, F)
        include 'th.fh'

        Real*8   F(nfmax)

        Real*8   Lp
        Real*8 Metric(3,nvmax)

        Integer  control(6), nEStar, iERR
        Real*8   Quality
        External crv_func

c  ===  generate metric (from SOL) optimal for the L_p norm
c       Lp = 0             ! maximum norm
        Lp = 1             ! L_1 norm
        Call Nodal2MetricVAR(F,
     &          vrt, nv, tri, nt, bnd, nb, Metric,
     &          MaxWr, rW, MaxWi, iW)

        If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


        Write(*,*) 'generate the adaptive mesh...'
c === generate the adaptive mesh to u
         nEStar = TRI_NUM
         control(1) = nEStar/10  ! MaxSkipE
         control(2) = nEStar*10  ! MaxQItr
         control(3) = 32+1    ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 0       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         Quality = 0.6

         Call mbaNodal(
c group (M)
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, crv_func,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &        nEStar, Quality, control, Metric,
c group (W)
     &        MaxWr, MaxWi, rW, iW, iERR)

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error if mbaNodal')
      write (*,*) 'mesh quality: ', Quality

        return
      end
