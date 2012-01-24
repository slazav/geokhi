c ======================================================================
      Program Stokes
c ======================================================================
c The adaptive solution for the following bvp:
c
c  -div grad u  + grad p = 0   in Omega
c        div u           = 0   in Omega
c
c                      u = u_0 on dOmega_1
c                      u = 0   on dOmega_2
c                  du/dn = 0   on dOmega_3
c
c where Omega is a domain with a circle hole,  dOmega_1 is the side 
c at x=-2, dOmega_3 is the side at x=2, and dOmega_2 is the rest of 
c the boundary. The non-homogeneous boundary condition is 
c
c    u_0 = { 1-4y^2), 0 }.
c
c We use the P2 finite elements for the velocity u and P1 finite elements 
c for the pressure p, which are known to be stable. The discretization 
c method results in a symmetric indefinite matrix.
c
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,nfmax,namax
      parameter(nvmax = 50 000, ntmax = 2*nvmax, nbmax = 5 000)
      parameter(nfmax = 200 000, namax = 2 000 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 4 000 000, MaxWi = 6 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)
      Real     ANI_CrvFunction

c ======================================================================
c Mesh definition
c ======================================================================

c mesh generator data specifying domain via in the segment format
      double precision vbr(2,nbmax)
      integer          Nbr

      double precision vrt(2,nvmax)
      integer          tri(3,ntmax), bnd(2,nbmax)
c ... AFT2D library function
      Integer   aft2dfront
      EXTERNAL  aft2dfront

      EXTERNAL  mkcell

      Integer  nv, nvfix, labelV(nvmax), fixedV(1)
      Integer  nb, nbfix, labelB(nbmax), fixedB(1)
      integer  nt, ntfix, labelT(ntmax), fixedT(1)
      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      DATA     nvfix/0/, fixedV/0/,  nbfix/0/, ntfix/0/, nc/0/


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), SOL(nfmax), RES(nfmax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(1)

      Integer  Dbc 
      EXTERNAL Dbc, FEM2Dext


c ======================================================================
c for library aniLU
c ======================================================================
      Integer  symbolic(2), numeric(2), sys
      Real*8   lucontrol(20), luinfo(90)


c ======================================================================
c for library aniLMR
c ======================================================================
      Real*8   Lp
      Real*8   Metric(3,nvmax)


c ======================================================================
c for library aniMBA
c ======================================================================
      Integer  control(6), nEStar
      Real*8   Quality
      EXTERNAL CrvFunction_user


c LOCAL VARIABLEs
      Integer  i, j, dummy, iLoop, nLOOPs, nRow, nCol, iERR
      Integer  nr, ibc, iv1,iv2
      Real*8   rmax, h, x, y, eBC(2)
      Character*30 file_mesh, file_velocity

c ======================================================================
c number of adaptive loops
      nLOOPs = 5

C Read input file that contains coordinates of boundary points
      call mkcell(Nbr, vbr, labelB)

      Write(*,*) Nbr

C Generate a mesh  starting  from boundary mesh
      ierr=aft2dfront(
     &           0, dummy, Nbr, vbr,   ! segment data
     &           nv, vrt,              ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
      If (ierr.ne.0) stop ' error in function aft2dfront'
      Write(*,5000) Nbr, nt, nv

        labelB(3) = 2
        labelB(8) = 3

C      labelb(2)=2
C      labelb(6)=3


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine Ddiff
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)

      write (*,*) 'VRT:'
      Do i=1, nv
        write (*,*) vrt(1,i), vrt(2,i), labelV(i)
      End do

      write (*,*) 'BRD:'
      Do i=1, nb
        write (*,*) bnd(1,i), bnd(2,i), labelB(i)
      End do

      write (*,*) 'TRI:'
      Do i=1, nt
        write (*,*) tri(1,i), tri(2,i), tri(3,i), labelT(i)
      End do


c === general sparse matrix in a 0-based CSC format used in UMFPACK
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &        nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c  ===   call the driver for LU factorization and solution
         Call CSC2CSC0(nCol, IA, JA)

c set up default control parameters & print only error messages
         Call umf4def(lucontrol)
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol, nCol, IA,JA,A, symbolic,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5001

c numeric factorization
         Call umf4num(IA,JA,A, symbolic,numeric,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5002

c free the symbolic analysis data
         Call umf4fsym(symbolic)

c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, SOL, RHS, numeric, lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5003

c free the numeric factorization data
         Call umf4fnum (numeric)


c check the residual
         Call mulAcsc0(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') 'LU:  maximal norm of residual: ', rmax


c === draw mesh and solution isolines as ps-figures (mid-edges are ignored)
         nr  = (nRow - 3*nv) / 2

c  PostScript file names must have extension .ps
c  demo graphics has been activated
c  Isolines of velocity components can be done with the following calls:
c  Call isolines(SOL(iux), nv,vrt, nt,tri, nb,bnd,'velocity_x.ps',20)
c  Call isolines(SOL(iuy), nv,vrt, nt,tri, nb,bnd,'velocity_y.ps',20)
         If(iLoop.EQ.1) Then
            Call isolines_demo(SOL(1), nv,vrt, nt,tri, nb,bnd,
     &           'streamlines_ini.ps', 20, '')
            Call graph_demo(nv,vrt, nt,tri,
     &           'mesh_ini.ps', 'Initial quasi-uniform mesh')
            Call draw_matrix(nRow, IA, JA, 'matrix_ini.ps')
         Else
            file_mesh     = 'mesh_final.ps'
            file_velocity = 'streamlines_fin.ps'
            Call isolines_demo(SOL(1), nv,vrt, nt,tri, nb,bnd,
     &           'streamlines_fin.ps', 40, '')
            Call graph_demo(nv,vrt, nt,tri,
     &           'mesh_final.ps', 'final mesh')
            Call draw_matrix(nRow, IA, JA, "matrix_fin.ps")
         End if

         If(iLoop.eq.nLOOPs) Goto 500


c  ===   generate metric (from SOL) optimal for the L_p norm
c        Lp = 0             ! maximum norm
         Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL(1),
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === generate the adaptive mesh to u_x
         nEStar = 20000
         control(1) = nEStar/10  ! MaxSkipE
         control(2) = nEStar*10  ! MaxQItr
         control(3) = 16+1    ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
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

 500     Continue
      End do


c === testing the results
      If(Quality.LT.0.5) Stop 911
      Stop 


c error messages
 5000 Format(
     & 'Program generates a mesh using the advanced front technique.',/,
     & '  The initial front has', I4, ' edges.',/,
     & '  The final mesh has', I4, ' triangles and', I4, ' vertices.',/,
     & 'Created Postscript figure $(ANIHOME)/bin/mesh_final.ps.',/)
 5001 Continue
      Write(*,*) 'Error occurred in umf4sym: ', luinfo(1)
      Stop 911

 5002 Continue
      Write(*,*) 'Error occurred in umf4num: ', luinfo(1)
      Stop 911

 5003 Continue
      Write(*,*) 'Error occurred in umf4sol: ', luinfo(1)
      Stop 911

      End

