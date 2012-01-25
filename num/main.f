c ======================================================================
      Program Th
c ======================================================================
      include 'th.fh'

c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), RES(nfmax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(1)

      Integer  Dbc 
      EXTERNAL Dbc, FEM2Dext
      Real*8   adapt_mesh
      EXTERNAL adapt_mesh


c ======================================================================
c for library aniLU
c ======================================================================
      Integer  symbolic(2), numeric(2), sys
      Real*8   lucontrol(20), luinfo(90)


c LOCAL VARIABLEs
      Integer  i, j, dummy, iLoop, nLOOPs, nRow, nCol
      Real*8   rmax, h, x, y, eBC(2)

      Real*8   Quality


        Real*8   Lp
        Real*8   Metric(3,nvmax)

        Integer  control(6), nEStar, iERR


c ======================================================================
c number of adaptive loops
      nLOOPs = 5

C Read input file that contains coordinates of boundary points
      call read_cfg
      call create_mesh

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine Ddiff
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)

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
         Call umf4sol(sys, SOL_U, RHS, numeric, lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5003

c free the numeric factorization data
         Call umf4fnum (numeric)


c check the residual
         Call mulAcsc0(nRow, IA, JA, A, SOL_U, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') 'LU:  maximal norm of residual: ', rmax


c === draw mesh and solution as ps-figures
         If(iLoop.EQ.1) Then
            Call isolines_demo(SOL_U(1), nv,vrt, nt,tri, nb,bnd,
     &           'sol_u_ini.ps', 50, '')
         Else
            Call isolines_demo(SOL_U(1), nv,vrt, nt,tri, nb,bnd,
     &           'sol_u_fin.ps', 50, '')
            Call graph_demo(nv,vrt, nt,tri, 'meshF.ps', '')
         End if
         If(iLoop.eq.nLOOPs) Goto 500

c         Quality = adapt_mesh()
c         Write (*,*) 'Quality: ', Quality

C =====================================================================

        Write(*,*) 'create metric from solution'
c  ===  generate metric (from SOL) optimal for the L_p norm
c       Lp = 0             ! maximum norm
        Lp = 1             ! L_1 norm
        Call Nodal2MetricVAR(SOL_U(1),
     &          vrt, nv, tri, nt, bnd, nb, Metric,
     &          MaxWr, rW, MaxWi, iW)

        If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


        Write(*,*) 'generate the adaptive mesh to u'
c === generate the adaptive mesh to u
         nEStar = 10000
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
C =====================================================================



  500    continue
      End do

c === testing the results
      If(Quality.LT.0.5) Stop 911
      Stop



c error messages
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

