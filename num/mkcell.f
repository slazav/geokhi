      subroutine mkcell(nbr, vbr, labelB)

        double precision vbr(2,*)
        integer nbr, labelB(*)

C       read cell dimensions
        real*8 L, D, Dw, Lw, Dd, Ld
        DATA L/1.0/, D/1.0/, Lw/0.3/, Dw/0.5/, Dd/0.02/, Ld/0.01/
        real*8 CFG_VAL
        character*64 CFG_KEY
        integer i

        open(54,FILE='cell.cfg')
 1001   read(54,*,END=1002) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'L') then
          L=CFG_VAL
        elseif (CFG_KEY.EQ.'D') then
          D=CFG_VAL
        elseif (CFG_KEY.EQ.'Lw') then
          Lw=CFG_VAL
        elseif (CFG_KEY.EQ.'Dw') then
          Dw=CFG_VAL
        elseif (CFG_KEY.EQ.'Ld') then
          Ld=CFG_VAL
        elseif (CFG_KEY.EQ.'Dd') then
          Dd=CFG_VAL

        else
          write(*,'(A,A20)')
     *     'warning: unknown parameter in cfg-file: ', CFG_KEY
        endif
        goto 1001
 1002   close(54)

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

        do i=1,8
          labelB(i) = 1
        enddo

        labelB(3) = 2
        labelB(8) = 3

        nbr = 9
        return
      end