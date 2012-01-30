c ====== Read configuration file
      subroutine read_cfg()
        include 'th.fh'
c === set default values
        DATA L/1.0/, D/1.0/, Lw/0.3/, Dw/0.5/,
     &       Ld/0.01/, Dd/0.02/,Dr/0.0/
        DATA U0/1.0/,
     &       INI_TRI_NUM/200/, FIN_TRI_NUM/8000/

        real*8 CFG_VAL
        character*64 CFG_KEY

        open(54,FILE='cell.cfg')
 1101   read(54,*,END=1102) CFG_KEY, CFG_VAL
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
        elseif (CFG_KEY.EQ.'Dr') then
          Dr=CFG_VAL

        elseif (CFG_KEY.EQ.'U0') then
          U0=CFG_VAL

        elseif (CFG_KEY.EQ.'INI_TRI_NUM') then
          INI_TRI_NUM=CFG_VAL
        elseif (CFG_KEY.EQ.'FIN_TRI_NUM') then
          FIN_TRI_NUM=CFG_VAL

        else
          write(*,'(A,A20)')
     *     'warning: unknown parameter in cfg-file: ', CFG_KEY
        endif
        goto 1101
 1102   close(54)
      end
