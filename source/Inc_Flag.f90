
      MODULE Inc_Flag

      IMPLICIT NONE

      LOGICAL(1) :: Flag_XSFB
      LOGICAL(1) :: Flag_THFB
      logical(1) :: flag_th_original=.true.
      logical(1) :: flag_th1d_pbp=.false.

      LOGICAL(1) :: Flag_InDET
      LOGICAL(1) :: Flag_ExDET

      LOGICAL(1) :: Flag_ARST = .false.

      LOGICAL(1) :: Flag_powhist = .false.
      LOGICAL(1) :: Flag_tmmhist= .false.

      LOGICAL(1) :: Flag_rod = .false.

      logical(1) :: flag_branch=.false.
      logical(1) :: flag_doing_branch=.false.
      logical(1) :: flag_tmfb=.true.
      logical(1) :: flag_tffb=.true.

      logical(1) :: flag_macroxs=.false.
      logical(1) :: flag_tr_macroxs=.false.
      logical(1) :: flag_tr_macroxs_inp=.false.
      logical(1) :: flag_tr_macroxs_1d=.true.
      logical(1) :: flag_microxs_1d=.false.
#ifdef tuan_fr
      logical(1) :: flag_miXS=.false.
#endif

#ifdef js_mpc
      logical(1) :: flag_pass_next = .false.
#endif

      logical(1) :: opt_pjumpin = .false.


      END MODULE Inc_Flag
