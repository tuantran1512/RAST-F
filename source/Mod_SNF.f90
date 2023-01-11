#ifdef siarhei_delete


      MODULE Mod_SNF
#ifdef JR_SNF
      USE Inc_SNF
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
      CONTAINS

      SUBROUTINE Write_SNF
      implicit none
#ifdef JR_SNF

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Write_SNF] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call GET_SNF
      return
#else
      return
#endif
      END SUBROUTINE Write_SNF

#ifdef JR_SNF
      SUBROUTINE GET_SNF
      IMPLICIT NONE
      INTEGER(4)         :: i, j
      INTEGER(4)         :: temp1


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [GET_SNF] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (.not.Flag_SNF) RETURN
      temp1 = 1
      IF (Flag_SNF_NFA) THEN
         DO i = 1, SNF_NFA
            SNF_PRINT_NN = SNF_PRINT_NN+SNF_NZ(i)
         ENDDO
      ENDIF
      Do i = 1, SNF_NFA
         IF (Flag_SNF_NFA) temp1 = (SNF_NZ(i))
         DO j = 1, temp1
            IF ((i==1) .and. (j==1)) THEN
                Flag_SNF_First = .TRUE.
            ELSE
                Flag_SNF_First = .FALSE.
            ENDIF
            IF (Flag_SNF_NFA) THEN
               SNF_Ixy_1N = SNF_Ixy_1N_M(i)
               SNF_K      = SNF_K_M(i,j)
               SNF_PRINT_INDEX = SNF_PRINT_INDEX+1
            ENDIF
            write(*,*) '============================================================'
            write(*,*) '++ SNF CALCULATION '
            if (flag_hri) then
               CALL SNF_HRI
            elseif (flag_hrst) THEN
               write(*,*) '03. PRINT HRST'
               CALL SNF_HRST
            ELSE
               write(*,*) '01. READ ND FILES '
               CALL SNF_READ
               write(*,*) '02. INTERPOLATION CALCULATION : LAGRANGE '
               CALL SNF_CAL
            ENDIF
         ENDDO
      ENDDO
      IF (flag_snf_decay) THEN
         call SNF_COOLING
      endif
      write(*,*) '============================================================'
      write(*,*) ' '

      RETURN
      END SUBROUTINE Get_SNF

      SUBROUTINE SNF_READ
      USE MOD_CHAREDIT,     ONLY: CHARTOUPPER
      USE Inc_Depletion,    only: Cycle_BU, N_BU
      USE Inc_History,      only: Hist_Burnup_FA
      USE Inc_Pinvar,       only: npin
      IMPLICIT NONE
      LOGICAL(1)                           :: FLAG_OPEN
      CHARACTER(15)                        :: SubTypeName
      CHARACTER(200)                       :: str1,str2
      CHARACTER(150)                        :: div_end
      REAL(8)                              :: temp1
      INTEGER(4)                           :: ii
      INTEGER(4)                           :: jj
      INTEGER(4)                           :: i
      INTEGER(4)                           :: mid_num
      INTEGER(4)                           :: mid_num_v
      INTEGER(4)                           :: mid_start
      INTEGER(4)                           :: mid_end
      REAL(8)                              :: temp2, temp3, temp4, temp5
      REAL(8)                              :: mid_index_num
      CHARACTER(200)                        :: str3,str4,str5
      CHARACTER(10)                        :: str7,str8,str9
      INTEGER(4)                           :: FLAG_MIDDLE_TMO_U,FLAG_MIDDLE_TMO_D
      INTEGER(4)                           :: FLAG_MIDDLE_TFU_U,FLAG_MIDDLE_TFU_D
      INTEGER(4)                           :: FLAG_MIDDLE_BOR_U,FLAG_MIDDLE_BOR_D
      div_end = '_______________________________________________________________&
              & ________________________________________'


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_READ] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      FLAG_MIDDLE_TMO_U = 0
      FLAG_MIDDLE_TMO_D = 0
      FLAG_MIDDLE_TFU_U = 0
      FLAG_MIDDLE_TFU_D = 0
      FLAG_MIDDLE_BOR_U = 0
      FLAG_MIDDLE_BOR_D = 0

      IF (SNF_N_BU==0) THEN
         SNF_N_END_CYBU = N_BU
      ELSE
         DO i = 1,N_BU
            IF(SNF_BU(SNF_N_BU) == Cycle_BU(i)) THEN
               SNF_N_END_CYBU = i
            END IF
         END DO
      endif

      ! [1-1] OPEN FILE ----------------------------------------- !
      DO WHILE (.TRUE.)
         INQUIRE( UNIT = R_ND, OPENED = Flag_Open )
         IF ( .NOT. Flag_Open ) THEN
            OPEN(R_ND, FILE = NAME_ND)
            EXIT
         ELSE
         R_ND = R_ND + 10
         END IF
      END DO

      ! INITIALIZATION OF BRANCH .ND ARRAYS
      IF (allocated(SNF_REF_ARRAY   )) deallocate(SNF_REF_ARRAY   )
      IF (allocated(SNF_TMO_ARRAY   )) deallocate(SNF_TMO_ARRAY   )
      IF (allocated(SNF_TFU_ARRAY   )) deallocate(SNF_TFU_ARRAY   )
      IF (allocated(SNF_BOR_ARRAY   )) deallocate(SNF_BOR_ARRAY   )
      IF (allocated(SNF_BRAN_BUHIST )) deallocate(SNF_BRAN_BUHIST )
      IF (allocated(SNF_BRAN_DAYHIST)) deallocate(SNF_BRAN_DAYHIST)
      !IF (allocated(SNF_PIN_POSITION)) deallocate(SNF_PIN_POSITION)
      !IF (allocated(SNF_PIN_RAD     )) deallocate(SNF_PIN_RAD)
      IF (allocated(SNF_LAMBDA_ARRAY)) deallocate(SNF_LAMBDA_ARRAY)
      IF (allocated(SNF_BRAN_FLUX   )) deallocate(SNF_BRAN_FLUX)
      IF (allocated(SNF_XS          )) deallocate(SNF_XS)
      !IF (allocated(SNF_BRAN_NDD    )) deallocate(SNF_BRAN_NDD)
      !IF (allocated(SNF_R2_NDD      )) deallocate(SNF_R2_NDD)
      SNF_BRAN_P    = 0d0
      mid_num       = 0
      mid_num_v     = 0
      mid_index_num = 0
      mid_start     = 0
      mid_end       = 0


      SubType: DO WHILE (.TRUE.)
666      READ(R_ND, *, END = 999) SubTypeName
         IF ((SubTypeName(1:5) .NE. 'BRANC') .and. (SubTypeName(1:3) .NE. 'For' )          &
              .and. (SubTypeName(1:5) .NE. 'TIT  ') .and. (SubTypeName(1:5) .NE. 'HIST|' ) &
              .and. (SubTypeName(1:5) .NE. 'Calcu') .and. (SubTypeName(1:5) .NE. 'Posit')  &
              .and. (SubTypeName(1:5) .NE. 'VOL|S') .and. (SubTypeName(1:5) .NE. 'NDEN2')) THEN
            CYCLE SubType
         END IF
         BACKSPACE(R_ND)
         CALL CharToUpper(SubTypeName)
         SELECT CASE (SubTypeName)
         CASE ('BRANCH')
            READ(R_ND,*,END=999) SubTypeName,str1,SNF_BRAN_OPT
         CASE ('FOR')
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) SubTypeName, SNF_BRAN_TFU
            READ(R_ND,*,END=999) SubTypeName, SNF_BRAN_TMO
            READ(R_ND,*,END=999) SubTypeName, SNF_BRAN_BOR
            READ(R_ND,*,END=999) str1, temp1
            READ(R_ND,*,END=999) SubTypeName, SNF_N_ISO
            READ(R_ND,*,END=999) SubTypeName, SNF_N_STEP
            READ(R_ND,*,END=999) SubTypeName, SNF_N_INDEXING ! Number of assembly
            READ(R_ND,*,END=999) SubTypeName, SNF_MAX_NR
            READ(R_ND,*,END=999) SubTypeName, SNF_MIN_MAT
            SNF_N_STEP = SNF_N_STEP+1  ! ADD ZERO STEP
         CASE ('POSITION')
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            IF (allocated(SNF_PIN_POSITION)) deallocate(SNF_PIN_POSITION)
            allocate(SNF_PIN_POSITION(int(npin)**2,2)) ! 1: ix, 2: iy
            SNF_PIN_POSITION = 0
            IF (allocated(SNF_PIN_MAT)     ) deallocate(SNF_PIN_MAT)
            allocate(SNF_PIN_MAT(int(npin)**2,SNF_MAX_NR))
            SNF_PIN_MAT = 0
            DO ii = 1,(int(npin)**2)
               READ(R_ND,*,END=999) temp1
               BACKSPACE(R_ND)
               READ(R_ND,*,END=999) temp2, temp3, temp4, temp5, SNF_PIN_POSITION(ii,1), SNF_PIN_POSITION(ii,2), &
                  SNF_PIN_MAT(ii,1:int(temp1))
            ENDDO
            SNF_MAX_MAT = maxval(SNF_PIN_MAT(:,1))
            SNF_MAX_MAT = int(npin**2d0)
            IF (.not.Flag_snf_pin) then
               SNF_MIN_MAT = 1
               SNF_MAX_MAT = 1
            endif
            !initialization
            IF (allocated(SNF_BRAN_BUHIST )) deallocate(SNF_BRAN_BUHIST )
            IF (allocated(SNF_BRAN_DAYHIST)) deallocate(SNF_BRAN_DAYHIST)
            allocate(SNF_BRAN_BUHIST (SNF_N_STEP))
            allocate(SNF_BRAN_DAYHIST(SNF_N_STEP))
            SNF_BRAN_BUHIST  = 0d0
            SNF_BRAN_DAYHIST = 0d0
            IF (.not.FLAG_SNF_PIN) SNF_N_INDEXING = 1    ! CHECK CALCULATION OPTION
            IF (allocated(SNF_ISO_ARRAY   )) deallocate(SNF_ISO_ARRAY)
            IF (allocated(SNF_REF_ARRAY   )) deallocate(SNF_REF_ARRAY)
            IF (allocated(SNF_TMO_ARRAY   )) deallocate(SNF_TMO_ARRAY)
            IF (allocated(SNF_TFU_ARRAY   )) deallocate(SNF_TFU_ARRAY)
            IF (allocated(SNF_BOR_ARRAY   )) deallocate(SNF_BOR_ARRAY)
            IF (allocated(SNF_LAMBDA_ARRAY)) deallocate(SNF_LAMBDA_ARRAY)
            IF (allocated(SNF_BRAN_FLUX   )) deallocate(SNF_BRAN_FLUX)
            IF (allocated(SNF_XS          )) deallocate(SNF_XS)
            IF (allocated(SNF_VOL_RATIO   )) deallocate(SNF_VOL_RATIO)
            IF (allocated(HNDEN2GCC       )) deallocate(HNDEN2GCC)
            allocate(SNF_ISO_ARRAY(SNF_N_ISO             ))
            allocate(SNF_REF_ARRAY(SNF_N_INDEXING,SNF_N_STEP,SNF_N_ISO  ))  ! ARRAY
            allocate(SNF_TMO_ARRAY(SNF_N_INDEXING,2,SNF_N_ISO,4))           ! ARRAY
            allocate(SNF_TFU_ARRAY(SNF_N_INDEXING,2,SNF_N_ISO,2))           ! ARRAY
            allocate(SNF_BOR_ARRAY(SNF_N_INDEXING,2,SNF_N_ISO,3))           ! ARRAY
            allocate(SNF_LAMBDA_ARRAY(SNF_N_ISO))
            allocate(SNF_BRAN_FLUX(SNF_N_INDEXING,SNF_N_STEP,SNF_N_ISO))
            allocate(SNF_XS(SNF_N_ISO))
            allocate(SNF_VOL_RATIO(SNF_N_INDEXING,SNF_N_STEP))
            allocate(HNDEN2GCC(SNF_N_ISO))
            SNF_REF_ARRAY = 0d0
            SNF_TMO_ARRAY = 0d0
            SNF_TFU_ARRAY = 0d0
            SNF_BOR_ARRAY = 0d0
            SNF_BRAN_FLUX = 0d0
            SNF_VOL_RATIO = 0d0
            go to 666
         CASE ('CALCULATION')
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            FLAG_SNF_BU_INTER = .TRUE.
            DO ii = 2, (SNF_N_STEP)
               READ(R_ND,*,END=999) temp1, SNF_BRAN_BUHIST(ii), SNF_BRAN_DAYHIST(ii), SNF_BRAN_P
               write(str3,'(f7.3)') SNF_BRAN_BUHIST(ii)
               write(str4,'(f7.3)') HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
               write(str5,'(f7.3)') SNF_BRAN_BUHIST(ii-1)
               IF (str3 == str4) THEN
                  FLAG_SNF_BU_INTER = .FALSE.
                  SNF_N_BU_UP   = (ii)
                  SNF_N_BU_DOWN = (ii)
               ELSE IF ((str5<str4) .and. (str3>str4)) THEN
                  SNF_N_BU_UP   = (ii)
                  SNF_N_BU_DOWN = (ii-1)
               END IF
            ENDDO
         CASE ('NDEN2GCC')
            READ(R_ND,*,END=999) str1
            DO ii = 1,SNF_N_ISO
               READ(R_ND,*,END=999) temp1, HNDEN2GCC(ii)
            ENDDO
         CASE ('VOL|SUM')
            READ(R_ND,*,END=999) str1
            mid_num_v=mid_num_v+1
            if (mid_num_v < (SNF_N_STEP+1)) then
               IF (flag_snf_pin) then
                  DO jj = SNF_MIN_MAT,SNF_MAX_MAT
                     READ(R_ND,*,END=999) temp1, SNF_VOL_RATIO(jj,mid_num_v)
                  ENDDO
               else
                  READ(R_ND,*,END=999) temp1, SNF_VOL_RATIO(1,mid_num_v)
               endif
            endif
         CASE ('TIT')
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) temp1, str1
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            READ(R_ND,*,END=999) str1
            mid_num = mid_num+1
            if (mid_num < (SNF_N_STEP+1)) then
               if (mid_num == 1) then
                  DO jj = SNF_MIN_MAT,SNF_MAX_MAT
                     DO ii = 1,SNF_N_ISO
                        READ(R_ND,*,END=999) SNF_ISO_ARRAY(ii),SNF_REF_ARRAY(jj,mid_num,ii),SNF_LAMBDA_ARRAY(ii),SNF_BRAN_FLUX(jj,mid_num,ii) ,SNF_XS(ii)
                     ENDDO
                  ENDDO
               else
                  DO jj = SNF_MIN_MAT,SNF_MAX_MAT
                     DO ii = 1,SNF_N_ISO
                        READ(R_ND,*,END=999) temp1,SNF_REF_ARRAY(jj,mid_num,ii),str1,SNF_BRAN_FLUX(jj,mid_num,ii),str2
                     ENDDO
                  ENDDO
               END IF
            endif
         CASE ('HIST|COE')
            BACKSPACE(R_ND)
            READ(R_ND, *) str1, str2
            IF     (str2 == 'TMO') THEN
               temp1 = 0d0
               READ(R_ND, * ,END=999) str1
               READ(R_ND, * ,END=999) temp1, str2
               IF (FLAG_SNF_BU_INTER) THEN
                  WRITE(str7,'(f7.3)') temp1
                  WRITE(str8,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  WRITE(str9,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_DOWN)
                  IF ( str7 == str8 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TMO_U=FLAG_MIDDLE_TMO_U+1
                     IF(FLAG_MIDDLE_TMO_U>4) go to 666
                     DO jj = SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TMO_ARRAY(jj,2,ii,FLAG_MIDDLE_TMO_U)
                        ENDDO
                     ENDDO
                  ELSEIF ( str7 == str9 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TMO_D=FLAG_MIDDLE_TMO_D+1
                     IF(FLAG_MIDDLE_TMO_D>4) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TMO_ARRAY(jj,1,ii,FLAG_MIDDLE_TMO_D)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE
                  WRITE(str1,'(f7.3)') temp1
                  WRITE(str2,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  IF ( str1 == str2 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TMO_U=FLAG_MIDDLE_TMO_U+1
                     IF(FLAG_MIDDLE_TMO_U>4) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TMO_ARRAY(jj,2,ii,FLAG_MIDDLE_TMO_U)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ELSEIF (str2 == 'BOR') THEN
               READ(R_ND, * ,END=999) str1
               READ(R_ND, * ,END=999) temp1, str2
               IF (FLAG_SNF_BU_INTER) THEN
                  WRITE(str1,'(f7.3)') temp1
                  WRITE(str2,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  WRITE(str3,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_DOWN)
                  IF ( str1 == str2 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_BOR_U=FLAG_MIDDLE_BOR_U+1
                     IF(FLAG_MIDDLE_BOR_U>3) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_BOR_ARRAY(jj,2,ii,FLAG_MIDDLE_BOR_U)
                        ENDDO
                     ENDDO
                  ELSEIF ( str1 == str3 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_BOR_D=FLAG_MIDDLE_BOR_D+1
                     IF(FLAG_MIDDLE_BOR_D>3) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_BOR_ARRAY(jj,1,ii,FLAG_MIDDLE_BOR_D)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE
                  WRITE(str1,'(f7.3)') temp1
                  WRITE(str2,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  IF ( str1 == str2 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_BOR_U=FLAG_MIDDLE_BOR_U+1
                     IF(FLAG_MIDDLE_BOR_U>3) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_BOR_ARRAY(jj,2,ii,FLAG_MIDDLE_BOR_U)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ELSEIF (str2 == 'TFU') THEN
               READ(R_ND, * ,END=999) str1
               READ(R_ND, * ,END=999) temp1, str2
               IF (FLAG_SNF_BU_INTER) THEN
                  WRITE(str1,'(f7.3)') temp1
                  WRITE(str2,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  WRITE(str3,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_DOWN)
                  IF ( str1 == str2 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TFU_U=FLAG_MIDDLE_TFU_U+1
                     IF(FLAG_MIDDLE_TFU_U>2) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TFU_ARRAY(jj,2,ii,FLAG_MIDDLE_TFU_U)
                        ENDDO
                     ENDDO
                  ELSEIF ( str1 == str3 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TFU_D=FLAG_MIDDLE_TFU_D+1
                     IF(FLAG_MIDDLE_TFU_D>2) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TFU_ARRAY(jj,1,ii,FLAG_MIDDLE_TFU_D)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE
                  WRITE(str1,'(f7.3)') temp1
                  WRITE(str2,'(f7.3)') SNF_BRAN_BUHIST(SNF_N_BU_UP)
                  IF ( str1 == str2 ) THEN
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     READ(R_ND, * ,END=999) str1
                     FLAG_MIDDLE_TFU_U=FLAG_MIDDLE_TFU_U+1
                     IF(FLAG_MIDDLE_TFU_U>2) go to 666
                     DO jj=SNF_MIN_MAT,SNF_MAX_MAT
                        DO ii = 1,SNF_N_ISO
                           READ(R_ND,*,END=999) temp1, SNF_TFU_ARRAY(jj,2,ii,FLAG_MIDDLE_TFU_U)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            END IF
         END SELECT
      END DO SubType

999   continue

      CLOSE(R_ND)

      RETURN
      END SUBROUTINE SNF_READ


      SUBROUTINE SNF_CAL
      USE Inc_History
      USE Inc_Utile,        only: Blank_Line
      USE Inc_Geometry
      USE Inc_Pinvar,       only: npin
      USE Inc_FA,           only: N_pin
#ifdef JR_SRCTRM
      USE MOD_END_SRC
!      USE Mod_SRC,          only: GET_SRCTRM
#endif
      IMPLICIT NONE
      CHARACTER(80)                        :: FMT
      INTEGER(4)                           :: pin_count
      INTEGER(4)                           :: temp_int
      INTEGER(4)                           :: Iz
      INTEGER(4)                           :: iy_1n
      INTEGER                              :: ii
      INTEGER                              :: jj
      INTEGER                              :: kk
      INTEGER                              :: rr
      INTEGER                              :: ip,jp
      CHARACTER(15)                        :: str1, str2
      INTEGER(4),DIMENSION(41)             :: R2_ISO_LIST
      INTEGER(4),DIMENSION(41)             :: R2_ISO_MATCH
      INTEGER(4)                           :: mid_iso
      INTEGER(4),DIMENSION(:),ALLOCATABLE  :: NN_MAT
      ! index parameters
      REAL(8)                              :: SUM_TIME_X_HIST
      REAL(8)                              :: dBU
      ! interpolation parameters
      REAL(8)                              :: L0, L1, L2, L3, L4
      REAL(8)                              :: tmo1, tmo2, tmo3, tmo4, tmo5, tmp_ref
      REAL(8)                              :: SUM_Middle
      REAL(8)                              :: tfu1, tfu2, tfu3
      REAL(8)                              :: bor1, bor2, bor3, bor4
      REAL(8)                              :: temp1, temp2, temp3, temp4, temp5
      REAL(8)                              :: mid
      REAL(8)                              :: L_UPPER, L_LOWER
#ifdef JR_SNF_INTERPOLATION
      REAL(8)                              :: B0, B1, B2, B3, B4
#endif
      ! interpolation results
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE     :: CAL_TMO_RESULT
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE     :: CAL_TFU_RESULT
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE     :: CAL_BOR_RESULT
      REAL(8),DIMENSION(:,:),ALLOCATABLE       :: SNF_RESULT_ARRAY_OLD


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_CAL] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_snf_pin) then
         IF (allocated(NN_MAT)) deallocate(NN_MAT)
         allocate(NN_MAT(SNF_MAX_MAT))
         NN_MAT = 0
         DO rr = 1, (int(npin)**2)
            !kk = SNF_PIN_MAT(rr,1)
            kk = rr
            NN_MAT(kk) = NN_MAT(kk)+1
         ENDDO
         DO kk = SNF_MIN_MAT,SNF_MAX_MAT
            IF (NN_MAT(kk)==0) NN_MAT(kk)=1
            DO jj = 1,2
               DO ii = 1,SNF_N_ISO
                  DO rr = 1, 4
                     SNF_TMO_ARRAY(kk,jj,ii,rr) = SNF_TMO_ARRAY(kk,jj,ii,rr)/real(NN_MAT(kk))
                  ENDDO
                  DO rr = 1,2
                     SNF_TFU_ARRAY(kk,jj,ii,rr) = SNF_TFU_ARRAY(kk,jj,ii,rr)/real(NN_MAT(kk))
                  ENDDO
                  DO rr = 1,3
                     SNF_BOR_ARRAY(kk,jj,ii,rr) = SNF_BOR_ARRAY(kk,jj,ii,rr)/real(NN_MAT(kk))
                  ENDDO
               ENDDO
            ENDDO
            DO jj = 1,SNF_N_STEP
               DO ii = 1,SNF_N_ISO
                  SNF_REF_ARRAY(kk,jj,ii)   = SNF_REF_ARRAY(kk,jj,ii)  /real(NN_MAT(kk))
               ENDDO
            ENDDO
         ENDDO
      endif
      ! 01. R2 HISTORY INDEX
      ! 01-1) TMO HISTORY INDEX
      SUM_TIME_X_HIST = 0d0
      IF (.not.flag_snf_pin .or. flag_hri) then
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_HRI_N
            dBU = SNF_HRI_BU(jj)-SNF_HRI_BU(jj-1)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*SNF_HRI_TMO(jj)
         ENDDO
         dBU = HIST_BURNUP_FA(1,SNF_Ixy_1N,SNF_K)-SNF_HRI_BU(SNF_HRI_N)
         SUM_TIME_X_HIST=SUM_TIME_X_HIST+Hist_T_Mod_FA(1,SNF_Ixy_1N,SNF_K)
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*Hist_T_Mod_FA(jj-1,SNF_Ixy_1N,SNF_K)
         ENDDO
         SNF_TMO_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
         ! 01-2) TFU HISTORY INDEX
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_HRI_N
            dBU = SNF_HRI_BU(jj)-SNF_HRI_BU(jj-1)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*SNF_HRI_TFU(jj)
         ENDDO
         dBU = HIST_BURNUP_FA(1,SNF_Ixy_1N,SNF_K)-SNF_HRI_BU(SNF_HRI_N)
         SUM_TIME_X_HIST=SUM_TIME_X_HIST+Hist_T_Fuel_FA(1,SNF_Ixy_1N,SNF_K)
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*Hist_T_Fuel_FA(jj-1,SNF_Ixy_1N,SNF_K)
         ENDDO
         SNF_TFU_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
         ! 01-3) BOR HISTORY INDEX
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_HRI_N
            dBU = SNF_HRI_BU(jj)-SNF_HRI_BU(jj-1)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*SNF_HRI_BOR(jj)
         ENDDO
         dBU = HIST_BURNUP_FA(1,SNF_Ixy_1N,SNF_K)-SNF_HRI_BU(SNF_HRI_N)
         SUM_TIME_X_HIST=SUM_TIME_X_HIST+hist_PPM(1)
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*hist_PPM(jj-1)
         ENDDO
         SNF_BOR_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
      else
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*Hist_T_Mod_FA(jj-1,SNF_Ixy_1N,SNF_K)
         ENDDO
         SNF_TMO_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
         ! 01-2) TFU HISTORY INDEX
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*Hist_T_Fuel_FA(jj-1,SNF_Ixy_1N,SNF_K)
         ENDDO
         SNF_TFU_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
         ! 01-3) BOR HISTORY INDEX
         SUM_TIME_X_HIST = 0d0
         DO jj = 2,SNF_N_END_CYBU
            dBU = HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K)-HIST_BURNUP_FA(jj-1,SNF_Ixy_1N,SNF_K)
            SUM_TIME_X_HIST=SUM_TIME_X_HIST+dBU*hist_PPM(jj-1)
         ENDDO
         SNF_BOR_INDEX = SUM_TIME_X_HIST/HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
      endif

      write(*,'(A,F10.2)') '     SNF_TMO_INDEX',SNF_TMO_INDEX
      write(*,'(A,F10.2)') '     SNF_TFU_INDEX',SNF_TFU_INDEX
      write(*,'(A,F10.2)') '     SNF_BOR_INDEX',SNF_BOR_INDEX
      ! 02. LAGRANGE INTERPOLATION
      ! Initialization
      IF (.not.allocated(CAL_TMO_RESULT)) allocate(CAL_TMO_RESULT(SNF_MAX_MAT,2,SNF_N_ISO)) ! 1: UP, 2: INTERPOLATION
      IF (.not.allocated(CAL_TFU_RESULT)) allocate(CAL_TFU_RESULT(SNF_MAX_MAT,2,SNF_N_ISO)) ! 1: UP, 2: INTERPOLATION
      IF (.not.allocated(CAL_BOR_RESULT)) allocate(CAL_BOR_RESULT(SNF_MAX_MAT,2,SNF_N_ISO)) ! 1: UP, 2: INTERPOLATION
      CAL_TMO_RESULT = 0d0
      CAL_TFU_RESULT = 0d0
      CAL_BOR_RESULT = 0d0
      IF(FLAG_SNF_BU_INTER) WRITE(*,'(A17)') 'Interpolation: T'
      ! 02-1) TMO
      WRITE(str1,'(f7.1)') SNF_BRAN_TMO
      WRITE(str2,'(f7.1)') SNF_TMO_INDEX
      IF (str1 .NE. str2) THEN
         write(*,*) '    TMO INTERPOLATION'
#ifdef JR_SNF_INTERPOLATION
         if (flag_calculation) then
#endif
#ifdef JR_SNF_TMO
         if (Flag_tmo_order==4) then
#endif
         tmo1 = (SNF_BRAN_TMO-40d0)
         tmo2 = (SNF_BRAN_TMO-20d0)
         tmo3 = (SNF_BRAN_TMO)
         tmo4 = (SNF_BRAN_TMO+20d0)
         tmo5 = (SNF_BRAN_TMO+40d0)
         tmp_ref = (SNF_TMO_INDEX)

         L_UPPER = (tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)*(tmp_ref-tmo5)
         L_LOWER = (tmo1-tmo2)*(tmo1-tmo3)*(tmo1-tmo4)*(tmo1-tmo5)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo3)*(tmp_ref-tmo4)*(tmp_ref-tmo5)
         L_LOWER = (tmo2-tmo1)*(tmo2-tmo3)*(tmo2-tmo4)*(tmo2-tmo5)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo4)*(tmp_ref-tmo5)
         L_LOWER = (tmo3-tmo1)*(tmo3-tmo2)*(tmo3-tmo4)*(tmo3-tmo5)
         L2      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo5)
         L_LOWER = (tmo4-tmo1)*(tmo4-tmo2)*(tmo4-tmo3)*(tmo4-tmo5)
         L3      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
         L_LOWER = (tmo5-tmo1)*(tmo5-tmo2)*(tmo5-tmo3)*(tmo5-tmo4)
         L4      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_TMO_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_TMO_ARRAY(jj, 2, ii, 2)
               temp3 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp4 = SNF_TMO_ARRAY(jj, 2, ii, 3)
               temp5 = SNF_TMO_ARRAY(jj, 2, ii, 4)
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               SUM_Middle = SUM_Middle+L3*temp4+L4*temp5
               CAL_TMO_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 1) - SNF_TMO_ARRAY(jj, 1, ii, 1)) + SNF_TMO_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 2) - SNF_TMO_ARRAY(jj, 1, ii, 2)) + SNF_TMO_ARRAY(jj, 1, ii, 2)
                  temp3 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp4 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 3) - SNF_TMO_ARRAY(jj, 1, ii, 3)) + SNF_TMO_ARRAY(jj, 1, ii, 3)
                  temp5 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 4) - SNF_TMO_ARRAY(jj, 1, ii, 4)) + SNF_TMO_ARRAY(jj, 1, ii, 4)
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3
                  SUM_Middle = SUM_Middle+L3*temp4+L4*temp5
                  CAL_TMO_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
#ifdef JR_SNF_TMO
         else if (Flag_tmo_order == 2)then
         write(*,*) "ORDER 2"
         tmo1 = (SNF_BRAN_TMO-20d0)
         tmo2 = (SNF_BRAN_TMO)
         tmo3 = (SNF_BRAN_TMO+20d0)
         tmp_ref = (SNF_TMO_INDEX)

         L_UPPER = (tmp_ref-tmo2)*(tmp_ref-tmo3)
         L_LOWER = (tmo1-tmo2)*(tmo1-tmo3)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo3)
         L_LOWER = (tmo2-tmo1)*(tmo2-tmo3)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)
         L_LOWER = (tmo3-tmo1)*(tmo3-tmo2)
         L2      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_TMO_ARRAY(jj, 2, ii, 2)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_TMO_ARRAY(jj, 2, ii, 3)
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               CAL_TMO_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 2) - SNF_TMO_ARRAY(jj, 1, ii, 2)) + SNF_TMO_ARRAY(jj, 1, ii, 2)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 3) - SNF_TMO_ARRAY(jj, 1, ii, 3)) + SNF_TMO_ARRAY(jj, 1, ii, 3)
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3
                  CAL_TMO_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
         elseif (Flag_tmo_order == 3) then
         write(*,*) "ORDER 3"
         if (tmp_ref < ((SNF_BRAN_TMO+20d0))) then
         tmo1 = (SNF_BRAN_TMO-40d0)
         tmo2 = (SNF_BRAN_TMO-20d0)
         tmo3 = (SNF_BRAN_TMO)
         tmo4 = (SNF_BRAN_TMO+20d0)
         tmp_ref = (SNF_TMO_INDEX)

         L_UPPER = (tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
         L_LOWER = (tmo1-tmo2)*(tmo1-tmo3)*(tmo1-tmo4)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
         L_LOWER = (tmo2-tmo1)*(tmo2-tmo3)*(tmo2-tmo4)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo4)
         L_LOWER = (tmo3-tmo1)*(tmo3-tmo2)*(tmo3-tmo4)
         L2      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)
         L_LOWER = (tmo4-tmo1)*(tmo4-tmo2)*(tmo4-tmo3)
         L3      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_TMO_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_TMO_ARRAY(jj, 2, ii, 2)
               temp3 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp4 = SNF_TMO_ARRAY(jj, 2, ii, 3)
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               SUM_Middle = SUM_Middle+L3*temp4
               CAL_TMO_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 1) - SNF_TMO_ARRAY(jj, 1, ii, 1)) + SNF_TMO_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 2) - SNF_TMO_ARRAY(jj, 1, ii, 2)) + SNF_TMO_ARRAY(jj, 1, ii, 2)
                  temp3 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp4 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 3) - SNF_TMO_ARRAY(jj, 1, ii, 3)) + SNF_TMO_ARRAY(jj, 1, ii, 3)
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3
                  SUM_Middle = SUM_Middle+L3*temp4
                  CAL_TMO_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
         endif
         else
         tmo1 = (SNF_BRAN_TMO-20d0)
         tmo2 = (SNF_BRAN_TMO)
         tmo3 = (SNF_BRAN_TMO+20d0)
         tmo4 = (SNF_BRAN_TMO+40d0)
         tmp_ref = (SNF_TMO_INDEX)

         L_UPPER = (tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
         L_LOWER = (tmo1-tmo2)*(tmo1-tmo3)*(tmo1-tmo4)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
         L_LOWER = (tmo2-tmo1)*(tmo2-tmo3)*(tmo2-tmo4)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo4)
         L_LOWER = (tmo3-tmo1)*(tmo3-tmo2)*(tmo3-tmo4)
         L2      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)
         L_LOWER = (tmo4-tmo1)*(tmo4-tmo2)*(tmo4-tmo3)
         L3      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_TMO_ARRAY(jj, 2, ii, 2)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_TMO_ARRAY(jj, 2, ii, 3)
               temp4 = SNF_TMO_ARRAY(jj, 2, ii, 4)
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               SUM_Middle = SUM_Middle+L3*temp4
               CAL_TMO_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 2) - SNF_TMO_ARRAY(jj, 1, ii, 2)) + SNF_TMO_ARRAY(jj, 1, ii, 2)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 3) - SNF_TMO_ARRAY(jj, 1, ii, 3)) + SNF_TMO_ARRAY(jj, 1, ii, 3)
                  temp4 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 4) - SNF_TMO_ARRAY(jj, 1, ii, 4)) + SNF_TMO_ARRAY(jj, 1, ii, 4)
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3
                  SUM_Middle = SUM_Middle+L3*temp4
                  CAL_TMO_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
         endif
#endif
#ifdef JR_SNF_INTERPOLATION
         else
         write(*,*) '    TMO INTERPOLATION >> Newton-Raphson polynomial'
         tmo1 = (SNF_BRAN_TMO-40d0)
         tmo2 = (SNF_BRAN_TMO-20d0)
         tmo3 = (SNF_BRAN_TMO)
         tmo4 = (SNF_BRAN_TMO+20d0)
         tmo5 = (SNF_BRAN_TMO+40d0)
         tmp_ref = (SNF_TMO_INDEX)
         DO jj = SNF_MIN_MAT, SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               temp1 = SNF_TMO_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_TMO_ARRAY(jj, 2, ii, 2)
               temp3 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp4 = SNF_TMO_ARRAY(jj, 2, ii, 3)
               temp5 = SNF_TMO_ARRAY(jj, 2, ii, 4)
               B0    = temp1
               B1    = (temp2 - temp1)/20d0
               B2    = (temp3 - B0 - B1 * (40d0))/(40d0*20d0)
               B3    = (temp4 - B0 - B1 * (60d0) - B2 * (60d0 * 40d0))/(60d0*40d0*20d0)
               B4    = (temp5 - B0 - B1 * (80d0) - B2 * (80d0 * 60d0) - B3 * (80d0 * 60d0 * 40d0))/(80d0*60d0*40d0*20d0)
               SUM_Middle = 0d0
               SUM_Middle = B0 + B1*(tmp_ref-tmo1) + B2*(tmp_ref-tmo1)*(tmp_ref-tmo2)
               SUM_Middle = SUM_Middle + B3*(tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)
               SUM_Middle = SUM_Middle + B4*(tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
               CAL_TMO_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 1) - SNF_TMO_ARRAY(jj, 1, ii, 1)) + SNF_TMO_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 2) - SNF_TMO_ARRAY(jj, 1, ii, 2)) + SNF_TMO_ARRAY(jj, 1, ii, 2)
                  temp3 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp4 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 3) - SNF_TMO_ARRAY(jj, 1, ii, 3)) + SNF_TMO_ARRAY(jj, 1, ii, 3)
                  temp5 = mid*(SNF_TMO_ARRAY(jj, 2, ii, 4) - SNF_TMO_ARRAY(jj, 1, ii, 4)) + SNF_TMO_ARRAY(jj, 1, ii, 4)
                  B0    = temp1
                  B1    = (temp2 - temp1)/20d0
                  B2    = (temp3 - B0 - B1 * (40d0))/(40d0*20d0)
                  B3    = (temp4 - B0 - B1 * (60d0) - B2 * (60d0 * 40d0))/(60d0*40d0*20d0)
                  B4    = (temp5 - B0 - B1 * (80d0) - B2 * (80d0 * 60d0) - B3 * (80d0 * 60d0 * 40d0))/(80d0*60d0*40d0*20d0)
                  SUM_Middle = 0d0
                  SUM_Middle = B0 + B1*(tmp_ref-tmo1) + B2*(tmp_ref-tmo1)*(tmp_ref-tmo2)
                  SUM_Middle = SUM_Middle + B3*(tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)
                  SUM_Middle = SUM_Middle + B4*(tmp_ref-tmo1)*(tmp_ref-tmo2)*(tmp_ref-tmo3)*(tmp_ref-tmo4)
                  CAL_TMO_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO
         ENDDO
         endif
#endif
      ELSE
         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               CAL_TMO_RESULT(jj,1,ii) = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               IF (FLAG_SNF_BU_INTER) THEN
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  CAL_TMO_RESULT(jj,2,ii) = temp1
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
      ENDIF

      ! 02-2) TFU
      WRITE(str1,'(f7.1)') SNF_BRAN_TFU
      WRITE(str2,'(f7.1)') SNF_TFU_INDEX
      IF (str1 .NE. str2) THEN
         write(*,*) '    TFU INTERPOLATION'
#ifdef JR_SNF_INTERPOLATION
         if (flag_calculation) then
#endif
         tfu1    = (SNF_BRAN_TMO-20d0)
         tfu2    = (SNF_BRAN_TFU)
         tfu3    = (1500d0)
         tmp_ref = (SNF_TFU_INDEX)

         L_UPPER = (tmp_ref-tfu2)*(tmp_ref-tfu3)
         L_LOWER = (tfu1-tfu2)*(tfu1-tfu3)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tfu1)*(tmp_ref-tfu3)
         L_LOWER = (tfu2-tfu1)*(tfu2-tfu3)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-tfu1)*(tmp_ref-tfu2)
         L_LOWER = (tfu3-tfu1)*(tfu3-tfu2)
         L2      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_TFU_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_TFU_ARRAY(jj, 2, ii, 2)
               SUM_Middle = 0
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               CAL_TFU_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TFU_ARRAY(jj, 2, ii, 1)-SNF_TFU_ARRAY(jj, 1, ii, 1))+SNF_TFU_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)-SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_TFU_ARRAY(jj, 2, ii, 2)-SNF_TFU_ARRAY(jj, 1, ii, 2))+SNF_TFU_ARRAY(jj, 1, ii, 2)
                  SUM_Middle = 0
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3
                  CAL_TFU_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
#ifdef JR_SNF_INTERPOLATION
         else
         write(*,*) '    TFU INTERPOLATION >> Newton-Raphson polynomial'
         tfu1    = (SNF_BRAN_TMO-20d0)
         tfu2    = (SNF_BRAN_TFU)
         tfu3    = (1500d0)
         tmp_ref = (SNF_TFU_INDEX)
         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               temp1 = SNF_TFU_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_TFU_ARRAY(jj, 2, ii, 2)
               B0    = temp1
               B1    = (temp2 - temp1)/(tfu2-tfu1)
               B2    = (temp3 - B0 - B1 * (tfu3-tfu1))/((tfu3-tfu1)*(tfu3-tfu2))
               SUM_Middle = 0d0
               SUM_Middle = B0 + B1*(tmp_ref-TFU1) + B2*(tmp_ref-TFU1)*(tmp_ref-TFU2)
               CAL_TFU_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_TFU_ARRAY(jj, 2, ii, 1) - SNF_TFU_ARRAY(jj, 1, ii, 1)) + SNF_TFU_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_TFU_ARRAY(jj, 2, ii, 2) - SNF_TFU_ARRAY(jj, 1, ii, 2)) + SNF_TFU_ARRAY(jj, 1, ii, 2)
                  B0    = temp1
                  B1    = (temp2 - temp1)/(tfu2-tfu1)
                  B2    = (temp3 - B0 - B1 * (tfu3-tfu1))/((tfu3-tfu1)*(tfu3-tfu2))
                  SUM_Middle = 0d0
                  SUM_Middle = B0 + B1*(tmp_ref-TFU1) + B2*(tmp_ref-TFU1)*(tmp_ref-TFU2)
                  CAL_TFU_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO
         ENDDO
         endif
#endif
      ELSE
         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               CAL_TFU_RESULT(jj,1,ii) = SNF_REF_ARRAY(jj,SNF_N_BU_UP, ii)
               IF (FLAG_SNF_BU_INTER) THEN
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  CAL_TFU_RESULT(jj,2,ii) = temp1
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
      ENDIF

      ! 02-3) BOR
      WRITE(str1,'(f7.1)') SNF_BRAN_BOR
      WRITE(str2,'(f7.1)') SNF_BOR_INDEX
      IF (str1 .NE. str2) THEN
         write(*,*) '    BOR INTERPOLATION'
#ifdef JR_SNF_INTERPOLATION
         if (flag_calculation) then
#endif
         bor1    = 0.1d0
         bor2    = (SNF_BRAN_BOR)
         bor3    = (SNF_BRAN_BOR*2d0)
         bor4    = (2400d0)
         tmp_ref = (SNF_BOR_INDEX)

         L_UPPER = (tmp_ref-bor2)*(tmp_ref-bor3)*(tmp_ref-bor4)
         L_LOWER = (bor1-bor2)*(bor1-bor3)*(bor1-bor4)
         L0      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-bor1)*(tmp_ref-bor3)*(tmp_ref-bor4)
         L_LOWER = (bor2-bor1)*(bor2-bor3)*(bor2-bor4)
         L1      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-bor1)*(tmp_ref-bor2)*(tmp_ref-bor4)
         L_LOWER = (bor3-bor1)*(bor3-bor2)*(bor3-bor4)
         L2      = (L_UPPER/L_LOWER)

         L_UPPER = (tmp_ref-bor1)*(tmp_ref-bor2)*(tmp_ref-bor3)
         L_LOWER = (bor4-bor1)*(bor4-bor2)*(bor4-bor3)
         L3      = (L_UPPER/L_LOWER)

         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               SUM_Middle = 0
               temp1 = SNF_BOR_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_BOR_ARRAY(jj, 2, ii, 2)
               temp4 = SNF_BOR_ARRAY(jj, 2, ii, 3)
               SUM_Middle = L0*temp1+L1*temp2+L2*temp3
               SUM_Middle = SUM_Middle+L3*temp4
               CAL_BOR_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 1)-SNF_BOR_ARRAY(jj, 1, ii, 1))+SNF_BOR_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)-SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 2)-SNF_BOR_ARRAY(jj, 1, ii, 2))+SNF_BOR_ARRAY(jj, 1, ii, 2)
                  temp4 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 3)-SNF_BOR_ARRAY(jj, 1, ii, 3))+SNF_BOR_ARRAY(jj, 1, ii, 3)
                  SUM_Middle = 0
                  SUM_Middle = L0*temp1+L1*temp2+L2*temp3+L3*temp4
                  CAL_BOR_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
#ifdef JR_SNF_INTERPOLATION
         else
         write(*,*) '    BOR INTERPOLATION >> Newton-Raphson polynomial'
         bor1    = 0.1d0
         bor2    = (SNF_BRAN_BOR)
         bor3    = (SNF_BRAN_BOR*2d0)
         bor4    = (2400d0)
         tmp_ref = (SNF_BOR_INDEX)
         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               temp1 = SNF_BOR_ARRAY(jj, 2, ii, 1)
               temp2 = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               temp3 = SNF_BOR_ARRAY(jj, 2, ii, 2)
               temp4 = SNF_BOR_ARRAY(jj, 2, ii, 3)
               B0    = temp1
               B1    = (temp2 - temp1)/(bor2-bor1)
               B2    = (temp3 - B0 - B1 * (bor3-bor1))/((bor3-bor1)*(bor3-bor2))
               B3    = (temp4 - B0 - B1 * (bor4-bor1) - B2*(bor4-bor1)*(bor4-bor2))/((bor4-bor1)*(bor4-bor2)*(bor4-bor3))
               SUM_Middle = 0d0
               SUM_Middle = B0 + B1*(tmp_ref-bor1) + B2*(tmp_ref-bor1)*(tmp_ref-bor2)
               SUM_Middle = SUM_Middle + B3*(tmp_ref-bor1)*(tmp_ref-bor2)*(tmp_ref-bor3)
               CAL_BOR_RESULT(jj,1,ii) = SUM_Middle
               IF (FLAG_SNF_BU_INTER) THEN
                  SUM_Middle = 0
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 1)-SNF_BOR_ARRAY(jj, 1, ii, 1))+SNF_BOR_ARRAY(jj, 1, ii, 1)
                  temp2 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)-SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  temp3 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 2)-SNF_BOR_ARRAY(jj, 1, ii, 2))+SNF_BOR_ARRAY(jj, 1, ii, 2)
                  temp4 = mid*(SNF_BOR_ARRAY(jj, 2, ii, 3)-SNF_BOR_ARRAY(jj, 1, ii, 3))+SNF_BOR_ARRAY(jj, 1, ii, 3)
                  B0    = temp1
                  B1    = (temp2 - temp1)/(bor2-bor1)
                  B2    = (temp3 - B0 - B1 * (bor3-bor1))/((bor3-bor1)*(bor3-bor2))
                  B3    = (temp4 - B0 - B1 * (bor4-bor1) - B2*(bor4-bor1)*(bor4-bor2))/((bor4-bor1)*(bor4-bor2)*(bor4-bor3))
                  SUM_Middle = 0d0
                  SUM_Middle = B0 + B1*(tmp_ref-bor1) + B2*(tmp_ref-bor1)*(tmp_ref-bor2)
                  SUM_Middle = SUM_Middle + B3*(tmp_ref-bor1)*(tmp_ref-bor2)*(tmp_ref-bor3)
                  CAL_BOR_RESULT(jj,2,ii) = SUM_Middle
               ENDIF
            ENDDO
         ENDDO
         endif
#endif
      ELSE
         DO jj = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               CAL_BOR_RESULT(jj,1,ii) = SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii)
               IF (FLAG_SNF_BU_INTER) THEN
                  mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                            /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                  temp1 = mid*(SNF_REF_ARRAY(jj, SNF_N_BU_UP, ii) - SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(jj, SNF_N_BU_DOWN, ii)
                  CAL_BOR_RESULT(jj,2,ii) = temp1
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
      ENDIF

      ! 03. CALCULATE power correction factor
      CALL SNF_DATA
      !CALL SNF_DATA_CORRECTION
      if (.not.FLAG_NO_PCF)  then
         CALL SNF_CAL_PCF
      else
         IF(allocated(PCF_RESULT)     ) deallocate(PCF_RESULT     )
         IF (Flag_SNF_PIN) THEN
            allocate(PCF_RESULT     (int(npin)**2,SNF_N_ISO))
         ELSE
            allocate(PCF_RESULT     (SNF_MAX_MAT,SNF_N_ISO))
         ENDIF
         PCF_RESULT      = 0d0
      endif
      PCF_96242 = 0d0
      rr=1
      kk = rr
      ! 04. UPDATE ND WITH HISTORY INDEX AND PCF
      ! 04-1) initalization
      IF (FLAG_NO_PCF) PCF_RESULT = 0d0
      If(allocated(SNF_RESULT_ARRAY)) deallocate(SNF_RESULT_ARRAY)
      If(allocated(SNF_RESULT_ARRAY_OLD)) deallocate(SNF_RESULT_ARRAY_OLD)
      IF (.not.Flag_SNF_PIN) THEN
         allocate(SNF_RESULT_ARRAY(SNF_N_INDEXING,SNF_N_ISO))
         allocate(SNF_RESULT_ARRAY_OLD(SNF_N_INDEXING,SNF_N_ISO))
         DO kk = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               IF (FLAG_SNF_BU_INTER) THEN
                  !write(*,*) SNF_ISO_ARRAY(ii), PCF_RESULT(1,ii)
                  IF (PCF_RESULT(kk,ii) .NE. 0) THEN
                     SNF_RESULT_ARRAY(kk,ii) = PCF_RESULT(kk,ii)*(CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))
                  ELSE
                     mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                               /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                     temp1 = mid*(SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii)-SNF_REF_ARRAY(kk, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(kk, SNF_N_BU_DOWN, ii)
                     SNF_RESULT_ARRAY(kk,ii) = (CAL_BOR_RESULT(kk,2,ii)+CAL_TMO_RESULT(kk,2,ii)+CAL_TFU_RESULT(kk,2,ii)- &
                                                2d0*temp1)
                  ENDIF
               ELSE
                  IF (PCF_RESULT(kk,ii) .NE. 0) THEN
                     SNF_RESULT_ARRAY(kk,ii) = PCF_RESULT(kk,ii)*(CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))
                  ELSE
                     SNF_RESULT_ARRAY(kk,ii) = (CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))
                  ENDIF
               ENDIF
               if (SNF_RESULT_ARRAY(kk,ii)<0d0) SNF_RESULT_ARRAY(kk,ii) = 0d0
            ENDDO
         ENDDO
      ELSE
         allocate(SNF_RESULT_ARRAY(int(npin)**2,SNF_N_ISO))
         allocate(SNF_RESULT_ARRAY_OLD(int(npin)**2,SNF_N_ISO))
         IF (allocated(NN_MAT)) deallocate(NN_MAT)
         allocate(NN_MAT(SNF_MAX_MAT))
         NN_MAT = 0
         DO rr = 1, (int(npin)**2)
            !kk = SNF_PIN_MAT(rr,1)
            kk = rr
            NN_MAT(kk) = NN_MAT(kk)+1
         ENDDO
         DO kk = SNF_MIN_MAT,SNF_MAX_MAT
            IF (NN_MAT(kk)==0) NN_MAT(kk)=1
         ENDDO
         !PCF_96242 = 0d0
         DO rr = 1, (int(npin)**2)
            kk = rr
            DO ii = 1, SNF_N_ISO
               IF (FLAG_SNF_BU_INTER) THEN
                  IF (PCF_RESULT(rr,ii) .NE. 0) THEN
                     SNF_RESULT_ARRAY(rr,ii) = PCF_RESULT(rr,ii)*(CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))!/real(NN_MAT(kk))
                  ELSE
                     mid   = (HIST_BURNUP_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN)) &
                                               /(SNF_BRAN_BUHIST(SNF_N_BU_UP)-SNF_BRAN_BUHIST(SNF_N_BU_DOWN))
                     temp1 = mid*(SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii)-SNF_REF_ARRAY(kk, SNF_N_BU_DOWN, ii))+SNF_REF_ARRAY(kk, SNF_N_BU_DOWN, ii)
                     SNF_RESULT_ARRAY(rr,ii) = (CAL_BOR_RESULT(kk,2,ii)+CAL_TMO_RESULT(kk,2,ii)+CAL_TFU_RESULT(kk,2,ii)- &
                                                2d0*temp1)!/real(NN_MAT(kk))
                  ENDIF
               ELSE
                  IF (PCF_RESULT(rr,ii) .NE. 0) THEN
                     SNF_RESULT_ARRAY(rr,ii) = PCF_RESULT(rr,ii)*(CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))!/real(NN_MAT(kk))
                  ELSE
                     SNF_RESULT_ARRAY(rr,ii) = (CAL_BOR_RESULT(kk,1,ii)+CAL_TMO_RESULT(kk,1,ii)+CAL_TFU_RESULT(kk,1,ii)- &
                                                             2d0*SNF_REF_ARRAY(kk, SNF_N_BU_UP, ii))!/real(NN_MAT(kk))
                     !SNF_RESULT_ARRAY(rr,ii) = SNF_RESULT_ARRAY(rr,ii)*15.33d0
                  ENDIF
               ENDIF
               if (SNF_RESULT_ARRAY(rr,ii)<0d0) SNF_RESULT_ARRAY(rr,ii) = 0d0
            ENDDO
         ENDDO
      ENDIF

      !(2*meshsize_x(1))**2*meshsize_z(SNF_K)*Hist_Power_FA (i_step,SNF_Ixy_1N,SNF_K)
      R2_ISO_LIST  = &
      (/53135, 54135, 60147, 60148, 60149, 61147, 61148, 61149, 61548, 62147, 62148, 62149, 64152, 64154, &
        64155, 64156, 64157, 64158, 64160, 92234, 92235, 92236, 92237, 92238, 93237, 93238, 93239, 94238, &
        94239, 94240, 94241, 94242, 94243, 95241, 95242, 95243, 95244, 95642, 96242, 96243, 96244 /)

      mid_iso = 1
      DO ii = 1, SNF_N_ISO
         DO WHILE (int(SNF_ISO_ARRAY(ii)) > R2_ISO_LIST(mid_iso))
            mid_iso = mid_iso+1
            if (mid_iso > 41) go to 888
         ENDDO
         IF(int(SNF_ISO_ARRAY(ii)) == R2_ISO_LIST(mid_iso)) THEN
            R2_ISO_MATCH(mid_iso) = ii
         ENDIF
      ENDDO

888   CONTINUE

      IF (.not.Flag_SNF_PIN) THEN
       DO kk = SNF_MIN_MAT,SNF_MAX_MAT
        mid_iso = 1
        ! Generally, 2*meshsize_x = Gridsize_x
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_I35_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Xe35_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Ps48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd52_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd54_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd55_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd56_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd57_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd58_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd60_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U34_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U35_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U36_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U37_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U38_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np37_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np38_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np39_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu38_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu39_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu40_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu41_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am41_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am44_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        ! SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_As42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        mid_iso=mid_iso+1  !! JUST FOR TEST
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm44_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24;mid_iso=mid_iso+1
       ENDDO
      !ELSE if(flag_snf_decay) then
      ELSE if(.true.) then
       pin_count = int(N_pin,4) !npin*npin
       if (flag_snf_coo) pin_count = (int(npin)**2)
       DO rr = 1,(int(npin)**2)
        kk = rr
        ip = SNF_PIN_POSITION(rr,1)
        jp = SNF_PIN_POSITION(rr,2)
        mid_iso = 1
        ! For TEST
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_I35_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Xe35_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Nd49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pm49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Ps48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm47_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm48_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Sm49_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd52_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd54_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd55_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd56_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd57_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd58_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Gd60_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U34_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U35_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U36_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U37_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_U38_FA (SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np37_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np38_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Np39_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu38_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu39_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu40_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu41_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Pu43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am41_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Am44_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        mid_iso=mid_iso+1  !! JUST FOR TEST
        !SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_As42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm42_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm43_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(kk,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*Hist_N_Cm44_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)*1d-24/real(pin_count);mid_iso=mid_iso+1

       ENDDO
      ELSE
       pin_count = 0
       pin_count = int(N_pin) !236
       if (flag_snf_decay) pin_count = npin*npin
       !write(*,*) 'pin_count', pin_count
       DO rr = 1,(int(npin)**2)
        kk = rr
        ip = SNF_PIN_POSITION(rr,1)
        jp = SNF_PIN_POSITION(rr,2)
        mid_iso = 1
        !write(*,*) ip,jp,SNF_Hist_N_U35 (SNF_N_END_CYBU,ip,jp)
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_I35 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Xe35(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Nd47(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Nd48(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Nd49(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pm47(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pm48(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pm49(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Ps48(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Sm47(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Sm48(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Sm49(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd52(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd54(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd55(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd56(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd57(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd58(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Gd60(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_U34 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_U35 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_U36 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_U37 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_U38 (SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Np37(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Np38(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Np39(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu38(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu39(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu40(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu41(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu42(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Pu43(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Am41(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Am42(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Am43(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Am44(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_As42(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Cm42(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        !IF ((rr==1)) PCF_96242=1d0/(PCF_96242/SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso)))
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Cm43(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
        SNF_RESULT_ARRAY(rr,R2_ISO_MATCH(mid_iso))=((2*meshsize_x(1))**2)*meshsize_z(SNF_K)*SNF_Hist_N_Cm44(SNF_N_END_CYBU,ip,jp)*1d-24/real(pin_count);mid_iso=mid_iso+1
       ENDDO

      ENDIF

      IF (Flag_SNF_PIN) THEN
         WRITE(W_SNF,'(A99)') '==================================================================================================='
         WRITE(W_SNF,'(A99)') '|+                                     RAST-K SNF SUMMARY                                        +|'
         WRITE(W_SNF,'(A99)') '==================================================================================================='
         WRITE(W_SNF,'(A7)') '|OPTION'
         WRITE(W_SNF,'(A42)')       '1. TYPE               : PIN-WISE[MAT-WISE]'
         WRITE(W_SNF,'(A23,1x,A)') '2. ND FILE NAME       :', NAME_SNF_INTER//'.nd'
         WRITE(W_SNF,'(A99)') '==================================================================================================='
         WRITE(W_SNF,'(A8)') '|PIN MAP'
         WRITE(W_SNF,'(A14,1x,I6,A6,I6,A5)')  '  FA POSITION:', SNF_Ixy_1N,' (Ixy)', SNF_k, ' (Iz)'
         WRITE(W_SNF,'(A9,I3)') '  NPIN : ',int(npin)
         WRITE(W_SNF,'(A9)') '- MAT-MAP'
         DO kk = 1, int(npin)
            DO jj = 1, int(npin)
               IF(((SNF_MIN_MAT-1)<SNF_PIN_MAT((kk-1)*int(npin)+jj,1)) &
                   .and. (SNF_PIN_MAT((kk-1)*int(npin)+jj,1)<(SNF_MAX_MAT+1))) THEN
                  WRITE(W_SNF,'(1x,I4)',advance='no') SNF_PIN_MAT((kk-1)*int(npin)+jj,1)
               ELSE
                  WRITE(W_SNF,'(1x,I4)',advance='no') 0
               ENDIF
            ENDDO
            WRITE(W_SNF,'(A)') ''
         ENDDO
         WRITE(W_SNF,'(A)') ' '
         WRITE(W_SNF,'(A9)') '- PIN-MAP'
         DO kk = 1, int(npin)
            DO jj = 1, int(npin)
               IF(((SNF_MIN_MAT-1)<SNF_PIN_MAT((kk-1)*int(npin)+jj,1)) &
                   .and. (SNF_PIN_MAT((kk-1)*int(npin)+jj,1)<(SNF_MAX_MAT+1))) THEN
                  WRITE(W_SNF,'(1x,I4)',advance='no') int((kk-1)*int(npin)+jj)
               ELSE
                  WRITE(W_SNF,'(1x,I4)',advance='no') 0
               ENDIF
            ENDDO
            WRITE(W_SNF,'(A)') ''
         ENDDO
         WRITE(W_SNF,'(A)') ' '
         WRITE(W_SNF,*) ' '
         WRITE(W_SNF,'(A99)') '==================================================================================================='
         IF(allocated(SNF_PRINT_SUM)) deallocate(SNF_PRINT_SUM)
         Allocate(SNF_PRINT_SUM(int(npin)**2,7))
         DO jj = 1, (int(npin)**2)
            kk = jj
            IF ((SNF_MIN_MAT> kk) .or. (kk>SNF_MAX_MAT)) cycle
            WRITE(W_SNF,*) ' '
#ifdef JR_SRCTRM
            CALL GET_SRCTRM(jj,SNF_N_ISO,SNF_ISO_ARRAY(:),SNF_RESULT_ARRAY(jj,:),1d-30*SNF_VOL_RATIO(kk,SNF_N_BU_UP) &
                                                                  ,.FALSE._1)
            CALL END_SRC
#endif
         ENDDO
      ELSE
         IF (Flag_SNF_First) THEN
            WRITE(W_SNF,'(A99)') '==================================================================================================='
            WRITE(W_SNF,'(A99)') '|+                                     RAST-K SNF SUMMARY                                        +|'
            WRITE(W_SNF,'(A99)') '==================================================================================================='
            WRITE(W_SNF,'(A7)') '|OPTION'
            WRITE(W_SNF,'(A31)')       '1. TYPE               : FA-WISE'
            WRITE(W_SNF,'(A23,1x,A)') '2. ND FILE NAME       :', NAME_SNF_INTER//'.nd'
            IF (.not.Flag_SNF_NFA) THEN
               WRITE(W_SNF,'(A23,1x,I6,A6,I6,A5)')  '3. FA POSITION        :', SNF_Ixy_1N,' (Ixy)', SNF_k, ' (Iz)'
            ELSE
               DO kk = 1, SNF_NFA
                  IF (kk==1) THEN
                     WRITE(W_SNF,'(A17,I5,A3,I3,A1,I6,A15,I4,A1)',advance='no')  '3. FA POSITION(N=',SNF_PRINT_NN,'):<', kk, '>', &
                                        SNF_Ixy_1N_M(kk),' (Ixy)     <NZ=', SNF_NZ(kk), '>'
                     DO jj = 1,SNF_NZ(kk)
                        IF (jj < SNF_NZ(kk)) THEN
                           WRITE(W_SNF,'(I6,A1)',advance='no') SNF_K_M(kk,jj),','
                        ELSE
                           WRITE(W_SNF,'(I6,A5)') SNF_K_M(kk,jj),' (IZ)'
                        ENDIF
                     ENDDO
                  ELSE
                     WRITE(W_SNF,'(A17,A5,A3,I3,A1,I6,A15,I4,A1)',advance='no')  ' ',' ','  <', kk, '>', &
                                        SNF_Ixy_1N_M(kk),' (Ixy)     <NZ=', SNF_NZ(kk), '>'
                     DO jj = 1,SNF_NZ(kk)
                        IF (jj < SNF_NZ(kk)) THEN
                           WRITE(W_SNF,'(I6,A1)',advance='no') SNF_K_M(kk,jj),','
                        ELSE
                           WRITE(W_SNF,'(I6,A5)') SNF_K_M(kk,jj),' (IZ)'
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               WRITE(W_SNF,'(A)') ' '
               WRITE(W_SNF,'(A33)') '                        1) RADIAL'
            ENDIF
            do iy_1n=1,ny_1n
               if (Ix_Start_y_1n(iy_1n)>0) then
                  if (Ix_Start_y_1n(Iy_1n)==1) then
                     write(FMT,*) "(A24, ", Ix_End_y_1n(iy_1n)-Ix_Start_y_1n(iy_1n)+1, " i5 )"
                     write(W_SNF,FMT) ' ', &
                        ixiy_1ntoixy_1n(Ix_Start_y_1n(iy_1n):Ix_End_y_1n(iy_1n), iy_1n)
                  else
                     write(FMT,*) "(A24, ", Ix_Start_y_1n(iy_1n)-1," a5, ",Ix_End_y_1n(iy_1n)-Ix_Start_y_1n(iy_1n)+1," i5 )"
                     write(W_SNF,FMT) ' ',  &
                        Blank_Line(1:Ix_Start_y_1n(iy_1n)-1), &
                        ixiy_1ntoixy_1n( Ix_Start_y_1n(iy_1n):Ix_End_y_1n(iy_1n), iy_1n)
                  endif
               endif
            enddo
            WRITE(W_SNF,'(A)') ' '
            WRITE(W_SNF,'(A91)') '                        2) AXIAL (R:Reflector, 0:SRC calculation off, 1:SRC calculation on)'
            WRITE(W_SNF,'(A28)',advance='no') '                        NODE'
            IF (.not.Flag_SNF_NFA) THEN
               WRITE(W_SNF,'(A5,I3,A1)',advance='no') '    <',1,'>'
               WRITE(W_SNF,*) ''
               Do Iz=1,Nz
                  IF ((Iz < IzFuelBot) .or. (Iz > IzFuelTop)) THEN
                     WRITE(W_SNF,'(24x,I4)',advance='no') Iz
                     DO kk = 1, SNF_NFA
                        WRITE(W_SNF,'(A8)',advance='no') 'R'
                     ENDDO
                     WRITE(W_SNF,*) ''
                  ELSE
                     WRITE(W_SNF,'(24x,I4,A8)',advance='no') Iz
                     temp_int   = 0
                     if (Iz == SNF_K) temp_int=1
                     WRITE(W_SNF,'(I8)',advance='no') temp_int
                     WRITE(W_SNF,*) ''
                  ENDIF
               ENDDO
            ELSE
               DO kk = 1, SNF_NFA
                  WRITE(W_SNF,'(A5,I3,A1)',advance='no') '    <',kk,'>'
               ENDDO
               WRITE(W_SNF,*) ''
               Do Iz=1,Nz
                  IF ((Iz < IzFuelBot) .or. (Iz > IzFuelTop)) THEN
                     WRITE(W_SNF,'(24x,I4)',advance='no') Iz
                     DO kk = 1, SNF_NFA
                        WRITE(W_SNF,'(A8)',advance='no') 'R'
                     ENDDO
                     WRITE(W_SNF,*) ''
                  ELSE
                     WRITE(W_SNF,'(24x,I4,A8)',advance='no') Iz
                     DO kk = 1, SNF_NFA
                        temp_int   = 0
                        DO jj = 1, SNF_NZ(kk)
                           if (Iz == SNF_K_M(kk,jj)) temp_int=1
                        ENDDO
                        WRITE(W_SNF,'(I8)',advance='no') temp_int
                     ENDDO
                     WRITE(W_SNF,*) ''
                  ENDIF
               ENDDO
            ENDIF
            WRITE(W_SNF,'(A)') ' '
            WRITE(W_SNF,'(A99)') '==================================================================================================='
            WRITE(W_SNF,'(A54)') '                                         PRINT INDEX=1'
         ELSE
            WRITE(W_SNF,'(A)') ' '
            WRITE(W_SNF,'(A99)') '==================================================================================================='
            WRITE(W_SNF,'(A53,I1)') '                                         PRINT INDEX=',SNF_PRINT_INDEX
         ENDIF
         WRITE(W_SNF,'(A99)') '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(W_SNF,'(A11)') 'FA POSITION'
         WRITE(W_SNF,'(A4,I6)') 'Ixy:', SNF_Ixy_1N
         WRITE(W_SNF,'(A4,I6)') 'Iz :', SNF_k

#ifdef JR_SRCTRM
         CALL GET_SRCTRM(0,SNF_N_ISO,SNF_ISO_ARRAY(:),SNF_RESULT_ARRAY(1,:),1d-30*SNF_VOL_RATIO(1,SNF_N_BU_UP) &
                                                                  ,.FALSE._1)
#endif
      ENDIF

#ifdef JR_SRCTRM
      CALL END_SRC
#endif
      CALL SNF_PRINT

      RETURN
      END SUBROUTINE SNF_CAL


      SUBROUTINE SNF_DATA_CORRECTION
      USE Inc_History
      IMPLICIT NONE
      INTEGER(4)                              :: ii
      INTEGER(4),DIMENSION(506)               :: R2_ISO_LIST
      INTEGER(4)                              :: mid_iso

      R2_ISO_LIST  = &
      (/   2006,   7016,  30073,  30074,  30075,  30076,  30077,  30078,  30079,  31073,  &
          31074,  31075,  31076,  31077,  31078,  31079,  31080,  31081,  31082,  31083,  &
          31084,  32075,  32077,  32078,  32079,  32080,  32081,  32082,  32083,  32084,  &
          32085,  32086,  32087,  32473,  32479,  32481,  33077,  33078,  33079,  33080,  &
          33081,  33082,  33083,  33084,  33085,  33086,  33087,  33088,  33089,  33482,  &
          34081,  34083,  34084,  34085,  34086,  34087,  34088,  34089,  34090,  34091,  &
          34479,  34481,  34483,  35082,  35083,  35084,  35085,  35086,  35087,  35088,  &
          35089,  35090,  35091,  35092,  35093,  35094,  35482,  35484,  36087,  36088,  &
          36089,  36090,  36091,  36092,  36093,  36094,  36095,  36096,  36097,  36098,  &
          36483,  36485,  37088,  37089,  37090,  37091,  37092,  37093,  37094,  37095,  &
          37096,  37097,  37098,  37099,  37100,  37490,  38091,  38092,  38093,  38094,  &
          38095,  38096,  38097,  38098,  38099,  38100,  38101,  38102,  39092,  39093,  &
          39094,  39095,  39096,  39097,  39098,  39099,  39100,  39101,  39102,  39103,  &
          39104,  39491,  39493,  39496,  39497,  39498,  40097,  40098,  40099,  40100,  &
          40101,  40102,  40103,  40104,  40105,  40106,  41096,  41097,  41098,  41099,  &
          41100,  41101,  41102,  41103,  41104,  41105,  41106,  41107,  41108,  41109,  &
          41497,  41498,  41499,  41500,  41502,  41504,  42099,  42101,  42102,  42103,  &
          42104,  42105,  42106,  42107,  42108,  42109,  42110,  42111,  43100,  43101,  &
          43102,  43103,  43104,  43105,  43106,  43107,  43108,  43109,  43110,  43111,  &
          43112,  43113,  43114,  43499,  43502,  44105,  44107,  44108,  44109,  44110,  &
          44111,  44112,  44113,  44114,  44115,  44116,  45107,  45108,  45109,  45110,  &
          45111,  45112,  45113,  45114,  45115,  45116,  45117,  45118,  45119,  45505,  &
          45506,  45508,  45510,  46109,  46111,  46112,  46113,  46114,  46115,  46116,  &
          46117,  46118,  46119,  46120,  46121,  46122,  46509,  47110,  47111,  47112,  &
          47113,  47114,  47115,  47116,  47117,  47118,  47119,  47120,  47121,  47122,  &
          47123,  47124,  47125,  47509,  47511,  47513,  47515,  47516,  47517,  47518,  &
          47520,  47522,  48115,  48117,  48118,  48119,  48120,  48121,  48122,  48123,  &
          48124,  48125,  48126,  48127,  48128,  48130,  48131,  48517,  48519,  48521,  &
          49116,  49117,  49118,  49119,  49120,  49121,  49122,  49123,  49124,  49125,  &
          49126,  49127,  49128,  49129,  49130,  49131,  49132,  49133,  49515,  49517,  &
          49519,  49520,  49521,  49522,  49523,  49524,  49525,  49526,  49527,  49528,  &
          49529,  49530,  49531,  50121,  50125,  50127,  50128,  50129,  50130,  50131,  &
          50132,  50133,  50134,  50135,  50523,  50525,  50527,  50528,  50529,  50530,  &
          50531,  51122,  51126,  51127,  51128,  51129,  51130,  51131,  51132,  51133,  &
          51134,  51135,  51136,  51137,  51138,  51526,  51528,  51530,  51532,  51534,  &
          52127,  52129,  52131,  52132,  52133,  52134,  52135,  52136,  52137,  52138,  &
          52139,  52140,  52531,  52533,  53130,  53131,  53132,  53133,  53134,  53136,  &
          53137,  53138,  53139,  53140,  53141,  53142,  53143,  53530,  53532,  53533,  &
          53534,  53536,  54133,  54137,  54138,  54139,  54140,  54141,  54142,  54143,  &
          54144,  54145,  54531,  54533,  54534,  54535,  55138,  55139,  55140,  55141,  &
          55142,  55143,  55144,  55145,  55146,  55147,  55535,  55536,  55538,  56139,  &
          56140,  56141,  56142,  56143,  56144,  56145,  56146,  56147,  56148,  56149,  &
          56150,  57140,  57141,  57142,  57143,  57144,  57145,  57146,  57147,  57148,  &
          57149,  57150,  57151,  57152,  57546,  58141,  58143,  58145,  58146,  58147,  &
          58148,  58149,  58150,  58151,  58152,  58153,  58154,  59143,  59145,  59146,  &
          59147,  59148,  59149,  59150,  59151,  59152,  59153,  59154,  59155,  59156,  &
          59544,  59548,  60151,  60152,  60153,  60154,  60155,  60156,  60157,  60158,  &
          61151,  61152,  61153,  61154,  61155,  61156,  61157,  61158,  61159,  61552,  &
          61554,  62153,  62155,  62156,  62157,  62158,  62159,  62160,  62161,  63157,  &
          63158,  63159,  63160,  63161,  63162,  63554,  64159,  64161,  64162,  64163,  &
          65161,  65162,  65163,  66165,  92239,  95644 /)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_DATA_CORRECTION] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      mid_iso = 1
      DO ii = 1, SNF_N_ISO
         DO WHILE (int(SNF_ISO_ARRAY(ii)) > R2_ISO_LIST(mid_iso))
            mid_iso = mid_iso+1
            if (mid_iso > 506) go to 888
         ENDDO
         IF(int(SNF_ISO_ARRAY(ii)) == R2_ISO_LIST(mid_iso)) THEN
            SNF_PCF_LIST(ii,1) = 1
         ENDIF
      ENDDO

888   continue

      RETURN
      END SUBROUTINE SNF_DATA_CORRECTION


      SUBROUTINE SNF_CAL_PCF
      USE Inc_History
      USE Inc_Geometry,     only: I_1Nto4N
      USE Inc_Pinvar,       only: npin
      IMPLICIT NONE
      INTEGER(4)                              :: ii
      INTEGER(4)                              :: jj
      INTEGER(4)                              :: kk
      INTEGER(4)                              :: rr
      INTEGER(4)                              :: jr_ii
      INTEGER(4)                              :: i
      INTEGER(4)                              :: m
      INTEGER(4)                              :: match_num
      REAL(8)                                 :: temp1, temp2
      ! Assembly-wise flux
      REAL(8),DIMENSION(:,:),ALLOCATABLE        :: snf_r2_flux
      INTEGER(4)                                :: mid_num
      ! PCF parameters
      REAL(8)                                 :: r_POW
      REAL(8)                                 :: LAMBDA
      REAL(8)                                 :: LAMBDA_1, LAMBDA_2
      REAL(8)                                 :: dT
      INTEGER(4)                              :: PCF_N_BU
      INTEGER(4)                              :: bu_r2
      INTEGER(4)                              :: bu_bran
      INTEGER(4)                              :: bu
      ! results
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE      :: FACTOR_PRE_R2
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE      :: FACTOR_PRE_BRAN
      REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE    :: FACTOR_DAU_R2
      REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE    :: FACTOR_DAU_BRAN
      REAL(8),DIMENSION(:),ALLOCATABLE    :: tmp_HRI_DAY
      ! daughter pcf
      LOGICAL(1)                              :: flag_pre
      REAL(8)                                 :: DELTA_C_R2
      REAL(8)                                 :: DELTA_C_BRAN
      REAL(8)                                 :: DECAY_C_R2
      REAL(8)                                 :: DECAY_C_BRAN
      REAL(8)                                 :: BEFORE_C_1
      REAL(8)                                 :: BEFORE_C_2
      INTEGER(4)                              :: SNF_MIN_MAT_old
      INTEGER(4)                              :: SNF_MAX_MAT_old
      REAL(8)                                 :: mid_min_cur
      REAL(8)                                 :: mid_min_val


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_CAL_PCF] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      SNF_MIN_MAT_old = SNF_MIN_MAT
      SNF_MAX_MAT_old = SNF_MAX_MAT
      !SNF_MIN_MAT = 1
      !SNF_MAX_MAT = 1

      ! BRANCH STEPS: SNF_N_BU_UP
      ! R2 STEPS    : SNF_N_END_CYBU
      ! relogical burnup step
      PCF_N_BU = 0
      if (flag_hri) then
      !if (.not.flag_snf_pin .and. flag_hri) then
         IF (SNF_N_BU_UP == (SNF_N_END_CYBU+SNF_HRI_N)) THEN
            PCF_N_BU = SNF_N_BU_UP
         ELSEIF (SNF_N_BU_UP>(SNF_N_END_CYBU+SNF_HRI_N)) THEN
            PCF_N_BU = SNF_N_END_CYBU+SNF_HRI_N
         ELSEIF (SNF_N_BU_UP<(SNF_N_END_CYBU+SNF_HRI_N)) THEN
            PCF_N_BU = SNF_N_BU_UP
         ENDIF
      else
         IF (SNF_N_BU_UP == SNF_N_END_CYBU) THEN
            PCF_N_BU = SNF_N_BU_UP
         ELSEIF (SNF_N_BU_UP>SNF_N_END_CYBU) THEN
            PCF_N_BU = SNF_N_END_CYBU
         ELSEIF (SNF_N_BU_UP<SNF_N_END_CYBU) THEN
            PCF_N_BU = SNF_N_BU_UP
         ENDIF
      endif

      ! initialization
      IF(allocated(FACTOR_PRE_R2)  ) deallocate(FACTOR_PRE_R2  )
      IF(allocated(FACTOR_DAU_R2)  ) deallocate(FACTOR_DAU_R2  )! 1: generated by pre 1, 2: generated by pre 2
      IF(allocated(FACTOR_PRE_BRAN)) deallocate(FACTOR_PRE_BRAN)
      IF(allocated(FACTOR_DAU_BRAN)) deallocate(FACTOR_DAU_BRAN)! 1: generated by pre 1, 2: generated by pre 2
      IF(allocated(snf_r2_flux)    ) deallocate(snf_r2_flux    )
      IF(allocated(PCF_RESULT)     ) deallocate(PCF_RESULT     )
      IF(allocated(tmp_HRI_DAY)    ) deallocate(tmp_HRI_DAY    )
      !IF (.FALSE.) THEN
      IF (Flag_SNF_PIN) THEN
         allocate(FACTOR_PRE_R2  (int(npin)**2,PCF_N_BU,SNF_N_ISO))
         allocate(FACTOR_DAU_R2  (int(npin)**2,2,PCF_N_BU,SNF_N_ISO)) ! 1: generated by pre 1, 2: generated by pre 2
         allocate(FACTOR_PRE_BRAN(int(npin)**2,PCF_N_BU,SNF_N_ISO))
         allocate(FACTOR_DAU_BRAN(int(npin)**2,2,PCF_N_BU,SNF_N_ISO)) ! 1: generated by pre 1, 2: generated by pre 2
         allocate(snf_r2_flux    (int(npin)**2,PCF_N_BU))
         allocate(PCF_RESULT     (int(npin)**2,SNF_N_ISO))
         !if (flag_hri) allocate(tmp_HRI_DAY(PCF_N_BU))
         allocate(tmp_HRI_DAY(PCF_N_BU))
      ELSE
         allocate(FACTOR_PRE_R2  (SNF_MAX_MAT,PCF_N_BU,SNF_N_ISO))
         allocate(FACTOR_DAU_R2  (SNF_MAX_MAT,2,PCF_N_BU,SNF_N_ISO)) ! 1: generated by pre 1, 2: generated by pre 2
         allocate(FACTOR_PRE_BRAN(SNF_MAX_MAT,PCF_N_BU,SNF_N_ISO))
         allocate(FACTOR_DAU_BRAN(SNF_MAX_MAT,2,PCF_N_BU,SNF_N_ISO)) ! 1: generated by pre 1, 2: generated by pre 2
         allocate(snf_r2_flux    (SNF_MAX_MAT,PCF_N_BU))
         allocate(PCF_RESULT     (SNF_MAX_MAT,SNF_N_ISO))
         !if (flag_hri) allocate(tmp_HRI_DAY(PCF_N_BU))
         allocate(tmp_HRI_DAY(PCF_N_BU))
      ENDIF
      FACTOR_PRE_R2   = 0d0
      FACTOR_PRE_BRAN = 0d0
      FACTOR_DAU_R2   = 0d0
      FACTOR_DAU_BRAN = 0d0
      snf_r2_flux     = 0d0
      PCF_RESULT      = 0d0
      if (flag_hri) tmp_HRI_DAY = 0d0
      ! r2 flux
      IF (.not.FLAG_SNF_PIN) THEN
         IF (flag_hri .and. (PCF_N_BU>SNF_N_END_CYBU)) then
            DO i = 1, PCF_N_BU-SNF_N_END_CYBU !SNF_N_END_CYBU+SNF_HRI_N
               bu = SNF_HRI_N-(PCF_N_BU-SNF_N_END_CYBU)+i
               snf_r2_flux(1,i)=SNF_HRI_FLUX(bu)
               tmp_hri_day(i) = snf_hri_day(bu)
            ENDDO
            DO i = 1,SNF_N_END_CYBU
                  bu = i
                  temp1 = 0d0
                  temp2 = 0d0
                  DO m = 1, 4
                     temp1 = temp1 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
                     temp2 = temp2 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
                  ENDDO
                  snf_r2_flux(1,PCF_N_BU-SNF_N_END_CYBU+i) = temp1/4d0 + temp2/4d0
                  tmp_hri_day(PCF_N_BU-SNF_N_END_CYBU+i) = snf_hri_last_day+hist_day(i)
            ENDDO
            correct_flux = snf_r2_flux(1,i)
         else
            DO i = 1, PCF_N_BU !SNF_N_END_CYBU
               bu = SNF_N_END_CYBU-PCF_N_BU+i
               temp1 = 0d0
               temp2 = 0d0
               DO m = 1, 4
                  temp1 = temp1 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
                  temp2 = temp2 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
               ENDDO
               snf_r2_flux(1,i) = temp1/4d0 + temp2/4d0
               !if (flag_hri) tmp_hri_day(i) = hist_day(bu)
               tmp_hri_day(i) = hist_day(bu)
            ENDDO
            correct_flux = snf_r2_flux(1,i)
         endif
      ELSE
         IF (flag_hri .and. (PCF_N_BU>SNF_N_END_CYBU)) then
            DO i = 1, PCF_N_BU-SNF_N_END_CYBU !SNF_N_END_CYBU+SNF_HRI_N
               bu = SNF_HRI_N-(PCF_N_BU-SNF_N_END_CYBU)+i
               snf_r2_flux(:,i)=SNF_HRI_FLUX(bu)
               tmp_hri_day(i) = snf_hri_day(bu)
            ENDDO
            DO i = 1,SNF_N_END_CYBU
                  bu = i
                  temp1 = 0d0
                  temp2 = 0d0
                  DO m = 1, 4
                     temp1 = temp1 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
                     temp2 = temp2 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
                  ENDDO
                  snf_r2_flux(:,PCF_N_BU-SNF_N_END_CYBU+i) = temp1/4d0 + temp2/4d0
                  tmp_hri_day(PCF_N_BU-SNF_N_END_CYBU+i) = snf_hri_last_day+hist_day(i)
            ENDDO
         ELSE
            mid_num=SNF_MIN_MAT
            mid_min_val = 10d+55
            mid_min_cur = 0d0
            DO jj = 1, (int(npin)**2)
               if (phihom_SNF(1,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),1) == 0) cycle
               DO i = 2, PCF_N_BU
                  bu = SNF_N_END_CYBU-PCF_N_BU+i
                  snf_r2_flux(jj,i)=phihom_SNF(bu,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),1)+ &
                                         phihom_SNF(bu,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),2)
               ENDDO
               correct_flux = phihom_SNF(SNF_N_END_CYBU,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),1)+ &
                              phihom_SNF(SNF_N_END_CYBU,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),2)
               mid_min_cur = phihom_SNF(bu,SNF_PIN_POSITION(jj,1),SNF_PIN_POSITION(jj,2),1)
               if ((mid_min_cur > 0d0) .and. (mid_min_val>mid_min_cur)) then
                   mid_min_val = mid_min_cur
                   SNF_MIN_VAL = real(jj)
               endif
            ENDDO
            DO i = 1, PCF_N_BU !SNF_N_END_CYBU
               bu = SNF_N_END_CYBU-PCF_N_BU+i
               temp1 = 0d0
               temp2 = 0d0
               DO m = 1, 4
                  temp1 = temp1 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
                  temp2 = temp2 + Hist_Flux (bu,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
               ENDDO
               snf_r2_flux(:,i) = temp1/4d0 + temp2/4d0
               !if (flag_hri) tmp_hri_day(i) = hist_day(bu)
               tmp_hri_day(i) = hist_day(bu)
            ENDDO
         ENDIF
      ENDIF

      !IF (.TRUE.) THEN
      IF (.not.Flag_SNF_PIN) THEN
         ! 01. PRECURSOR
         DO kk = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               bu_r2   = SNF_N_END_CYBU-PCF_N_BU+2
               if (flag_hri) bu_r2 = 2
               bu_bran = SNF_N_BU_UP-PCF_N_BU+2
               ! R2
               r_POW  = SNF_BRAN_DAYHIST(bu_bran) / tmp_hri_day(bu_r2) ! SAME TIME STEP
               LAMBDA = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
               dT     = tmp_hri_day(bu_r2)*60d0*60d0*24d0
               FACTOR_PRE_R2(kk,2,ii) = r_POW*(1-EXP(-LAMBDA*dT))
               ! BRANCH
               dT     = SNF_BRAN_DAYHIST(bu_bran)*60d0*60d0*24d0
               LAMBDA = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(kk,bu_bran,ii)
               FACTOR_PRE_BRAN(kk,2,ii) = (1-EXP(-LAMBDA*dT))
               DO i = 3, PCF_N_BU
                  bu_r2   = SNF_N_END_CYBU-PCF_N_BU+i
                  if (flag_hri) bu_r2 = i
                  bu_bran = SNF_N_BU_UP-PCF_N_BU+i
                  ! R2
                  r_POW   = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1)) / &
                                                (tmp_hri_day(bu_r2)-tmp_hri_day(bu_r2-1))
                  LAMBDA  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                  dT      = (tmp_hri_day(bu_r2)-tmp_hri_day(bu_r2-1))*60d0*60d0*24d0
                  FACTOR_PRE_R2(kk,i,ii) = FACTOR_PRE_R2(kk,i-1,ii)*EXP(-LAMBDA*dT) + &
                                        r_POW*(1-EXP(-LAMBDA*dT))
                  ! BRANCH
                  dT      = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1))*60d0*60d0*24d0
                  LAMBDA  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(kk,bu_bran,ii)
                  FACTOR_PRE_BRAN(kk,i,ii) = FACTOR_PRE_BRAN(kk,i-1,ii)*EXP(-LAMBDA*dT) + &
                                          (1-EXP(-LAMBDA*dT))
               ENDDO
               IF ((FACTOR_PRE_BRAN(kk,PCF_N_BU,ii) .NE. 0 ) .and. (SNF_PCF_LIST(ii,1)==1))THEN
                  PCF_RESULT(kk,ii) = FACTOR_PRE_R2(kk,PCF_N_BU,ii)/FACTOR_PRE_BRAN(kk,PCF_N_BU,ii)
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
         ! 02 . DAUGHTER
         DO kk = SNF_MIN_MAT,SNF_MAX_MAT
            DO ii = 1, SNF_N_ISO
               IF (SNF_PCF_LIST(ii,1) > 1) THEN
                  flag_pre  = .false.
                  match_num = 0
                  DO jj = 2, SNF_PCF_LIST(ii,1)
                     DO jr_ii = 1, SNF_N_ISO
                        IF (int(SNF_ISO_ARRAY(ii)) == SNF_PCF_LIST(ii,jj)) THEN
                           flag_pre  = .true.
                           match_num = jr_ii
                        ENDIF
                     ENDDO
                     IF (flag_pre) THEN
                        ! ++ INDEX --------------- !
                        ! 1 : PRECURSOR            !
                        ! 2 : DAUGHTER             !
                        ! ------------------------ !
                        DELTA_C_R2     = 0d0
                        DELTA_C_BRAN   = 0d0
                        bu_r2   = SNF_N_END_CYBU-PCF_N_BU+2
                        if (flag_hri) bu_r2 = 2
                        bu_bran = SNF_N_BU_UP-PCF_N_BU+2
                        ! R2
                        LAMBDA_1 = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*snf_r2_flux(kk,bu_r2)
                        LAMBDA_2 = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                        r_POW  = SNF_BRAN_DAYHIST(bu_bran) / tmp_hri_day(bu_r2)
                        dT     = tmp_hri_day(bu_r2)*60d0*60d0*24d0
                        IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                           DELTA_C_R2 = r_POW*LAMBDA_1/LAMBDA_2*(1+                      &
                                         LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)- &
                                         LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                        ENDIF
                        ! BRANCH
                        dT       = SNF_BRAN_DAYHIST(bu_bran)*60d0*60d0*24d0
                        LAMBDA_1 = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*SNF_BRAN_FLUX(kk,bu_bran,match_num)
                        LAMBDA_2 = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(kk,bu_bran,ii)
                        IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                           DELTA_C_BRAN = LAMBDA_1/LAMBDA_2*(1+                          &
                                         LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)- &
                                         LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                        ENDIF
                        ! UPDATE PCF_RESULTS
                        FACTOR_DAU_R2(kk, jj-1, 2, ii )   = DELTA_C_R2
                        FACTOR_DAU_BRAN(kk, jj-1, 2, ii ) = DELTA_C_BRAN
                        DO i = 3, PCF_N_BU
                           bu_r2   = SNF_N_END_CYBU-PCF_N_BU+i
                           if (flag_hri) bu_r2 = i
                           bu_bran = SNF_N_BU_UP-PCF_N_BU+i
                           ! initialization
                           DELTA_C_R2 = 0d0
                           DECAY_C_R2 = 0d0
                           DELTA_C_BRAN = 0d0
                           DECAY_C_BRAN = 0d0
                           ! R2
                           r_POW   = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1)) / &
                                                         (tmp_hri_day(bu_r2)-tmp_hri_day(bu_r2-1))
                           LAMBDA_1  = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*snf_r2_flux(kk,bu_r2)
                           LAMBDA_2  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                           IF (LAMBDA_1 .NE. LAMBDA_2 ) THEN
                              dT        = (tmp_hri_day(bu_r2)-tmp_hri_day(bu_r2-1))*60d0*60d0*24d0
                              DELTA_C_R2 = r_POW*LAMBDA_1/LAMBDA_2*(1+                            &
                                                 LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)-  &
                                                 LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                              BEFORE_C_1 = FACTOR_PRE_R2(kk,i-1,match_num)
                              DECAY_C_R2 = LAMBDA_1/(LAMBDA_2-LAMBDA_1)*BEFORE_C_1*(EXP(-LAMBDA_1*dT)-EXP(-LAMBDA_2*dT))
                              BEFORE_C_2 = FACTOR_DAU_R2(kk,jj-1,i-1,ii)
                              FACTOR_DAU_R2(kk,jj-1,i,ii) = BEFORE_C_2*EXP(-LAMBDA_2*dT)+DECAY_C_R2+DELTA_C_R2
                           ENDIF

                           ! BRANCH
                           LAMBDA_1  = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*SNF_BRAN_FLUX(kk,bu_bran,match_num)
                           LAMBDA_2  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(kk,bu_bran,ii)
                           IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                              dT      = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1))*60d0*60d0*24d0
                              DELTA_C_BRAN = LAMBDA_1/LAMBDA_2*(1+                                  &
                                                 LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)-  &
                                                 LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                              BEFORE_C_1 = FACTOR_PRE_BRAN(kk,i-1,match_num)
                              DECAY_C_BRAN = LAMBDA_1/(LAMBDA_2-LAMBDA_1)*BEFORE_C_1*(EXP(-LAMBDA_1*dT)-EXP(-LAMBDA_2*dT))
                              BEFORE_C_2 = FACTOR_DAU_BRAN(kk,jj-1,i-1,ii)
                              FACTOR_DAU_BRAN(kk,jj-1,i,ii) = BEFORE_C_2*EXP(-LAMBDA_2*dT)+DECAY_C_BRAN+DELTA_C_BRAN
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO ! SNF_PCF_LIST(ii,1)
                  IF (SNF_PCF_LIST(ii,1)==2) THEN
                     IF (FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) .NE. 0) THEN
                        PCF_RESULT(kk,ii) = FACTOR_DAU_R2(kk,1,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii)
                     ENDIF
                  ELSE
                     IF ((FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) .NE. 0) .and. (FACTOR_DAU_BRAN(kk,2,PCF_N_BU,ii) .NE. 0)) THEN
                        PCF_RESULT(kk,ii) = (FACTOR_DAU_R2(kk,1,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) + &
                                          FACTOR_DAU_R2(kk,2,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,2,PCF_N_BU,ii))/2d0
                     ENDIF
                  ENDIF
               ENDIF ! SNF_PCF_LIST(ii,1) > 2
            ENDDO ! ISO
         ENDDO ! MAT
      ELSE
         ! 01. PRECURSOR
         DO kk = 1, (int(npin)**2)
            !rr = SNF_PIN_MAT(kk,1)
            rr = kk
            DO ii = 1, SNF_N_ISO
               bu_r2   = SNF_N_END_CYBU-PCF_N_BU+2
               if (flag_hri) bu_r2 = 2
               bu_bran = SNF_N_BU_UP-PCF_N_BU+2
               ! R2
               r_POW  = SNF_BRAN_DAYHIST(bu_bran) / tmp_Hri_DAY(bu_r2) ! SAME TIME STEP
               LAMBDA = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
               dT     = tmp_Hri_DAY(bu_r2)*60d0*60d0*24d0
               FACTOR_PRE_R2(kk,2,ii) = r_POW*(1-EXP(-LAMBDA*dT))
               ! BRANCH
               dT     = SNF_BRAN_DAYHIST(bu_bran)*60d0*60d0*24d0
               LAMBDA = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(rr,bu_bran,ii)
               FACTOR_PRE_BRAN(kk,2,ii) = (1-EXP(-LAMBDA*dT))
               DO i = 3, PCF_N_BU
                  bu_r2   = SNF_N_END_CYBU-PCF_N_BU+i
                  if (flag_hri) bu_r2 = i
                  bu_bran = SNF_N_BU_UP-PCF_N_BU+i
                  ! R2
                  r_POW   = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1)) / &
                                                (tmp_hri_DAY(bu_r2)-tmp_hri_DAY(bu_r2-1))
                  LAMBDA  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                  dT      = (tmp_hri_DAY(bu_r2)-tmp_hri_DAY(bu_r2-1))*60d0*60d0*24d0
                  FACTOR_PRE_R2(kk,i,ii) = FACTOR_PRE_R2(kk,i-1,ii)*EXP(-LAMBDA*dT) + &
                                        r_POW*(1-EXP(-LAMBDA*dT))
                  ! BRANCH
                  dT      = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1))*60d0*60d0*24d0
                  LAMBDA  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(rr,bu_bran,ii)
                  FACTOR_PRE_BRAN(kk,i,ii) = FACTOR_PRE_BRAN(kk,i-1,ii)*EXP(-LAMBDA*dT) + &
                                          (1-EXP(-LAMBDA*dT))
               ENDDO
               IF ((FACTOR_PRE_BRAN(kk,PCF_N_BU,ii) .NE. 0 ) .and. (SNF_PCF_LIST(ii,1)==1))THEN
                  PCF_RESULT(kk,ii) = FACTOR_PRE_R2(kk,PCF_N_BU,ii)/FACTOR_PRE_BRAN(kk,PCF_N_BU,ii)
               ENDIF
            ENDDO ! ISO
         ENDDO ! MAT
         ! 02 . DAUGHTER
         DO kk = 1, (int(npin)**2)
            !rr = SNF_PIN_MAT(kk,1)
            rr = kk
            DO ii = 1, SNF_N_ISO
               IF (SNF_PCF_LIST(ii,1) > 1) THEN
                  flag_pre  = .false.
                  match_num = 0
                  DO jj = 2, SNF_PCF_LIST(ii,1)
                     DO jr_ii = 1, SNF_N_ISO
                        IF (int(SNF_ISO_ARRAY(ii)) == SNF_PCF_LIST(ii,jj)) THEN
                           flag_pre  = .true.
                           match_num = jr_ii
                        ENDIF
                     ENDDO
                     IF (flag_pre) THEN
                        ! ++ INDEX --------------- !
                        ! 1 : PRECURSOR            !
                        ! 2 : DAUGHTER             !
                        ! ------------------------ !
                        DELTA_C_R2     = 0d0
                        DELTA_C_BRAN   = 0d0
                        bu_r2   = SNF_N_END_CYBU-PCF_N_BU+2
                        if (flag_hri) bu_r2 = 2
                        bu_bran = SNF_N_BU_UP-PCF_N_BU+2
                        ! R2
                        LAMBDA_1 = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*snf_r2_flux(kk,bu_r2)
                        LAMBDA_2 = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                        r_POW  = SNF_BRAN_DAYHIST(bu_bran) / tmp_hri_DAY(bu_r2)
                        dT     = tmp_hri_DAY(bu_r2)*60d0*60d0*24d0
                        IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                           DELTA_C_R2 = r_POW*LAMBDA_1/LAMBDA_2*(1+                      &
                                         LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)- &
                                         LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                        ENDIF
                        ! BRANCH
                        dT       = SNF_BRAN_DAYHIST(bu_bran)*60d0*60d0*24d0
                        LAMBDA_1 = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*SNF_BRAN_FLUX(rr,bu_bran,match_num)
                        LAMBDA_2 = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(rr,bu_bran,ii)
                        IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                           DELTA_C_BRAN = LAMBDA_1/LAMBDA_2*(1+                          &
                                         LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)- &
                                         LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                        ENDIF
                        ! UPDATE PCF_RESULTS
                        FACTOR_DAU_R2(kk, jj-1, 2, ii )   = DELTA_C_R2
                        FACTOR_DAU_BRAN(kk, jj-1, 2, ii ) = DELTA_C_BRAN
                        DO i = 3, PCF_N_BU
                           bu_r2   = SNF_N_END_CYBU-PCF_N_BU+i
                           if (flag_hri) bu_r2 = i
                           bu_bran = SNF_N_BU_UP-PCF_N_BU+i
                           ! initialization
                           DELTA_C_R2 = 0d0
                           DECAY_C_R2 = 0d0
                           DELTA_C_BRAN = 0d0
                           DECAY_C_BRAN = 0d0
                           ! R2
                           r_POW   = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1)) / &
                                                         (tmp_hri_DAY(bu_r2)-tmp_hri_DAY(bu_r2-1))
                           LAMBDA_1  = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*snf_r2_flux(kk,bu_r2)
                           LAMBDA_2  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*snf_r2_flux(kk,bu_r2)
                           IF (LAMBDA_1 .NE. LAMBDA_2 ) THEN
                              dT        = (tmp_hri_DAY(bu_r2)-tmp_hri_DAY(bu_r2-1))*60d0*60d0*24d0
                              DELTA_C_R2 = r_POW*LAMBDA_1/LAMBDA_2*(1+                            &
                                                 LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)-  &
                                                 LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                              BEFORE_C_1 = FACTOR_PRE_R2(kk,i-1,match_num)
                              DECAY_C_R2 = LAMBDA_1/(LAMBDA_2-LAMBDA_1)*BEFORE_C_1*(EXP(-LAMBDA_1*dT)-EXP(-LAMBDA_2*dT))
                              BEFORE_C_2 = FACTOR_DAU_R2(kk,jj-1,i-1,ii)
                              FACTOR_DAU_R2(kk,jj-1,i,ii) = BEFORE_C_2*EXP(-LAMBDA_2*dT)+DECAY_C_R2+DELTA_C_R2
                           ENDIF

                           ! BRANCH
                           LAMBDA_1  = SNF_LAMBDA_ARRAY(match_num) + SNF_XS(match_num)*SNF_BRAN_FLUX(rr,bu_bran,match_num)
                           LAMBDA_2  = SNF_LAMBDA_ARRAY(ii) + SNF_XS(ii)*SNF_BRAN_FLUX(rr,bu_bran,ii)
                           IF (LAMBDA_1 .NE. LAMBDA_2) THEN
                              dT      = (SNF_BRAN_DAYHIST(bu_bran)-SNF_BRAN_DAYHIST(bu_bran-1))*60d0*60d0*24d0
                              DELTA_C_BRAN = LAMBDA_1/LAMBDA_2*(1+                                  &
                                                 LAMBDA_1/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_2*dT)-  &
                                                 LAMBDA_2/(LAMBDA_2-LAMBDA_1)*EXP(-LAMBDA_1*dT))
                              BEFORE_C_1 = FACTOR_PRE_BRAN(kk,i-1,match_num)
                              DECAY_C_BRAN = LAMBDA_1/(LAMBDA_2-LAMBDA_1)*BEFORE_C_1*(EXP(-LAMBDA_1*dT)-EXP(-LAMBDA_2*dT))
                              BEFORE_C_2 = FACTOR_DAU_BRAN(kk,jj-1,i-1,ii)
                              FACTOR_DAU_BRAN(kk,jj-1,i,ii) = BEFORE_C_2*EXP(-LAMBDA_2*dT)+DECAY_C_BRAN+DELTA_C_BRAN
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO ! SNF_PCF_LIST(ii,1)
                  IF (SNF_PCF_LIST(ii,1)==2) THEN
                     IF (FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) .NE. 0) THEN
                        PCF_RESULT(kk,ii) = FACTOR_DAU_R2(kk,1,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii)
                     ENDIF
                  ELSE
                     IF ((FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) .NE. 0) .and. (FACTOR_DAU_BRAN(kk,2,PCF_N_BU,ii) .NE. 0)) THEN
                        PCF_RESULT(kk,ii) = (FACTOR_DAU_R2(kk,1,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,1,PCF_N_BU,ii) + &
                                          FACTOR_DAU_R2(kk,2,PCF_N_BU,ii)/FACTOR_DAU_BRAN(kk,2,PCF_N_BU,ii))/2d0
                     ENDIF
                  ENDIF
               ENDIF ! SNF_PCF_LIST(ii,1) > 2
            ENDDO ! ISO
         ENDDO ! MAT
      ENDIF

      SNF_MIN_MAT = SNF_MIN_MAT_old
      SNF_MAX_MAT = SNF_MAX_MAT_old

      RETURN
      END SUBROUTINE SNF_CAL_PCF


      SUBROUTINE SNF_PRINT
      USE Inc_Geometry,     only: meshsize_z, meshsize_x
      USE Inc_Pinvar,       only: npin
      IMPLICIT NONE
      INTEGER(4)                              :: ii
      INTEGER(4)                              :: kk
      REAL(8)                                 :: temp_vol


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_PRINT] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (Flag_snf_pin) THEN
         WRITE(W_INTER,'(A3, 2x, I2, 2x, I7, 2x, I7)') 'SNF', 0, (SNF_MAX_MAT-SNF_MIN_MAT+1), SNF_N_ISO
         temp_vol = ((2*meshsize_x(1))**2)*meshsize_z(SNF_K)/npin/npin
         DO kk = 1,int(npin)**2
            DO ii = 1, SNF_N_ISO
               WRITE(W_INTER,'(I7, 2x, es13.5, 2x, I7, 2x, es13.5, 2x, I7)') SNF_ISO_ARRAY(ii), SNF_RESULT_ARRAY(kk,ii), kk, temp_vol, ii
            ENDDO
         ENDDO
      ELSE
         WRITE(W_INTER,'(A3, 2x, I2, 2x, es13.5, 2x, I7)') 'SNF', 0, ((2*meshsize_x(1))**2)*meshsize_z(SNF_K), SNF_N_ISO
         DO ii = 1, SNF_N_ISO
            WRITE(W_INTER,'(I7, 2x, es13.5, 2x, 7x, 2x, I7)') SNF_ISO_ARRAY(ii), SNF_RESULT_ARRAY(1,ii), ii
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SNF_PRINT


      SUBROUTINE SNF_HRI
      USE Inc_History
      USE MOD_CHAREDIT,     ONLY: CHARTOUPPER
      IMPLICIT NONE
      INTEGER(4)                           :: ii
      REAL(8)                              :: temp1
      LOGICAL(1)                           :: FLAG_OPEN
      CHARACTER(15)                        :: SubTypeName


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_HRI] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO WHILE (.TRUE.)
         INQUIRE( UNIT = R_HRI, OPENED = Flag_Open )
         IF ( .NOT. Flag_Open ) THEN
            OPEN(R_HRI, FILE = NAME_SNF_HRI)
            EXIT
         ELSE
         R_HRI = R_HRI + 10
         END IF
      END DO
      write(*,*) NAME_SNF_HRI

      SubType: DO WHILE (.TRUE.)
         READ(R_HRI, *, END = 999) SubTypeName
         BACKSPACE(R_HRI)
         CALL CharToUpper(SubTypeName)
         SELECT CASE (SubTypeName)
         CASE ('HRST')
            READ(R_HRI,*,END=999) SubTypeName, SNF_HRI_N, SNF_HRI_LAST_DAY
            IF (Allocated(SNF_HRI_BU))   deallocate(SNF_HRI_BU)
            IF (Allocated(SNF_HRI_DAY))   deallocate(SNF_HRI_DAY)
            IF (Allocated(SNF_HRI_TMO))  deallocate(SNF_HRI_TMO)
            IF (Allocated(SNF_HRI_TFU))  deallocate(SNF_HRI_TFU)
            IF (Allocated(SNF_HRI_BOR))  deallocate(SNF_HRI_BOR)
            IF (Allocated(SNF_HRI_FLUX)) deallocate(SNF_HRI_FLUX)
            Allocate(SNF_HRI_BU  (SNF_HRI_N))
            Allocate(SNF_HRI_DAY (SNF_HRI_N))
            Allocate(SNF_HRI_TMO (SNF_HRI_N))
            Allocate(SNF_HRI_TFU (SNF_HRI_N))
            Allocate(SNF_HRI_BOR (SNF_HRI_N))
            Allocate(SNF_HRI_FLUX(SNF_HRI_N))
            !write(*,*) 'SNF_HRI_N', SNF_HRI_N
            DO ii = 1,SNF_HRI_N
               READ(R_HRI,*,END=999) temp1, SNF_HRI_DAY(ii), SNF_HRI_BU(ii), SNF_HRI_TMO(ii), SNF_HRI_TFU(ii), SNF_HRI_BOR(ii), SNF_HRI_FLUX(ii)
            ENDDO
         END SELECT
      END DO SubType

999   continue

      RETURN
      END SUBROUTINE SNF_HRI


      SUBROUTINE SNF_HRST
      USE Inc_Geometry,     only: I_1Nto4N
      USE Inc_History
      USE Inc_Depletion,    only: N_BU
      USE Inc_Pinvar,       only: npin
      IMPLICIT NONE
      INTEGER(4)                              :: i1
      INTEGER(4)                              :: jj
      INTEGER(4)                              :: m
      REAL(8)                                 :: temp1, temp2
      REAL(8)                                 :: temp_flux


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_HRST] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (Flag_snf_pin) THEN
         WRITE(W_HRST,'(A4,2x,I6)') 'HRST', N_BU
         DO jj = 1, N_BU
            temp1 =0d0
            temp2 =0d0
            DO m = 1,4
               temp1 = temp1 + Hist_Flux (jj,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
               temp2 = temp2 + Hist_Flux (jj,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
            ENDDO
            temp_flux = temp1/4d0+temp2/4d0
            WRITE(W_HRST,'(I5, 2x, es15.3, 2x, es13.5, 2x, es13.5, 2x, es13.5, 2x, es13.5)',advance='no') jj, HIST_DAY(jj), &
                         HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K), &
                         Hist_T_Mod_FA(jj,SNF_Ixy_1N,SNF_K), Hist_T_Fuel_FA(jj,SNF_Ixy_1N,SNF_K), hist_PPM(jj)!,           &
                         !temp_flux
            do i1 =  1, int(npin*npin)
               write(W_HRST,'(2x, es13.5)',advance='no') temp_flux
            enddo
            write(W_HRST,*) ''
         ENDDO
      ELSE
         WRITE(W_HRST,'(A4,2x,I6)') 'HRST', N_BU
         DO jj = 1, N_BU
            temp1 =0d0
            temp2 =0d0
            DO m = 1,4
               temp1 = temp1 + Hist_Flux (jj,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,1)
               temp2 = temp2 + Hist_Flux (jj,I_1Nto4N(SNF_Ixy_1N, m),SNF_k,2)
            ENDDO
            temp_flux = temp1/4d0+temp2/4d0
            WRITE(W_HRST,'(I5, 2x, es15.3, 2x, es13.5, 2x, es13.5, 2x, es13.5, 2x, es13.5, 2x, es13.5)') jj, HIST_DAY(jj), &
                         HIST_BURNUP_FA(jj,SNF_Ixy_1N,SNF_K), &
                         Hist_T_Mod_FA(jj,SNF_Ixy_1N,SNF_K), Hist_T_Fuel_FA(jj,SNF_Ixy_1N,SNF_K), hist_PPM(jj),           &
                         temp_flux
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SNF_HRST


      SUBROUTINE SNF_DATA
      USE Inc_History
      IMPLICIT NONE
      INTEGER(4)                 :: ii
      INTEGER(4)                 :: mid_iso
      INTEGER(4),DIMENSION(1567) :: SNF_ISO_LIST
      INTEGER(4),DIMENSION(1567) :: SNF_ISO_TYPE
      INTEGER(4),DIMENSION(1567) :: SNF_PRE_1
      INTEGER(4),DIMENSION(1567) :: SNF_PRE_2

      SNF_ISO_LIST = &
      (/1001, 1002, 1003, 2003, 2004, 2006, 3006, 3007, 3008, 3009, 4007, 4008, 4009, 4010, &
      4011, 5010, 5011, 5012, 5013, 6012, 6013, 6014, 6015, 7013, 7014, 7015, 7016, 7017,  &
      8015, 8016, 8017, 8018, 8019, 9018, 9019, 9020, 10020, 10021, 10022, 10023, 11022,  &
      11023, 11024, 11424, 12024, 12025, 12026, 13027, 13029, 14028, 14029, 14030, 14031, &
      14032, 15031, 15032, 15033, 15034, 16032, 16033, 16034, 16035, 16036, 16037, 17034, &
      17035, 17036, 17037, 17038, 17040, 17438, 18036, 18037, 18038, 18039, 18040, 18041, &
      18042, 19038, 19039, 19040, 19041, 19042, 19043, 19044, 20040, 20041, 20042, 20043, &
      20044, 20045, 20046, 20047, 20048, 21044, 21045, 21046, 21047, 21048, 21049, 21050, &
      21445, 21446, 22045, 22046, 22047, 22048, 22049, 22050, 22051, 23049, 23050, 23051, &
      23052, 23053, 23054, 24049, 24050, 24051, 24052, 24053, 24054, 24055, 24066, 24067, &
      25053, 25054, 25055, 25056, 25057, 25058, 25066, 25067, 25068, 25069, 26053, 26054, &
      26055, 26056, 26057, 26058, 26059, 26061, 26065, 26066, 26067, 26068, 26069, 26070, &
      26071, 26072, 27057, 27058, 27059, 27060, 27061, 27062, 27064, 27065, 27066, 27067, &
      27068, 27069, 27070, 27071, 27072, 27073, 27074, 27075, 27458, 27460, 28057, 28058, &
      28059, 28060, 28061, 28062, 28063, 28064, 28065, 28066, 28067, 28068, 28069, 28070, &
      28071, 28072, 28073, 28074, 28075, 28076, 28077, 28078, 29062, 29063, 29064, 29065, &
      29066, 29067, 29068, 29069, 29070, 29071, 29072, 29073, 29074, 29075, 29076, 29077, &
      29078, 29079, 29080, 29081, 29468, 29470, 30062, 30063, 30064, 30065, 30066, 30067, &
      30068, 30069, 30070, 30071, 30072, 30073, 30074, 30075, 30076, 30077, 30078, 30079, &
      30080, 30081, 30082, 30083, 30469, 30471, 31066, 31067, 31068, 31069, 31070, 31071, &
      31072, 31073, 31074, 31075, 31076, 31077, 31078, 31079, 31080, 31081, 31082, 31083, &
      31084, 31085, 31086, 31472, 31474, 32068, 32069, 32070, 32071, 32072, 32073, 32074, &
      32075, 32076, 32077, 32078, 32079, 32080, 32081, 32082, 32083, 32084, 32085, 32086, &
      32087, 32088, 32089, 32471, 32473, 32475, 32477, 32479, 32481, 33071, 33072, 33073, &
      33074, 33075, 33076, 33077, 33078, 33079, 33080, 33081, 33082, 33083, 33084, 33085, &
      33086, 33087, 33088, 33089, 33090, 33091, 33092, 33482, 34073, 34074, 34075, 34076, &
      34077, 34078, 34079, 34080, 34081, 34082, 34083, 34084, 34085, 34086, 34087, 34088, &
      34089, 34090, 34091, 34092, 34093, 34094, 34473, 34477, 34479, 34481, 34483, 35075, &
      35077, 35078, 35079, 35080, 35081, 35082, 35083, 35084, 35085, 35086, 35087, 35088, &
      35089, 35090, 35091, 35092, 35093, 35094, 35095, 35096, 35097, 35477, 35479, 35480, &
      35482, 35484, 36077, 36078, 36079, 36080, 36081, 36082, 36083, 36084, 36085, 36086, &
      36087, 36088, 36089, 36090, 36091, 36092, 36093, 36094, 36095, 36096, 36097, 36098, &
      36099, 36100, 36479, 36481, 36483, 36485, 37079, 37081, 37083, 37084, 37085, 37086, &
      37087, 37088, 37089, 37090, 37091, 37092, 37093, 37094, 37095, 37096, 37097, 37098, &
      37099, 37100, 37101, 37102, 37486, 37490, 38083, 38084, 38085, 38086, 38087, 38088, &
      38089, 38090, 38091, 38092, 38093, 38094, 38095, 38096, 38097, 38098, 38099, 38100, &
      38101, 38102, 38103, 38104, 38105, 38485, 38487, 39085, 39087, 39088, 39089, 39090, &
      39091, 39092, 39093, 39094, 39095, 39096, 39097, 39098, 39099, 39100, 39101, 39102, &
      39103, 39104, 39105, 39106, 39107, 39108, 39487, 39489, 39490, 39491, 39493, 39496, &
      39497, 39498, 40087, 40088, 40089, 40090, 40091, 40092, 40093, 40094, 40095, 40096, &
      40097, 40098, 40099, 40100, 40101, 40102, 40103, 40104, 40105, 40106, 40107, 40108, &
      40109, 40110, 40489, 40490, 41089, 41090, 41091, 41092, 41093, 41094, 41095, 41096, &
      41097, 41098, 41099, 41100, 41101, 41102, 41103, 41104, 41105, 41106, 41107, 41108, &
      41109, 41110, 41111, 41112, 41113, 41491, 41493, 41494, 41495, 41497, 41498, 41499, &
      41500, 41502, 41504, 42091, 42092, 42093, 42094, 42095, 42096, 42097, 42098, 42099, &
      42100, 42101, 42102, 42103, 42104, 42105, 42106, 42107, 42108, 42109, 42110, 42111, &
      42112, 42113, 42114, 42115, 42493, 43095, 43096, 43097, 43098, 43099, 43100, 43101, &
      43102, 43103, 43104, 43105, 43106, 43107, 43108, 43109, 43110, 43111, 43112, 43113, &
      43114, 43115, 43116, 43117, 43118, 43495, 43497, 43499, 43502, 44095, 44096, 44097, &
      44098, 44099, 44100, 44101, 44102, 44103, 44104, 44105, 44106, 44107, 44108, 44109, &
      44110, 44111, 44112, 44113, 44114, 44115, 44116, 44117, 44118, 44119, 44120, 45099, &
      45101, 45102, 45103, 45104, 45105, 45106, 45107, 45108, 45109, 45110, 45111, 45112, &
      45113, 45114, 45115, 45116, 45117, 45118, 45119, 45120, 45121, 45122, 45123, 45501, &
      45502, 45503, 45504, 45505, 45506, 45508, 45510, 46101, 46102, 46103, 46104, 46105, &
      46106, 46107, 46108, 46109, 46110, 46111, 46112, 46113, 46114, 46115, 46116, 46117, &
      46118, 46119, 46120, 46121, 46122, 46123, 46124, 46125, 46126, 46507, 46509, 46511, &
      47103, 47105, 47106, 47107, 47108, 47109, 47110, 47111, 47112, 47113, 47114, 47115, &
      47116, 47117, 47118, 47119, 47120, 47121, 47122, 47123, 47124, 47125, 47126, 47127, &
      47128, 47129, 47130, 47505, 47506, 47507, 47508, 47509, 47510, 47511, 47513, 47515, &
      47516, 47517, 47518, 47520, 47522, 48105, 48106, 48107, 48108, 48109, 48110, 48111, &
      48112, 48113, 48114, 48115, 48116, 48117, 48118, 48119, 48120, 48121, 48122, 48123, &
      48124, 48125, 48126, 48127, 48128, 48129, 48130, 48131, 48132, 48511, 48513, 48515, &
      48517, 48519, 48521, 49109, 49111, 49112, 49113, 49114, 49115, 49116, 49117, 49118, &
      49119, 49120, 49121, 49122, 49123, 49124, 49125, 49126, 49127, 49128, 49129, 49130, &
      49131, 49132, 49133, 49134, 49135, 49512, 49513, 49514, 49515, 49516, 49517, 49518, &
      49519, 49520, 49521, 49522, 49523, 49524, 49525, 49526, 49527, 49528, 49529, 49530, &
      49531, 49616, 49618, 49620, 49622, 49630, 50111, 50112, 50113, 50114, 50115, 50116, &
      50117, 50118, 50119, 50120, 50121, 50122, 50123, 50124, 50125, 50126, 50127, 50128, &
      50129, 50130, 50131, 50132, 50133, 50134, 50135, 50136, 50137, 50513, 50517, 50519, &
      50521, 50523, 50525, 50527, 50528, 50529, 50530, 50531, 51115, 51117, 51118, 51119, &
      51120, 51121, 51122, 51123, 51124, 51125, 51126, 51127, 51128, 51129, 51130, 51131, &
      51132, 51133, 51134, 51135, 51136, 51137, 51138, 51139, 51518, 51520, 51522, 51524, &
      51526, 51528, 51529, 51530, 51532, 51534, 51624, 51626, 52118, 52119, 52120, 52121, &
      52122, 52123, 52124, 52125, 52126, 52127, 52128, 52129, 52130, 52131, 52132, 52133, &
      52134, 52135, 52136, 52137, 52138, 52139, 52140, 52141, 52142, 52521, 52523, 52525, &
      52527, 52529, 52531, 52533, 53121, 53122, 53123, 53124, 53125, 53126, 53127, 53128, &
      53129, 53130, 53131, 53132, 53133, 53134, 53135, 53136, 53137, 53138, 53139, 53140, &
      53141, 53142, 53143, 53144, 53145, 53530, 53532, 53533, 53534, 53536, 54122, 54123, &
      54124, 54125, 54126, 54127, 54128, 54129, 54130, 54131, 54132, 54133, 54134, 54135, &
      54136, 54137, 54138, 54139, 54140, 54141, 54142, 54143, 54144, 54145, 54146, 54147, &
      54525, 54527, 54529, 54531, 54533, 54534, 54535, 55127, 55128, 55129, 55130, 55131, &
      55132, 55133, 55134, 55135, 55136, 55137, 55138, 55139, 55140, 55141, 55142, 55143, &
      55144, 55145, 55146, 55147, 55148, 55149, 55150, 55151, 55534, 55535, 55536, 55538, &
      56128, 56129, 56130, 56131, 56132, 56133, 56134, 56135, 56136, 56137, 56138, 56139, &
      56140, 56141, 56142, 56143, 56144, 56145, 56146, 56147, 56148, 56149, 56150, 56151, &
      56152, 56153, 56531, 56533, 56535, 56536, 56537, 57133, 57134, 57135, 57136, 57137, &
      57138, 57139, 57140, 57141, 57142, 57143, 57144, 57145, 57146, 57147, 57148, 57149, &
      57150, 57151, 57152, 57153, 57154, 57155, 57546, 58134, 58135, 58136, 58137, 58138, &
      58139, 58140, 58141, 58142, 58143, 58144, 58145, 58146, 58147, 58148, 58149, 58150, &
      58151, 58152, 58153, 58154, 58155, 58156, 58157, 58537, 58539, 59139, 59140, 59141, &
      59142, 59143, 59144, 59145, 59146, 59147, 59148, 59149, 59150, 59151, 59152, 59153, &
      59154, 59155, 59156, 59157, 59158, 59159, 59542, 59544, 59548, 60140, 60141, 60142, &
      60143, 60144, 60145, 60146, 60147, 60148, 60149, 60150, 60151, 60152, 60153, 60154, &
      60155, 60156, 60157, 60158, 60159, 60160, 60161, 60541, 61141, 61143, 61144, 61145, &
      61146, 61147, 61148, 61149, 61150, 61151, 61152, 61153, 61154, 61155, 61156, 61157, &
      61158, 61159, 61160, 61161, 61162, 61163, 61548, 61552, 61554, 61652, 62143, 62144, &
      62145, 62146, 62147, 62148, 62149, 62150, 62151, 62152, 62153, 62154, 62155, 62156, &
      62157, 62158, 62159, 62160, 62161, 62162, 62163, 62164, 62165, 63147, 63149, 63150, &
      63151, 63152, 63153, 63154, 63155, 63156, 63157, 63158, 63159, 63160, 63161, 63162, &
      63163, 63164, 63165, 63166, 63167, 63552, 63554, 63652, 64149, 64150, 64151, 64152, &
      64153, 64154, 64155, 64156, 64157, 64158, 64159, 64160, 64161, 64162, 64163, 64164, &
      64165, 64166, 64167, 64168, 64169, 65151, 65153, 65155, 65156, 65157, 65158, 65159, &
      65160, 65161, 65162, 65163, 65164, 65165, 65166, 65167, 65168, 65169, 65170, 65171, &
      65556, 65558, 66154, 66155, 66156, 66157, 66158, 66159, 66160, 66161, 66162, 66163, &
      66164, 66165, 66166, 66167, 66168, 66169, 66170, 66171, 66172, 66565, 67159, 67160, &
      67161, 67162, 67163, 67164, 67165, 67166, 67167, 67168, 67169, 67170, 67171, 67172, &
      67559, 67561, 67562, 67563, 67564, 67566, 67570, 68160, 68161, 68162, 68163, 68164, &
      68165, 68166, 68167, 68168, 68169, 68170, 68171, 68172, 68567, 69165, 69166, 69167, &
      69168, 69169, 69170, 69171, 69172, 69173, 70166, 70167, 70168, 70169, 70170, 70171, &
      70172, 70173, 70174, 70175, 70176, 70177, 70569, 70575, 71169, 71171, 71172, 71173, &
      71174, 71175, 71176, 71177, 71178, 71179, 71180, 71569, 71571, 71572, 71576, 71577, &
      72171, 72172, 72174, 72175, 72176, 72177, 72178, 72179, 72180, 72181, 72182, 72577, &
      72580, 73179, 73180, 73181, 73182, 73183, 73580, 74180, 74181, 74182, 74183, 74184, &
      74186, 74583, 75185, 75187, 77191, 77193, 78197, 78198, 78199, 78201, 79197, 79198, &
      79199, 79200, 79201, 79202, 79204, 80196, 80198, 80199, 80200, 80201, 80202, 80203, &
      80204, 80205, 80206, 81202, 81203, 81204, 81205, 81206, 81207, 81208, 81209, 81210, &
      82202, 82203, 82204, 82205, 82206, 82207, 82208, 82209, 82210, 82211, 82212, 82213, &
      82214, 82607, 83207, 83208, 83209, 83210, 83211, 83212, 83213, 83214, 83215, 83610, &
      84210, 84211, 84212, 84213, 84214, 84215, 84216, 84217, 84218, 84611, 85215, 85216, &
      85217, 85218, 85219, 86216, 86217, 86218, 86219, 86220, 86221, 86222, 86223, 87219, &
      87220, 87221, 87222, 87223, 88220, 88221, 88222, 88223, 88224, 88225, 88226, 88227, &
      88228, 88229, 89223, 89224, 89225, 89226, 89227, 89228, 89229, 90225, 90226, 90227, &
      90228, 90229, 90230, 90231, 90232, 90233, 90234, 90235, 91229, 91230, 91231, 91232, &
      91233, 91234, 91235, 91634, 92230, 92231, 92232, 92233, 92234, 92235, 92236, 92237, &
      92238, 92239, 92240, 92241, 92242, 92635, 93233, 93234, 93235, 93236, 93237, 93238, &
      93239, 93240, 93241, 93242, 93636, 93640, 94234, 94235, 94236, 94237, 94238, 94239, &
      94240, 94241, 94242, 94243, 94244, 94245, 94246, 94247, 94637, 95239, 95240, 95241, &
      95242, 95243, 95244, 95245, 95246, 95247, 95642, 95644, 95646, 96239, 96240, 96241, &
      96242, 96243, 96244, 96245, 96246, 96247, 96248, 96249, 96250, 96251, 97247, 97248, &
      97249, 97250, 97251, 98247, 98248, 98249, 98250, 98251, 98253, 98254, 99253, 99254, &
      99255, 100255/)

      SNF_ISO_TYPE = &
      (/ 0,  0,  0,  2,  0,  1,  2,  3,  1,  1,  0,  1,  2,  0,  1,  &
       3,  0,  1,  1,  0,  0,  0,  1,  1,  3,  0,  1,  1,  1,  0,  &
       0,  0,  1,  1,  0,  1,  0,  0,  3,  1,  3,  3,  1,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  2,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  3,  &
       2,  1,  1,  1,  1,  0,  1,  2,  1,  1,  1,  1,  1,  1,  1,  &
       0,  0,  1,  2,  1,  3,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  0,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  0,  1,  0,  1,  0,  2,  2,  2,  0,  2,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  &
       1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  0,  1,  3,  1,  2,  2,  0,  1,  0,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  0,  1,  0,  2,  3,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  2,  2,  0,  0,  1,  &
       0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  2,  1,  0,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  &
       2,  2,  2,  0,  0,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  1,  0,  1,  &
       0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  0,  2,  0,  2,  2,  0,  0,  &
       2,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  0,  1,  1,  1,  0,  2,  0,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       3,  3,  3,  2,  0,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  3,  &
       0,  0,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  3,  1,  3,  2,  2,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  3,  0,  &
       0,  0,  3,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  3,  1,  2,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  0,  3,  0,  0,  3,  3,  3,  2,  2,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  3,  1,  1,  2,  1,  2,  2,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  &
       1,  1,  1,  1,  2,  2,  3,  3,  0,  2,  2,  1,  2,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  2,  1,  2,  1,  2,  1,  2,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  &
       1,  3,  2,  3,  2,  2,  0,  3,  1,  0,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  &
       1,  1,  1,  3,  1,  2,  1,  3,  1,  2,  2,  2,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  3,  3,  3,  2,  3,  3,  &
       0,  3,  0,  1,  0,  1,  0,  1,  2,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  &
       1,  0,  1,  2,  3,  1,  3,  2,  3,  1,  3,  2,  2,  2,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  3,  1,  3,  2,  0,  &
       2,  1,  2,  1,  0,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  2,  1,  2,  0,  &
       3,  1,  3,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  2,  0,  2,  &
       0,  3,  2,  3,  3,  1,  0,  2,  0,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  &
       1,  1,  1,  2,  2,  0,  1,  2,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  &
       2,  0,  0,  0,  2,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  2,  2,  &
       1,  3,  0,  0,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  2,  2,  0,  &
       1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  2,  2,  1,  1,  2,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  3,  1,  2,  3,  &
       0,  0,  0,  2,  0,  1,  2,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  3,  1,  1,  1,  2,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  &
       1,  0,  1,  3,  3,  0,  2,  2,  2,  2,  2,  0,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  3,  1,  0,  0,  0,  3,  2,  &
       2,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  &
       1,  3,  3,  3,  3,  1,  3,  2,  2,  0,  0,  1,  2,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  1,  0,  0,  3,  1,  &
       2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  &
       0,  1,  3,  2,  2,  2,  2,  3,  0,  2,  1,  1,  1,  1,  1,  &
       1,  1,  1,  3,  2,  1,  1,  0,  1,  2,  1,  1,  1,  1,  1,  &
       1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  &
       2,  0,  0,  2,  0,  1,  1,  1,  1,  2,  2,  1,  2,  1,  0,  &
       1,  1,  1,  1,  3,  2,  2,  2,  2,  3,  3,  3,  3,  3,  1,  &
       1,  3,  2,  2,  1,  1,  3,  0,  3,  1,  3,  3,  1,  1,  1,  &
       1,  1,  3,  0,  0,  1,  2,  2,  3,  2,  2,  0,  2,  0,  1,  &
       1,  0,  2,  2,  0,  2,  3,  0,  3,  2,  3,  0,  0,  0,  0,  &
       0,  0,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  1,  0,  0,  &
       0,  0,  0,  0,  1,  0,  1,  0,  1,  2,  1,  2,  2,  0,  0,  &
       1,  0,  0,  3,  2,  1,  2,  2,  2,  1,  2,  0,  0,  1,  0,  &
       1,  1,  1,  2,  0,  0,  0,  1,  0,  0,  1,  0,  0,  3,  1,  &
       1,  0,  0,  1,  0,  1,  0,  1,  1,  0,  0,  1,  1,  1,  0,  &
       0,  1,  0,  1,  1,  1,  2,  1,  0,  3,  1,  1,  0,  0,  3,  &
       0,  1,  1,  1,  1,  1,  2,  1,  2,  1,  1,  1,  2,  0,  0,  &
       0,  0,  1,  0,  1,  0,  1,  1,  1,  2,  1,  0,  2,  1,  0,  &
       2,  1,  0,  3,  3,  2,  3,  1,  0,  1,  2,  1,  1,  1,  3,  &
       1,  2,  2,  2,  2,  1,  1,  1,  1,  1,  1,  3,  1,  3,  2,  &
       3,  3,  0,  3,  3,  1,  2,  1,  2,  1,  1,  1,  1,  2,  1,  &
       0,  1,  2,  0,  1,  0,  1,  1,  1,  3,  3,  1,  0,  0,  0,  &
       2,  2,  3,  1,  0,  1,  0,  0,  2,  1,  1,  0,  0,  0,  3,  &
           0,      0,      0,      0,      0,      0,      0/)

      SNF_PRE_1 = &
      (/  0,  0,  0,   1003,  0,  0,   2006,   2008,  0,  0,  0,  0,   3009,  0,  0, &
       4010,  0,  0,  0,  0,  0,  0,  0,  0,   6014,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,   9022,  0,  12022,  10023,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  21050,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  25052, &
      23053,  0,  0,  0,  0,  0,  0,  24055,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  25056,  0,  25058,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  26059,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  27460,  27061,  29062,  0,  29064,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  30063, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  29064,  0,  31066,  29067,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  31472,  30073,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  33071,  33072,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  34073,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  33074, &
      35075,  33076,  33077,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  36077,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  37079,  0,  37081,  35082,  0,  0, &
      36485,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  37486,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      37084,  39085,  37086,  39487,  0,  0,  37490,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  38089,  38090,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  39091, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  41491,  0,  40093,  41494,  40095,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  42493,  0, &
          0,  0,  41097,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  43497,  0,  42099,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  45097,  0,  0,  43100,  43101,  43502,  43103,  43104, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  45499,  0,  0,  44103,  0,  44105,  44106,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  45102,  47103,  45504,  47105,  0,  45107,  45508,  0,  45510,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  48105,  0,  48107,  0,  46109,  0,  46511,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  47108,  49109,  47510,  47111,  47112,  0,  49114,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  49511,  0,  49513,  0,  48115,  0,  49517,  48118,  49519,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  49112,  50513,  49114,  49515,  48116,  49117, &
          0,  49519,  0,  0,  0,  0,  0,  0,  49527,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  52118,  51519,  0,  50121,  51522,  50123,  0,  50125,  51526,  50127,  51528,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  52521,  0,  52523,  51124,  0, &
      51126,  0,  51128,  0,  0,  0,  51532,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  54123,  0,  54125,  0, &
      52127,  0,  52129,  53530,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  54525,  0,  55127, &
          0,  54529,  53130,  53131,  53132,  0,  0,  53135,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  56128,  56129, &
          0,  0,  0,  54133,  55534,  0,  0,  54137,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  55130,  56531, &
      55132,  0,  0,  0,  55136,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  58134,  58135, &
          0,  58137,  0,  0,  56140,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  57140,  57141,  0, &
          0,  57144,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  60140,  58141,  0,  0,  58144,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  61140,  0,  59142,  59143, &
          0,  0,  0,  59147,  0,  0,  59150,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  62143,  0,  0,  0,  60147,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  61146,  61147,  0,  61149,  61150,  61151,  63552,  61153,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  64147,  0,  0,  0,  0,  62153,  63554, &
      62155,  62156,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  65149,  65550,  65151,  63152,  0,  65154,  65155,  63156,  0,  0,  0,  63160,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  66155,  0,  0,  0,  64159,  0, &
      64161,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  67558,  67159,  65160,  65161,  67562,  65163,  0,  65165,  0,  0,  0,  0,  0, &
          0,  0,  0,  67559,  68160,  0,  0,  0,  0,  68165,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      67166,  0,  0,  67169,  0,  0,  0,  0,  0,  70166,  70167,  0,  68169,  0,  0, &
          0,  0,  0,  0,  71168,  71169,  69170,  71171,  69172,  69173,  69174,  69175,  69176,  69177,  0, &
          0,  71569,  72171,  71572,  0,  0,  70175,  0,  70177,  0,  70179,  70180,  0,  0,  0, &
          0,  0,  72571,  0,  0,  0,  71576,  71177,  71578,  73179,  72580,  0,  74186,  0,  0, &
          0,  0,  74181,  72182,  0,  73580,  74580,  0,  75182,  73183,  73184,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  80203,  0,  82205,  83610,  0,  0, &
          0,  0,  0,  82603,  81204,  0,  84210,  83207,  81208,  0,  81210,  0,  0,  0,  0, &
          0,  0,  0,  82209,  0,  0,  0,  0,  0,  0,  0,  0,  0,  83212,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
          0,  0,  0,  0,  0,  0,  89225,  0,  0,  87220,  0,  0,  0,  0,  87225, &
          0,  0,  0,  0,  0,  0,  88225,  0,  88227,  0,  0,  0,  92230,  0,  0, &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  93235,  0,  0,  91634,  0,  0, &
      91230,  0,  0,  91233,  91634,  92635,  93236,  0,  0,  0,  94244,  0,  0,  0,  94233, &
          0,  94235,  95240,  92237,  95242,  0,  0,  0,  0,  0,  0,  95234,  0,  93236,  94637, &
      95238,  93239,  0,  93241,  93642,  0,  96248,  0,  96250,  0,  0,  0,  0,  96241,  0, &
          0,  0,  97249,  0,  0,  0,  0,  0,  0,  97240,  97241,  0,  0,  0,  0, &
      95246,  95247,  95248,  0,  0,  0,  0,  0,  96249,  0,  0,  0,  0,  0,  97250, &
          0,  0,  0,  0,  0,  0,  0/)

      SNF_PRE_2 = &
      (/  0,  0,  0,  0,  0,  0,  0,   4007,  0,  0,  0,  0,  0,  0,  0,  &
       6010,  0,  0,  0,  0,  0,  0,  0,  0,   8014,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  11022,  0,  13023,  12023,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  25452,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  25458,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  31064,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  30473,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      39484,  39485,  39486,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  41091,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  42091,  0,  41493,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  43093,  0,  &
          0,  0,  43497,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  44097,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  45497,  0,  0,  45100,  45501,  45102,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  46099,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  47104,  47505,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  49508,  0,  49510,  0,  0,  0,  49514,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  50111,  0,  0,  0,  48515,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  51112,  51113,  51114,  0,  51116,  50517,  &
          0,  51119,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  52119,  0,  50521,  0,  50523,  0,  50525,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  53121,  0,  53123,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
      54127,  0,  52529,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  55129,  0,  55131,  55132,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  58537,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  61540,  0,  0,  61143,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  62543,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  64150,  64151,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  65151,  0,  0,  0,  0,  64153,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  65549,  66154,  65551,  65552,  0,  65154,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  66159,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  67158,  0,  0,  0,  0,  67163,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  68159,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  71568,  0,  0,  0,  0,  71173,  71574,  70575,  71576,  70577,  0,  &
          0,  72169,  0,  0,  0,  0,  72175,  0,  71577,  0,  71579,  71580,  0,  0,  0,  &
          0,  0,  73171,  0,  0,  0,  0,  0,  73178,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  76184,  0,  76186,  0,  75584,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  83203,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  84612,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  90224,  0,  0,  0,  0,  90229,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  93233,  93234,  0,  93636,  0,  0,  0,  0,  0,  0,  0,  95237,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  96238,  0,  93636,  0,  &
      96242,  96243,  0,  96245,  95242,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
          0,  0,  0,  0,  0,  0,  0,  0,  0,  98244,  98245,  0,  0,  0,  0,  &
          0,  0,  97648,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 100254,  &
          0,  0,  0,  0,  0,  0,  0/)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_DATA] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (.not.allocated(SNF_PCF_LIST)) allocate(SNF_PCF_LIST(SNF_N_ISO,3))
      SNF_PCF_LIST = 0
      mid_iso = 1
      DO ii = 1, SNF_N_ISO
         DO WHILE (int(SNF_ISO_ARRAY(ii))>int(SNF_ISO_LIST(mid_iso)))
             mid_iso = mid_iso+1
             if (mid_iso > 1567) go to 777
         END DO

         IF(int(SNF_ISO_ARRAY(ii)) == int(SNF_ISO_LIST(mid_iso))) THEN
            SNF_PCF_LIST(ii,1) = SNF_ISO_TYPE(mid_iso)
            SNF_PCF_LIST(ii,2) = SNF_PRE_1(mid_iso)
            SNF_PCF_LIST(ii,3) = SNF_PRE_2(mid_iso)
         ENDIF
      ENDDO

777   CONTINUE

      RETURN
      END SUBROUTINE SNF_DATA


      SUBROUTINE SNF_COOLING
      IMPLICIT NONE

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_COOLING] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      CALL SNF_MK_INPUT
      END SUBROUTINE SNF_COOLING


      SUBROUTINE SNF_MK_INPUT
      USE Inc_File,         only: Len_INP, Name_inp
      USE Inc_Geometry,     only: meshsize_z !, meshsize_x !, GridSize_x
      USE Inc_History,      only: Hist_T_Mod_FA, Hist_T_Fuel_FA, hist_PPM
      USE Inc_Pinvar,       only: npin
      USE Inc_FA,           only: N_pin
      IMPLICIT NONE
      INTEGER(4)                                  :: ii, kk
      REAL(8)                                     :: mid, corr
      LOGICAL(1)                                  :: FLAG_OPEN
      INTEGER(4)                                  :: W_COOLING_INPUT=9006
      REAL(8)                                     :: dT
      CHARACTER(300)                              :: NAME_cooling_input
      REAL(8)                                     :: vol_coeff = 6.36173E-01 ! per cm


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SNF_MK_INPUT] in Mod_SNF'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      dT = SNF_DECAY_YEAR*3600d0*24d0*365d0 ! cooling : year
      Name_cooling_input = Name_INP( 1 : Len_INP )
      Name_cooling_input( Len_INP + 1 : Len_INP + 10 ) = "_sttmp.inp"
      DO WHILE (.TRUE.)
         INQUIRE( UNIT = W_COOLING_INPUT, OPENED = Flag_Open )
         IF ( .NOT. Flag_Open ) THEN
            OPEN(W_COOLING_INPUT, FILE = NAME_cooling_input)
            EXIT
         ELSE
         W_COOLING_INPUT = W_COOLING_INPUT + 10
         END IF
      END DO

      write(W_COOLING_INPUT,'(A5)') 'SYM 8'
      write(W_COOLING_INPUT,'(A12)') 'SRCTRM MERGE'
      write(W_COOLING_INPUT,'(A8,F9.5)') 'DECAY 1 ',SNF_DECAY_YEAR
      if (SNF_ENDF == 0) THEN
         write(W_COOLING_INPUT,'(A11)') 'DEPLIB E6.8'
         write(W_COOLING_INPUT,'(A11)') 'XSLIB  E6.8'
      elseif (SNF_ENDF==1) THEN
         write(W_COOLING_INPUT,'(A11)') 'DEPLIB E7.0'
         write(W_COOLING_INPUT,'(A11)') 'XSLIB  E7.0'
      elseif (SNF_ENDF==2) THEN
         write(W_COOLING_INPUT,'(A11)') 'DEPLIB E7.1'
         write(W_COOLING_INPUT,'(A11)') 'XSLIB  E7.1'
      elseif (SNF_ENDF==3) THEN
         write(W_COOLING_INPUT,'(A14)') 'DEPLIB E71JD40'
         write(W_COOLING_INPUT,'(A14)') 'XSLIB  E71JD40'
      endif
      write(W_COOLING_INPUT,'(A10)') 'ON RESRING'
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A22)') 'DEP 0.0008 1 1 w/g DAY '
      write(W_COOLING_INPUT,'(A1)') '1'
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A5,F7.2,F7.2)') 'TEMP ',Hist_T_Fuel_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K), Hist_T_Mod_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
      write(W_COOLING_INPUT,'(A5,F7.2,F13.2)') 'H2O  ',Hist_T_Mod_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K), hist_PPM(SNF_N_END_CYBU)
      write(W_COOLING_INPUT,'(A6,1x,F9.2)') 'HEIGHT', meshsize_z(SNF_K)
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A9,F7.3,A5,F7.2)') 'FA 1 1 1 ', 1.5d0, ' MOD ', Hist_T_Mod_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
      write(W_COOLING_INPUT,'(A1)') '1'
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A5)') 'Pin 1'
      write(W_COOLING_INPUT,'(F8.4,2x,A4)') 0.4500d0, 'FUE1'
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A)') 'MAT FUE1 BURN'
      ! VOL per 1cm :: 6.36173E-01
      corr = 1d0
      if (.not.flag_snf_vol) then
         if (flag_snf_pin) then
            DO ii = 1, SNF_N_ISO
               kk = SNF_MIN_MAT
               mid = 0d0
               do kk = 1,int(npin*npin)
                  if (SNF_VOL_RATIO(kk,1) .EQ. 0) cycle
                  mid = mid + SNF_RESULT_ARRAY(kk,ii)
               enddo
               if (flag_snf_coo) corr = npin*npin/N_pin
               write(W_COOLING_INPUT,'(I7, ES13.5)') SNF_ISO_ARRAY(ii), corr*mid /(vol_coeff * meshsize_z(SNF_K))  !/n_pin ! unit is barn-cm
            ENDDO
         else
            DO ii = 1, SNF_N_ISO
               kk = SNF_MIN_MAT
               !write(W_COOLING_INPUT,'(I7, ES13.5)') SNF_ISO_ARRAY(ii), SNF_RESULT_ARRAY(kk,ii) /SNF_VOL_RATIO(kk,SNF_N_BU_UP) !/n_pin ! unit is barn-cm
               write(W_COOLING_INPUT,'(I7, ES13.5)') SNF_ISO_ARRAY(ii), SNF_RESULT_ARRAY(kk,ii) / (vol_coeff * meshsize_z(SNF_K))
            ENDDO
         endif
      else
         if (flag_snf_pin) then
            DO ii = 1, SNF_N_ISO
               kk = SNF_MIN_MAT
               mid = 0d0
               do kk = 1,int(npin*npin)
                  mid = mid + SNF_RESULT_ARRAY(kk,ii)
               enddo
               write(W_COOLING_INPUT,'(I7, ES13.5)') SNF_ISO_ARRAY(ii), mid /SNF_VOL !/n_pin ! unit is barn-cm
            ENDDO
         else
            DO ii = 1, SNF_N_ISO
               kk = SNF_MIN_MAT
               write(W_COOLING_INPUT,'(I7, ES13.5)') SNF_ISO_ARRAY(ii), SNF_RESULT_ARRAY(kk,ii) /SNF_VOL !/n_pin ! unit is barn-cm
            ENDDO
         endif
      endif
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A)') 'MAT MOD'
      write(W_COOLING_INPUT,'(A3,F7.2)') 'COO',Hist_T_Mod_FA(SNF_N_END_CYBU,SNF_Ixy_1N,SNF_K)
      write(W_COOLING_INPUT,'(A)') ''
      write(W_COOLING_INPUT,'(A)') 'MAT GAP'
      write(W_COOLING_INPUT,'(A17)') '2004 2.68714E-05'
!      stop

      RETURN
      END SUBROUTINE SNF_MK_INPUT
#endif

      END MODULE Mod_SNF

#endif
