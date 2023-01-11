
    MODULE Inc_H2O

    implicit none
    REAL(8)                                     :: CORE_MPA=15.5d0
    LOGICAL(1)                                  :: if_h2o_Table_done=.false.

    INTEGER(4),PARAMETER                        :: n_T=1000
    REAL(8),PARAMETER                           :: min_T=273.15d0                       ! K
    REAL(8),PARAMETER                           :: max_T=640.00d0                       ! K
    REAL(8),PARAMETER                           :: del_t=(max_t-min_t)/real((n_t-1), 8 )! K
    REAL(8),PARAMETER                           :: del_to=1.0d0/del_t

    REAL(8),ALLOCATABLE                         :: T_grid(:)                            ! K
    REAL(8),ALLOCATABLE                         :: den(:)                               ! g/cc
    REAL(8),ALLOCATABLE                         :: h(:)                                 ! J/kg
    REAL(8),ALLOCATABLE                         :: k(:)                                 ! W/m-K
    REAL(8),ALLOCATABLE                         :: mu(:)                                ! Pa-s
    REAL(8),ALLOCATABLE                         :: cp(:)                                ! J/kg-K

    INTEGER(4),PARAMETER                        :: n_h=1000
    REAL(8)                                     :: min_h=0.0d0                          ! J/kg
    REAL(8)                                     :: max_h=0.0d0                          ! J/kg
    REAL(8)                                     :: del_h=0.0d0                          ! J/kg
    REAL(8)                                     :: del_ho=0.0d0                         ! 1/(J/kg)
    REAL(8),ALLOCATABLE                         :: h_grid(:)                            ! J/kg
    REAL(8),ALLOCATABLE                         :: h_T(:)                               ! K
    REAL(8),ALLOCATABLE                         :: h_den(:)                             ! g/cc

    END MODULE Inc_H2O
