#ifdef tuan_fr_crm
MODULE Mod_CRM_HEX

    USE Inc_INP
    USE Inc_Constant
    USE Inc_Geometry
    USE Inc_RP
    USE Mod_Alloc

    IMPLICIT NONE
    CONTAINS
    
    SUBROUTINE CRM_HEX
    USE Inc_CR
    USE Inc_RP 
    USE Inc_XS_File 
    use Inc_Geometry, only: I_1Nto4N, I_4Nto1N
    use Inc_Option,        only: n_group
    use Inc_TPEN
    use Inc_maXS, only: maXS_R_3D, maXS_SCAT_3D, maXs_s_3D

    IMPLICIT NONE
    
    IF (.not. Flag_Card_CR) RETURN
    
    
    
    
    
    
    END SUBROUTINE CRM_HEX

END MODULE Mod_CRM_HEX

      
    
    
#endif