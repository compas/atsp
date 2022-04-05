      PROGRAM ZEEMAN
      IMPLICIT NONE
      INTEGER I, J, N_CFG, N_JVALUE, N_EIGVEC
      INTEGER, DIMENSION(:), ALLOCATABLE :: DOMINANT_CSF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G_LSJ, WEIGHT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: L, S, JQ
      CHARACTER*132, DIMENSION(:), ALLOCATABLE :: CONFIG
      CHARACTER*24 NAME, CFILE, JFILE, GJFILE
      !
      ! Open required files
      !
      WRITE(*,*) ' Name of the file'
      READ(*,'(A)') NAME
      J = INDEX(NAME,' ')
      IF (J==1) THEN
         WRITE(*,*) ' Name may not start with blanks'
         STOP
      END IF
      
      CFILE = NAME(1:J-1)//'.c'
      JFILE = NAME(1:J-1)//'.j'
      GJFILE = NAME(1:J-1)//'.gj'
      
      OPEN(UNIT=11,FILE=CFILE,STATUS='OLD')
      OPEN(UNIT=12,FILE=JFILE,STATUS='OLD')
      OPEN(UNIT=13,FILE=GJFILE,STATUS='UNKNOWN')
      
      CALL NEIGVEC(12,N_CFG,N_JVALUE,N_EIGVEC)                
      ! Determine number of configurations, J-values and eigenvectors
      ALLOCATE (L(N_CFG), S(N_CFG),JQ(N_EIGVEC))              
      ! Allocate memory for the L and S values of the CSFs in name.c          
      ALLOCATE (G_LSJ(N_CFG,N_EIGVEC))                        
      ! Allocate memory for the g factors of the LSJ coupled CSFs
      ALLOCATE (WEIGHT(N_CFG,N_EIGVEC))                       
      ! Allocate memory for the expansion coefficients (WEIGHT)
      ALLOCATE (CONFIG(N_EIGVEC))                             
      ! Allocate memory for the string containing the ......
      ALLOCATE (DOMINANT_CSF(N_EIGVEC))                       
      ! Allocate memory for the dominant CSF in each eigenvector 
      
      CALL READWEIGHT(12,WEIGHT,JQ,CONFIG,DOMINANT_CSF)       
      ! Read the expansioncoefficienst for the eigenvectors
     
      CALL READLS(11,N_CFG,L,S)                               
      ! Read the LS values of the CSFs in name.c
      CALL GLSJ(N_CFG,N_EIGVEC,L,S,JQ,G_LSJ)                  
      ! Compute the g factors for the CSFs
      CALL GJ(13,N_CFG,N_EIGVEC,JQ,WEIGHT,G_LSJ,CONFIG,DOMINANT_CSF)       
      ! Compute and print the g factor for the eigenvectors
      DEALLOCATE(WEIGHT,G_LSJ,L,S,JQ,DOMINANT_CSF)
      
      CONTAINS 
      
      !
      ! This subroutine determines the number of configurations (N_CFG), 
      ! J-values (N_JVALUE) 
      ! and eigenvectors (N_EIGVEC)
      !
      SUBROUTINE NEIGVEC(IU,N_CFG,N_JVALUE,N_EIGVEC)
      IMPLICIT NONE
      INTEGER :: IU, I, J, J2, N_CFG, NJ,  N_JVALUE, N_EIGVEC, N_LINES, NVEC
      INTEGER :: STATUS = 0
      N_JVALUE = 0
      N_EIGVEC = 0
      READ(IU,'(40X,I6)') N_CFG
      N_LINES = (N_CFG-1)/7+1
      DO 
         READ(IU,'(//8X,I4,10X,I4)',IOSTAT=STATUS) J2,NJ
         IF (STATUS<0) EXIT 
         N_JVALUE = N_JVALUE + 1   
         N_EIGVEC = N_EIGVEC + NJ
         DO I = 1,NJ
            READ(IU,'(A)')
            READ(IU,'(A)')
            DO J = 1,N_LINES
               READ(IU,'(A)')
            END DO
         END DO
      END DO
      REWIND(UNIT=IU)
      END SUBROUTINE NEIGVEC
      
      ! 
      ! This subroutine reads the configuration weights for all eigenvectorse
      !
      SUBROUTINE READWEIGHT(IU,WEIGHT,JQ,CONFIG,DOMINANT_CSF)
      IMPLICIT NONE
      INTEGER :: IU, I, K, J2, N_CFG, NJ
      INTEGER :: J = 0, STATUS = 0
      INTEGER, DIMENSION(:) :: DOMINANT_CSF
      DOUBLE PRECISION, DIMENSION(:,:) :: WEIGHT
      DOUBLE PRECISION, DIMENSION(:) :: JQ
      CHARACTER*132, DIMENSION(:) :: CONFIG
      READ(IU,'(40X,I6)') N_CFG
      DO 
         READ(IU,'(//8X,I4,10X,I4)',IOSTAT=STATUS) J2,NJ
         IF (STATUS<0) EXIT 
         DO I = 1,NJ
            J = J + 1
            JQ(J) = DBLE(J2)/2.D0                            
      ! J-value of eigvect number J
            READ(IU,'(A)')
            READ(IU,'(I6,18X,A)') DOMINANT_CSF(J),CONFIG(J)
            READ(IU,'(7F11.8)') (WEIGHT(K,J), K=1,N_CFG) 
         END DO
      END DO
      END SUBROUTINE READWEIGHT
      
      SUBROUTINE READLS(IU,N_CFG,L,S)
      IMPLICIT NONE
      INTEGER :: IU, N_CFG, L_STRING
      DOUBLE PRECISION :: MAP(68:83)
      DOUBLE PRECISION, DIMENSION(:) :: L, S
      CHARACTER*132 LINE
      !
      ! Define a map from the ASCII code of the spectroscopic notation 
      ! to the L quantum number
      !
      MAP(83) = 0.D0          ! S -> L = 0
      MAP(80) = 1.D0          ! P -> L = 1
      MAP(68) = 2.D0          ! D -> L = 2
      MAP(70) = 3.D0          ! F -> L = 3
      MAP(71) = 4.D0          ! G -> L = 4
      MAP(72) = 5.D0          ! H -> L = 5
      MAP(73) = 6.D0          ! I -> L = 6
      MAP(75) = 7.D0          ! K -> L = 7
      MAP(76) = 8.D0          ! L -> L = 8
      MAP(77) = 9.D0          ! M -> L = 9
      MAP(78) = 10.D0         ! N -> L = 10
      
      READ(IU,'(A)')
      READ(IU,'(A)')
      DO I = 1,N_CFG
         READ(IU,'(/A)') LINE
         L_STRING = LEN_TRIM(LINE)
         IF (IACHAR(LINE(L_STRING-1:L_STRING-1))>57) L_STRING = L_STRING - 1   
      ! To deal with configurations with only one group of equiv elec.
         S(I) = DBLE(IACHAR(LINE(L_STRING-1:L_STRING-1))-49)/2.D0
         L(I) = MAP(IACHAR(LINE(L_STRING:L_STRING)))
      END DO
      END SUBROUTINE READLS
      
      SUBROUTINE GLSJ(N_CFG,N_EIGVEC,L,S,JQ,G_LSJ)
      IMPLICIT NONE
      INTEGER :: N_CFG, N_EIGVEC, I, K
      DOUBLE PRECISION :: G_S = 2.002319304386D0
      DOUBLE PRECISION, DIMENSION(:,:) :: G_LSJ
      DOUBLE PRECISION, DIMENSION(:) :: L, S, JQ
      DO I = 1,N_EIGVEC
         IF (JQ(I)>0.25D0) THEN       
      ! Compute the g_j factor for all eigenvektors with J > 0
            DO K = 1,N_CFG
               G_LSJ(K,I) = 1.D0 + (G_S-1.D0)*(JQ(I)*(JQ(I)+1.D0)+S(K)*(S(K)+1.D0)-L(K)*(L(K)+1.D0))/(2.D0*JQ(I)*(JQ(I)+1.D0))
            END DO
         END IF
      END DO
      END SUBROUTINE GLSJ
      
      SUBROUTINE GJ(IU,N_CFG,N_EIGVEC,JQ,WEIGHT,G_LSJ,CONFIG,DOMINANT_CSF)
      IMPLICIT NONE
      INTEGER :: IU, N_CFG, N_EIGVEC, I, K
      INTEGER, DIMENSION(:) :: DOMINANT_CSF
      DOUBLE PRECISION :: G_J
      DOUBLE PRECISION, DIMENSION(:) :: JQ
      DOUBLE PRECISION, DIMENSION(:,:) :: WEIGHT, G_LSJ
      CHARACTER*132, DIMENSION(:) :: CONFIG
      WRITE(IU,*) '     g_J           g_J (LS)      J           E        dominant CSF'
      WRITE(IU,*)
      DO I = 1,N_EIGVEC
         G_J = 0.D0 
         IF (JQ(I)>0.25D0) THEN
            DO K = 1,N_CFG
               G_J = G_J + WEIGHT(K,I)*WEIGHT(K,I)*G_LSJ(K,I)
            END DO
            WRITE(IU,'(F11.8,F17.8,3X,F4.1,5X,A,F11.8)') G_J, G_LSJ(DOMINANT_CSF(I),I), JQ(I), TRIM(CONFIG(I))
         END IF
      END DO
      END SUBROUTINE GJ
      
      
      END PROGRAM ZEEMAN
      
            
