program mtc_simula
use NETCDF
!use OMP_LIB
implicit none
integer :: I, J, K, Z, M, N, ICASE, TIME, LS               ! LOOP COUNTERS
integer :: CASEN, CASEDN, REG
integer :: CLDERA = 0 , MILAN = 0, ADV=0, CNV=0  ! Flags
integer :: LH = 1 , XK0 = 0 , A20 = 0 , A30 = 0 , TM01 = 1 ! Flags
integer,parameter :: DATEINI  = 199701, DATEEND =199701 !DATEINI  = 194802, DATEEND =201512
integer,parameter :: YEARINI = INT(DATEINI/100)
integer,parameter :: YEAREND = INT(DATEEND/100)
integer,parameter :: MESINI  = DATEINI-YEARINI*100
integer,parameter :: MESEND  = DATEEND-YEAREND*100
integer,parameter :: NTIME = (YEAREND-YEARINI)*12+(MESEND-MESINI)+1
integer :: CURRENTDATE, CURRENTMES, CURRENTYEAR
integer, dimension(NTIME) :: MNTH
real 	:: t1, t2, TMX, TSPX, TPMX, OCEAN, FP
real 	:: F30, F30p, F31, F33, F32, F32P, F38, F34p, F35, F36, F34! LW func.
character (len = 7),  parameter :: OF='output/'
character (len = 15), parameter :: FORMAT1 = "(X,'|',23X,14A)"
character (len = 38), parameter :: FORMAT2 = "(1X,'+',16('-'),3X,A23,3X,16('-'),'+')"
character (len = 15), parameter :: FORMAT3 = "(X,'|',61X,'|')"
character (len = 79), parameter :: FORMAT4 = "(X,'|'2X,'MNTHI = ',I2,3X,'MNTH = ',I2,3X,'SEASON = ',I2,3X,' N = ',I3,12X,'|')"
character (len = 22), parameter :: FORMAT5 = "(X,A,4X,39A,18X,A)"
character (len = 79), parameter :: FORMAT6 = "(X,'|',2X,'-',X,A53,4X,'|')"
character (len = 79), parameter :: FORMAT7 = "(X,'|'2X,'YEAR = ',I4,2X,'REC = ',I3,37X,'|')"
integer, dimension(12) :: PREV = (/ 12,1,2,3,4,5,6,7,8,9,10,11 /)
integer, dimension(12) :: SEA  = (/ 4,4,1,1,1,2,2,2,3,3,3,4 /)
real, dimension(12) :: CYEAR  = (/ 31,59,90,120,151,181,212,243,273,304,334,360 /)
!============================ Constants parameters =============================
real, parameter :: B       = 33.8639 * 1000.* .0378
real, parameter :: C       = 299792458.        ! [M S^-1] Light speed
real, parameter :: CS      = 4.1855e3          ! [J Kg^-1 K^-1] Spec. Heat water
real, parameter :: CV      = 1.004e3           ! [J Kg^-1 K^-1]
real, parameter :: CC      = 0.8374e3          ! [J Kg^-1 K^-1] dry soil
real, parameter :: CPV     = 1.884e3           ! [J Kg^-1 K^-1]
real, parameter :: CVD     = 0.718e3           ! [J Kg^-1 K^-1]
real, parameter :: CVV     = 1.422e3           ! [J Kg^-1 K^-1]
real, parameter :: D2      = 6.959e-3          ! [W^-1] Mendoza & Oda (1998)
real, parameter :: GAMMA   = 6.5e-3            ! [DegC M^-1] Lapse rate
real, parameter :: G       = 9.81              ! [M S^-2] Grav. ac.
real, parameter :: H0      = 11e3              ! [M] Height at 266MB
real, parameter :: H5      = 5.5e3             ! [M] Height at 500MB
real, parameter :: H7      = 3e3               ! [M] Height at 700MB
real, parameter :: HC      = 2.5               ! [M] Continental soil deep
real, parameter :: HS0     = 60.               ! [M] Ocean mixed layer deep
real, parameter :: HPL     = 6.62606896e-34    ! [J S] Planck constant
real, parameter :: K1      = 4.56080611E-41
real, parameter :: K2      = 26.8380051
real, parameter :: K3      = -CV*0.003*1.!26.8380051
real, parameter :: K4      = 0.7*0.003*1.!4.05013077E-02
real, parameter :: KB      = 1.3806503e-23     ! [J K^-1] Boltzmann constant
real, parameter :: KAUST   = 5e6               ! [K S^-1] Austausch coeficient
real, parameter :: P0      = 2.66e4            ! [Pa] Pressure at 266MB
real, parameter :: P5      = 5.00e4            ! [Pa] Pressure at 500MB
real, parameter :: P7      = 7.00e4            ! [Pa] Pressure AT 700MB
real, parameter :: PI	   = 3.141592
real, parameter :: R       = 6.37e6            ! [M] Earth's radius
real, parameter :: RD      = 287.053           ! Gas constant for dry air.
real, parameter :: RW      = 461.50            ! Gas constant for water vapor.
real, parameter :: RHOC    = 1.43e3            ! [Kg M^-3] soil density
real, parameter :: RHOS    = 1.e3              ! [Kg M^-3] Water density
real, parameter :: S0      = 1367.             ! [W M^-2] Solar constant
real, parameter :: SIGMA   = 5.6704e-8         ! [W M^-2 K^-4] Steffan-Boltzmann
real, parameter :: T0      = 216.5             ! [K] Temperature AT 266MB
real, parameter :: TS0     = 288.
real, parameter :: TC0     = 263.5              ! [K] Cloud Temperature
!real, parameter :: DT      = 3600*24*31        ! [Seg] Time step
real, parameter :: DT      = 2.68e6            ! [Seg] Time step
!real, parameter :: UM      = 20.
!real, parameter :: VM      = 0.
!=========================== Calculated parameters =============================
!real, parameter :: A       = 5.425*H0/(4.*T0+GAMMA*H0)    ! [18.4331799]
real, parameter :: A       = 2.*H0/(4.*T0+GAMMA*H0)
real, parameter :: C1      = 2*PI*HPL*C**2                 ! [3.742E-16 W M^2]
real, parameter :: C2      = (HPL*C/KB)                    ! [1.4385E-2 M K]
real, parameter :: DS      = Cs * rhoS * hs0
real, parameter :: EPS     = 23.4*PI/180.
real, parameter :: EX      = 0.017
real, parameter :: GRB     = g / ( Rd * gamma )            ! [5.25767279]
real, parameter :: RHO0    = P0/ ( RD * T0 )               ! [0.47 KG M^-3]
!real, parameter :: TS0     = GAMMA*H0 + T0                ! [288 K] SURF. TEMP.
!real, parameter :: TM0     = ( T0 + TS0 ) / 2.            ! [258.75 K]
real, parameter :: TM0     = T0 + GAMMA*H0/2
!real, parameter :: TA0     = GAMMA*H0 + T0             ! [288 K] AIR SURF. TEMP.
real, parameter :: TA0     = TS0
real, parameter :: F82     = (CV*P0/G)*((TA0/T0)**(GRB) - 1)! [7295820.50] Cv*a0
real, parameter :: F83     = - F82
real, parameter :: TAP2TMP = 1. + A*GAMMA/2.           ! [1.0599] T'A to T'M
real, parameter :: TP2TMP  = 1. - A*GAMMA/2.           ! [0.94] T' to T'M
!============================ Grid parameters ==================================
integer, parameter :: NLON = 192, NLAT = 94  ! data resolution
real, dimension(NLON) :: LON				           ! data longitude
real, dimension(NLAT) :: LAT				           ! data latitude
character (len = 1 )  :: MONTH
real, dimension(:,:), allocatable :: LAND, TOPO
real, dimension(NLON*NLAT) :: LANDMASK , OCEANMASK
real, dimension(NLON*NLAT) :: CLD, QQ0, QQ, S, ALBS, A1, A2, A3, VAN, UM, VM, XK
real, dimension(NLON*NLAT) :: G2N, G3N, G5N, G2, G3, G5, G2DN, G3DN, G5DN, CLDDN, HUMN
real, dimension(NLON*NLAT) :: H700, H500, H250, ADV_U, ADV_V, CONV
real, dimension(NLON*NLAT) :: TPN, TPPN, TSPN, TSPP, TMPN, TMPP, TMPDN, TSPDN,TMN
real, dimension(NLON*NLAT) :: TFRONT, TMDN
real, dimension(NLON*NLAT) :: F71, F72, F73, F74, F75, F76, F77, F78
real, dimension(NLON*NLAT) :: F81, F89, F90
real, dimension(NLON*NLAT) :: F41, F42, F46, F48, F49
real, dimension(NLON,NLAT) :: F84, F85, COSFP, SINFP, TPM, TMNOR
real, dimension(NLON,NLAT) :: F81A, F89A, F90A
real, dimension(NLON,NLAT,NTIME) :: F81B, F89B, F90B
real, dimension(NLON,NLAT,NTIME) :: F78B, F77B, F76B, F75B, F73B, F72B, F71B
real, dimension(NLON,NLAT,NTIME) :: F78A, F77A, F76A, F75A, F73A, F72A, F71A
real, dimension(NLON,NLAT,NTIME) :: F48A, A1A, G5A
real, dimension(NLON,NLAT,NTIME) :: IT, INSO, ESDN, ETDN, ABSDN
real, dimension(NLON,NLAT,NTIME) :: H500N , H500A , H500B , H250N , H250A
real, dimension(NLON,NLAT,NTIME) :: T500, TPMN , TPMDN , TSN , TSDN, TPMAB, TSUP
real, dimension(NLON,NLAT,NTIME) :: AUX1, AUX2, AUX3, AUX4
real, dimension(NLON,NLAT,NTIME) :: AUX5, AUX6, AUX7, AUX8, AUX9, AUX10, AUX11
real, dimension(NLON,NLAT,NTIME) :: AUX12, AUX13, AUX14, AUX15, AUX16, P500
!! snow
INTEGER :: loop, flag
real :: eps1
real, dimension(NLON,NLAT,NTIME) :: TSDN2, ALBSNEW
!!
call WELLCOME
call cpu_time( t1 )
!============================= Loading climatic fields =========================
call CLIM_FILES
!============================= Calculate LW functions ==========================
call RADIATIVE(T0,TA0,TS0,F30,F30p,F31,F32,F32p,F33,F34,F34p,F35,F36,F38)
!================================= Initial parameters ==========================
write(*,FORMAT2) '                 '
CURRENTDATE = YEARINI*100+MESINI
CURRENTMES = MESINI
CURRENTYEAR = YEARINI
!write(*,*)  YEAREND-YEARINI, YEAREND,YEARINI,NTIME
!write(*,*)
!return
do N = 1 , NTIME
    REG = (YEARINI-1948)*12+MESINI+N-2
    if ( REG <1 ) then
    write(*,FORMAT2) 'Invalid Record Number'
    call exit(0)
    endif
    write(*,FORMAT3)
    write(*,FORMAT7) CURRENTYEAR, REG
    write(*,FORMAT4) PREV(CURRENTMES) , CURRENTMES , SEA(CURRENTMES) , N
    call CLIM_DATA( N , CURRENTMES , PREV(CURRENTMES) )
    do ICASE = 1 , 1
        CASEN  = 2. - ICASE                  ! Normal case flag
        CASEDN = ICASE - 1.                  ! Departure to normal case flag
        if ( N > 1 ) then
            if ( ICASE == 1 ) then
	        m=1
    	    do j = 1 , NLAT
    	        do i = 1 , NLON
    	        	TMN(M) = TPMN(i,j,N-1)+ TM0
    	        	m=m+1      	
    	        enddo
    	    enddo
            write(*,*) 'Usando Tm output MTC de MNTH=',PREV(CURRENTMES)
            endif
            if ( ICASE == 2 ) then
	    	    m=1
    		    do j = 1 , NLAT
    		        do i = 1 , NLON
    		        	TMDN(M)=TPMDN(i,j,N-1)
    		        	m=m+1
    		        enddo
    		    enddo
    		    write(*,*) 'Usando Tm anomalia del MTC del MNTH=',PREV(CURRENTMES)
    	    endif
        endif
        loop=0
200     continue
        do M = 1 , NLON*NLAT
!                  TPN(M)  = TPN(M) + GAMMA*(H250(M)-H0)
!                  TPMX    =  TP(M)*CASEN !+ T1(M)*CASEDN
            OCEAN   =  LANDMASK(M)!1-OCEANMASK(M)
            if ( N == 1 .and. CASEDN==1 ) then
			 TMDN(M)=0.
            if (M==1) write(*,*) 'Usando TM Anomalia nula'
            endif
            TMX     =  TMN(M)  + TMDN(M) *CASEDN
            TSPX    =  TSPN(M) + TSPDN(M)*CASEDN
            TSPP(M) =  TSPX-TS0  ! TSP prima = TSP - TS0
            !if ( N == 1 ) then
            TMPP(M) =  TMX- TM0  !TMP prima = T prima/(1-A*gamma/2)
            !if (M==1) write(*,*) 'Usando Tm previa Observada'
            !else
            !TMPP(M) = TMDN(M)
            !endif
            A1(M)   =  QQ0(M)*(1.-(1.-XK(M))*CLD(M))*(1.-ALBS(M))
            G2(M)   =  G2N(M) + G2DN(M)*CASEDN    !F44
            G3(M)   =  G3N(M) !+ G3DN(M)*CASEDN
            G5(M)   =  G5N(M) + ( ADV_U(M) + ADV_V(M) )*ADV*CASEN + CNV*CASEN*CONV(M)!+ G5DN(M)*CASEDN
            F41(M)  =  0.*OCEAN*VAN(M)*0.981*CASEDN*K4*B
            F42(M)  =  0.*HUMN(M)/0.981*F41(M)!+ (1-OCEAN)*CASEDN*(1-D7)*K3*VAN(M)*(1-HUMN(M))
            F46(M)  =  0.!( OCEAN*K2+(1-OCEAN)*K3 )*VAN(M)*CASEDN
            F48(M)  =  F34 + F34P*CLD(M) + A1(M)           ! Es
            F49(M)   =  -QQ0(M)*(1.-XK(M))*(1.-ALBS(M)) + F34P
            F71(M)  =  DS*OCEAN/DT - F36 + F41(M) + F46(M)
            F72(M)  =  ( TSPP(M)*DS*OCEAN/DT - G2(M) - G3(M) + F48(M) ) / F71(M)!Ts
            F73(M)  =  ( F35 - F42(M) + F46(M) ) / F71(M)           ! [-]
            F75(M)  =  F49(M) / F71(M)
            F76(M)  =  ( F41(M) + F46(M) ) / F71(M)          ! [-]
            F77(M)  =  ( F42(M) - F46(M) ) / F71(M)
            F78(M)  =  F32 + CLD(M)*F32P + F46(M) !+ ( F41(M) ) !+ D2*()     ! [W m^-2 K^-1]
            F81(M)  =  F30 + G2(M) + G3(M)+ A2(M)*S(M) + (OCEAN)*F78(M)*F72(M)+ G5(M) + &
                       CLD(M)*(F30P + A3(M)*S(M)) + ( (OCEAN)*F78(M)*F76(M)- F46(M)-F41(M) )  ! [W m^-2]
            F89(M)  =  F82/F83/DT - ( F78(M)*F73(M) + F38 ) / F83! [seg^-1]
            !F89(M)  =  F82/F83/DT - (F78(M)*F73(M)*TAP2TMP+F38)/F83!
            F90(M)  = - TMPP(M)*F82/F83/DT - F81(M)/F83     ! [K seg^-1]
        enddo
        !===============================================================================
        if ( ICASE==1 ) then
        M = 1
        do j = 1 , NLAT
            do i = 1 , NLON
                AUX7(i,j,N) = TMPP(M) +TM0 -273.15
!                F90B(I,J,N) =  F90(M)
!                F89B(I,J,N) =  F89(M)       ! Convert vector into matrix
!                F81B(I,J,N) =  F81(M)
!                F78B(I,J,N) =  F78(M)       ! "    "
!                F77B(I,J,N) =  F77(M)
!                F76B(I,J,N) =  F76(M)
!                F75B(I,J,N) =  F75(M)
!                F73B(I,J,N) =  F73(M)       ! Convert vector into matrix
!                F72B(I,J,N) =  F72(M)       ! "    "
!                F71B(I,J,N) =  F71(M)       ! Convert vector into matrix
!                !F48A(I,J,N) =  F48(M)       ! "    "
                m=m+1
            enddo
        enddo
        else
        m=1
        do j = 1 , NLAT
            do i = 1 , NLON
                AUX8(i,j,N) = TMDN(M)
                m=m+1
            enddo
        enddo
        endif
        M = 1
        do j = 1 , NLAT
            do i = 1 , NLON
                AUX11(i,j,N)=  CLD(M)
                AUX12(i,j,N)=  ALBS(M)
                AUX1 (i,j,N)=  A2(M)*S(M)
                AUX2 (i,j,N)=  A3(M)*S(M)
                aux5 (i,j,N)= LANDMASK(M)
                !aux6 (i,j,k)=  S(M)
                H500N(i,j,N)=  H500(M)
                H250N(i,j,N)=  H250(M)
                A1A(i,j,N)  =  A1(M)
                G5A(i,j,N)  =  G5(M)
                if ( ICASE==2 ) then
                    AUX4(i,j,N)=TSPP(M)

                endif
                F90A(I,J) =  F90(M)       ! "    "
                F89A(I,J) =  F89(M)       ! Convert vector into matrix
                F81A(I,J) =  F81(M)       ! Convert vector into matrix
                F78A(I,J,N) =  F78(M)       ! "    "
                F77A(I,J,N) =  F77(M)
                F76A(I,J,N) =  F76(M)
                F75A(I,J,N) =  F75(M)
                F73A(I,J,N) =  F73(M)       ! Convert vector into matrix
                F72A(I,J,N) =  F72(M)       ! "    "
                F71A(I,J,N) =  F71(M)       ! Convert vector into matrix
                F48A(I,J,N) =  F48(M)       ! "    "
                FP          = -LAT(J)*pi/180
                COSFP(I,J)  =  SQRT(1. - SIN( FP )**2 )
                SINFP(I,J)  =  SIN( FP )
                F84(I,J)    = -UM(M)*(1-ADV)
                F85(I,J)    = -VM(M)*(1-ADV)-KAUST*SINFP(I,J)/COSFP(I,J)/R
                M           =  M + 1
            enddo
        enddo
! === CALL LAMBDA ==============================================================
        call LAMBDA( TPM , F90A , F89A , NLON , NLAT )
! === CALL RELAX ===============================================================
        call RELAX ( ICASE, TPM, F90A, F89A, F84, F85, COSFP, NLON, NLAT )
        if ( ICASE == 1 ) then
        !write(*,*) 'ICASE = ' , ICASE
        do j = 1 , NLAT
            do i = 1 , NLON
                TPMN(i,j,N)= TPM(i,j)
            enddo
        enddo
        else
        !write(*,*) 'ICASE = ' , ICASE
        do j = 1 , NLAT
            do i = 1 , NLON
                H250N(i,j,N)= H250N(i,j,N) + A*( TPM(i,j) - TPMN(i,j,N) )
                TPMAB(i,j,N) = TPM(i,j)
            enddo
        enddo
        endif
        call HGT(TPM+TM0,H250N(:,:,N),H500A(:,:,N),GAMMA,G,RD,P0,H0,P5,H5)
        !write(*,*) F72(9024) , F73(9024) , F76(9024) , F77(9024) , F75(9024)
        !write(*,*) '         ', '1' , '      ', TMDN(9024), TSPN(9024), TPMN(86,46,N), G5DN(9024)-G5N(9024)
        M = 1
        do J = 1 , NLAT
            do I = 1 , NLON
                if ( ICASE == 1 ) then
                T500(i,j,N)  = (TPM(i,j)+TM0)-GAMMA*(H500A(i,j,N)-H0/2.)
                G5DN(M)  = G5N(M) ! +  CLAPP*(BP(N)*(TPM(i,j) - TMN(i,j)) +     &
                !CP(M)*(TPM(i+1,j) - TPM(i-1,j) - TMN(i+1,j) + TMN(i-1,j)) +    &
                !DP(M)*(TPM(i,j+1) - TPM(i,j-1) - TMN(i,j+1) + TMN(i,j-1)) +    &
                !EPP(M)*(TPM(i,j+1) + TPM(i,j-1) + TPM(i+1,j) + TPM(i-1,j) -    &
                !4.*TPM(i,j) - TMN(i,j+1) - TMN(i,j-1) - TMN(i+1,j) - TMN(i-1,j) + &
                !4.*TMN(i,j)))
                IF( G5DN(M) < 0. ) G5DN(M) = 0.
                TMDN(M)  = TPM(i,j)!*CASEDN
                !write(*,*) 'F71   *',TPMN(i,j,N)
                !F71(M)  =  DS*(OCEANMASK(M))/DT-F36+F41(M)+F46(M)
                !F72(M)  =  (TSPP(M)*DS*(OCEANMASK(M))/DT-G2(M)-G3(M)+F48(M))/F71(M)
                !F73(M)  =  (F35-F42(M)+F46(M))/F71(M)           ! [-]
                !F75(M)  =  F49(M) /F71(M)
                !F76(M)  =  ( F41(M) + F46(M) ) / F71(M)          ! [-]
                TSDN(i,j,N)   = F72(M) + F73(M) * TMDN(M) !+ F76(M)*(TSPN(M)-ts0)*CASEDN  &
                                !+ F77(M)*TPMN(i,j,N)*CASEDN + F75(M) * D2 * ( G5DN(M)-G5N(M) ) *CASEDN

                !TSUP(i,j,N)   =  TSDN(i,j,N) - gamma*(TOPO(i,j))
                !G5DN(M)  =  G5DN(M) !+ APP(M)*( TSDN(i,j,N) - TSN(i,j,N) )
                !CLDDN(M) = CLD(M) + D2 * ( G5DN(M)-G5N(M) )*CASEDN
                !CLDDN(M) = -1.26e-2 *
                A1(M)    = QQ0(M) * ( 1. - (1.-XK(M))*CLDDN(M))*(1.-ALBS(M))
                ESDN(i,j,N)  = -A1(M) + F34 + F34P*CLD(M) + F35*TMDN(M) + &
                               F36*TSDN(i,j,N)
                ETDN(i,j,N)  = F30 + F30P*CLDDN(M) + F31*TMDN(M) + ( F32  + &
                               F32P*CLDDN(M) ) * TSN(i,j,N) + S(M) * &
                               ( CLD(M)*A3(M)+A2(M) )
                !G2DN(M)  =  G2N(M) + VAN(M)* (TSDN(i,j,N)-TSN(i,j,N)-TMDN(M)+TPMN(i,j,N))
                !G3DN(M)  =  0.981*(TSDN(i,j,N)-TSN(i,j,N)) - UN(M)*(TMDN(M) - TPMN(i,j,N)))  !! Falta
                ABSDN(i,j,N) =  A1(M) + S(M) * (CLDDN(M)*A3(M)+A2(M))
                M            = M + 1
                else
                TPMDN(i,j,N) = TPMAB(i,j,N) -  TPMN(i,j,N)
                !TSDN(i,j,N) = F72(M) + F73(M)*TPMDN(i,j,N) + F76*TSN(i,j,N) + &
                !               F77*TPMN(i,j,N) + F75*( G5(M)-G5N(M) )
                !M            = M + 1
                endif
            enddo
        enddo

!!!!!!! SUBROUTINE SNOW
        !loop=loop+1
            eps1=0.10
            m=1
            DO J = 1 , NLAT
                DO I = 1 , NLON
                    if ( loop==0 ) TSDN2(i,j,LOOP+1)=TSDN(i,j,N)
                    if ( TSDN(i,j,N) + TS0 -273.15 < 0 ) then
                        ALBSNEW(i,j,LOOP+1)=0.8
                    ELSE
                        ALBSNEW(i,j,LOOP+1)=AUX12(i,j,N)
                    ENDIF
                    ALBS(M)=ALBSNEW(i,j,LOOP+1)
                    m = m + 1
                    if ( loop > 1 ) then
                        if ( TSDN(i,j,N) - TSDN2(i,j,LOOP+1) < eps1) flag=1
                    endif
                ENDDO
            ENDDO
            write(*,*) 'LOOP' , loop , 'Subrutina SNOW'
            loop=loop+1
            if (loop<5 .or. flag==0) go to 200
call MAKEDAT3D(OF//'TSNLOOP.NC' , 'TS' ,  TSDN + TS0 - 273.15 , &
               'Monthly Long-Term Mean of Surface temperature' , &
               'degC' , 'Surface temperature' , 'MTCG output' ,  &
               'Surface' , LON , LAT , real( (/(i, i=1,LOOP-1) /) ) , NLON , NLAT , LOOP-1 )
call MAKEDAT3D(OF//'ALBLOOP.NC' , 'ALBS' ,  ALBSNEW , &
               'Monthly Long-Term Mean of Surface albedo' , &
               'degC' , 'Surface albedo' , 'MTCG output' ,  &
               'Surface' , LON , LAT , real( (/(i, i=1,LOOP-1) /) ) , NLON , NLAT , LOOP-1 )






!!!!!!! END SUBRUTINE SNOW
        if ( ICASE == 1 ) then
            do j = 1 , NLAT
                do i = 1 , NLON
                    TSN(i,j,N)  = TSDN(i,j,N)
                enddo
            enddo
        endif

    enddo
MNTH(N)=CURRENTDATE
CURRENTMES=CURRENTMES+1
CURRENTDATE=CURRENTYEAR*100+CURRENTMES
if (mod(CURRENTMES,13)==0) then
CURRENTYEAR=CURRENTYEAR+1
CURRENTMES=1
CURRENTDATE=CURRENTYEAR*100+CURRENTMES
endif

enddo
! Temperatures
call MAKEDAT3D(OF//'TPMN.NC' , 'TPMN' ,  TPMN + TM0 - 273.15 , &
               'Monthly Long-Term Mean of Mid-troposphere temperature' , &
               'degC' , 'Mid-troposphere temperature' , 'MTCG output' ,  &
               'Troposphere' , LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
call MAKEDAT3D(OF//'TPMAB.NC' , 'TPMAB' ,  TPMAB+ TM0 - 273.15 , &
               'Monthly Mean of Mid-troposphere temperature' , &
               'degC' , 'Mid-troposphere temperature' , 'MTCG output' ,  &
               'Troposphere' , LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
call MAKEDAT3D(OF//'TPMDN.NC' , 'TPMDN' , TPMDN , &
               'Monthly Anomaly of Mid-troposphere Temperature',&
               'degC' , 'Mid-troposphere Temperature Anomaly' , 'MTCG output', &
               'Troposphere' , LON , LAT , real(MNTH) , NLON , NLAT , NTIME )


call MAKEDAT3D(OF//'TMPP.NC' , 'TMPP' , AUX7 , &
               'Monthly Long-Term Mean of input TM' , &
               'K' , 'Monthly Long-Term Mean of input TM' , 'MTCG input' , 'Atmosphere' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
call MAKEDAT3D(OF//'TMDN.NC' , 'TMDN' , AUX8 , &
               'Monthly Long-Term Mean of input TM anomaly' , &
               'K' , 'Monthly Long-Term Mean of input TM anomaly' , 'MTCG input' , 'Atmosphere' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )


!call MAKEDAT3D(OF//'T500.NC' , 'T500' , T500 - 273.15 , &
!               'Monthly Long-Term Mean of 500mb temperature' , &
!               'degC' , '500mb temperature' , 'MTCG output' , '500mb' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
call MAKEDAT3D(OF//'ALBS.NC' , 'ALBS' , AUX12 , &
               'Monthly Long-Term Mean of Surface Albedo' , &
               'n/a' , 'Surface Albedo' , 'MTCG output' , 'Surface' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
call MAKEDAT3D(OF//'ALBSNEW.NC' , 'ALBSNEW' , ALBSNEW , &
               'Monthly Long-Term Mean of Surface Albedo' , &
               'n/a' , 'Surface Albedo' , 'MTCG output' , 'Surface' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )

!
!call MAKEDAT3D(OF//'A2.NC' , 'A2' , AUX1 , &
!               'Monthly Long-Term Mean of Net Solar Radiation in Atmosphere' , &
!               'W m^-2' , 'Net Solar Radiation in Atmosphere' , 'MTCG output' , 'Atmosphere' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
!call MAKEDAT3D(OF//'A3.NC' , 'A3' , AUX2 , &
!               'Monthly Long-Term Mean of Net Solar Radiation in Clouds' , &
!               'W m^-2' , 'Net Solar Radiation in Clouds' , 'MTCG output' , 'Clouds' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
!
!
call MAKEDAT3D(OF//'TSN.NC' , 'TSN' , TSDN + TS0 - 273.15 , &
               'Monthly Long-Term Mean of Surface temperature' , &
               'degC' , 'Surface temperature' , 'MTCG output' , 'Suface' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
call MAKEDAT3D(OF//'TSDNIN.NC' , 'TSDN' , AUX3 , &
               'Monthly Long-Term Mean of Surface temperature' , &
               'degC' , 'Surface temperature' , 'MTCG input' , 'Suface' , &
               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
!call MAKEDAT3D(OF//'TSAB.NC' , 'TSAB' , AUX4 + TS0 , &
!'Monthly Long-Term Mean of Surface temperature' , &
!'degC' , 'Surface temperature' , 'MTCG input' , 'Suface' , &
!LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
call MAKEDAT3D(OF//'LAND.NC' , 'land' , AUX5 , &
'Monthly Long-Term Mean of Land Mask' , &
'N/A' , 'Land Mask' , 'MTCG input' , 'Suface' , &
LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
!call MAKEDAT3D(OF//'G5.NC' , 'G5' , G5A , &
!               'Monthly Long-Term Mean of Latent Heat of condensation' , &
!               'W m^-2' , 'Latent Heat of condensation' , 'MTCG input' , 'Troposphere' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!
!call MAKEDAT3D(OF//'A1.NC' , 'A1' , A1A , &
!               'Monthly Long-Term Mean of Net Surface Radiation' , &
!               'W m^-2' , 'Net Surface Radiation' , 'MTCG output' , 'Suface' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!call MAKEDAT3D(OF//'CLD.NC' , 'CLD' , AUX11 , &
!               'Monthly Long-Term Mean of Cloud cover' , &
!               'W m^-2' , 'Cloud cover' , 'MTCG output' , 'Troposphere' , &
!               LON , LAT , real(MNTH) , NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/P500.NC' , 'P500' , P500 , 'Pressure' , 'Pa' , &
!               'Monthly Long-Term Mean of Tropospheric Pressure' , &
!               'MTCG calculated' , 'Troposphere' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/H250.NC' , 'H250' , H250N , 'Geopotential Height' , &
!               'meters' ,'Monthly Long-Term Mean of Geopotential Height' , &
!               'MTCG input' , '250mb' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/H500B.NC' , 'H500B' , H500A , 'Geopotential Height' , &
!               'meters' ,'Monthly Long-Term Mean of Geopotential Height' , &
!               'MTCG calculated' , '500mb' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/ESDN.NC' , 'ESDN' , ESDN , 'Net surface radiation anomaly' , &
!               'W m^-2' ,'Monthly Long-Term Mean of Net surface radiation anomaly' , &
!               'MTCG calculated' , 'surface' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/ETDN.NC' , 'ETDN' , ETDN , 'Net atmospheric radiation anomaly' , &
!               'W m^-2' ,'Monthly Long-Term Mean of Net atmospheric radiation anomaly' , &
!               'MTCG calculated' , 'surface' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )

!call MAKEDAT3D( 'output/F71.NC' , 'F71' , F71A , 'F71' , &
!               'n/a' ,'F71' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!
!call MAKEDAT3D( 'output/F72.NC' , 'F72' , F72A , 'F72' , &
!               'n/a' ,'F72' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F73.NC' , 'F73' , F73A , 'F73' , &
!               'n/a' ,'F73' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F75.NC' , 'F75' , F75A , 'F75' , &
!               'n/a' ,'F75' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F48.NC' , 'F48' , F48A , 'F48' , &
!               'n/a' ,'F48' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!
!call MAKEDAT3D( 'output/F76.NC' , 'F76' , F76A , 'F76' , &
!               'n/a' ,'F76' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!
!call MAKEDAT3D( 'output/F77.NC' , 'F77' , F77A , 'F77' , &
!               'n/a' ,'F77' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F78.NC' , 'F78' , F78A , 'F78' , &
!               'n/a' ,'F78' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F81.NC' , 'F81' , F81A , 'F81' , &
!               'n/a' ,'F81' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F89.NC' , 'F89' , F89A , 'F89' , &
!               'n/a' ,'F89' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F90.NC' , 'F90' , F90A , 'F90' , &
!               'n/a' ,'F90' , &
!               'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )


!call MAKEDAT3D( 'output/F71N.NC' , 'F71' , F71B , 'F71' , &
!'n/a' ,'F71' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F72N.NC' , 'F72' , F72B , 'F72' , &
!'n/a' ,'F72' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F73N.NC' , 'F73' , F73B , 'F73' , &
!'n/a' ,'F73' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F75N.NC' , 'F75' , F75B , 'F75' , &
!'n/a' ,'F75' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!!call MAKEDAT3D( 'output/F48N.NC' , 'F48' , F48B , 'F48' , &
!!'n/a' ,'F48' , &
!!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F76N.NC' , 'F76' , F76B , 'F76' , &
!'n/a' ,'F76' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F77N.NC' , 'F77' , F77B , 'F77' , &
!'n/a' ,'F77' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F78N.NC' , 'F78' , F78B , 'F78' , &
!'n/a' ,'F78' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F81N.NC' , 'F81' , F81B , 'F81' , &
!'n/a' ,'F81 normal' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F89N.NC' , 'F89' , F89B , 'F89' , &
!'n/a' ,'F89 normal' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/F90N.NC' , 'F90' , F90B , 'F90' , &
!'n/a' ,'F90 normal' , &
!'MTCG calculated' , 'n/a' , LON , LAT , real(MNTH) , &
!NLON , NLAT , NTIME )
!
!call MAKEDAT3D( 'output/S.NC' , 'S' , AUX1 , 'S' , &
!               'W m^2' ,'S' , &
!               'MTCG calculated' , 'TOA' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )
!call MAKEDAT3D( 'output/TSUP.NC' , 'Ts' , TSUP , 'Surface temperature (with toporgaphy)' , &
!               'K' ,'Monthly Long-Term Mean of Surface temperature (with toporgaphy)' , &
!               'MTCG calculated' , 'surface' , LON , LAT , real(MNTH) , &
!               NLON , NLAT , NTIME )


! === Close input files, write output files and end program. ===================
do i = 31 , 52
      close(i)
enddo
call cpu_time( t2 )
write(*,FORMAT3)
write(*,FORMAT2) '_________________'
write(*,FORMAT2) 'END OF EXPERIMENT'
write(*,*) 'in',t2-t1,'seconds'
call exit(0)
! ==============================================================================
! ==========================  end MAIN PROGRAM  ================================
! ==============================================================================
contains
!===============================================================================
! === MILANK subroutine ========================================================
! ===== Soubrutine to compute the latitudinal solar irradiance at the top of ===
! ===== the atmosphere (TOA) using milankovitch formula on a given day. ========
! ==============================================================================
!    IT = s0 (a0/r)^2 [ w1 sin(phi) sin(delta) + sin(w1) cos(phi) cos(delta) ]
!===============================================================================
subroutine MILANK( LATS, NLAT, S0, EPSI, EXEN, LS , IT)
implicit none
integer :: i, j
real :: a0r, delta, w1, cosz, yr
real, parameter  :: pi  = 3.141592
real, parameter  :: lsp = 90.*pi/180.
real, intent(in) :: S0, EPSI, EXEN
integer, intent(in) :: LS
integer, intent(in) :: NLAT
real, dimension(NLAT) :: PHI
real, dimension(NLAT), intent(in)  :: LATS
real, dimension(NLAT), intent(out) :: IT
yr    = ls
a0r   = (1+EXEN*COS(yr*pi/180.))/(1-EXEN**2)
delta = ASIN(SIN(EPSI)*SIN((yr-LSP)*pi/180.))
do j = 1 , NLAT
      PHI(J) = LATS(J)*pi/180.
      w1    = ACOS( - TAN(PHI(j)) * TAN(delta) )
      cosz  = (w1*SIN(PHI(j))*SIN(delta) + SIN(w1)*COS(PHI(j))*COS(delta))
      IT(j) = (s0/pi) * cosz*(a0r)**2
      ! at poles
      if ( ABS(PHI(j)) > pi/2.-ABS(delta) ) then
      IT(j) = 0.0                                             ! noche polar
          if ( lat(j)*delta>0 ) then
              IT(j) = s0*SIN(lat(j))*SIN(delta)*(a0r)**2      ! dia polar
          endif
      endif
enddo
end subroutine MILANK
!===============================================================================
! === RADIAT subroutine ========================================================
! ===== subroutine to display the coeficients from Long-wave radiation =========
! ===== balance in MTC. ( F30, F30p, F32, F34, F34p, F35, F36, F38 ).  =========
!===============================================================================
subroutine RADIATIVE(T0,TA0,TS0,F30,F30p,F31,F32,F32p,F33,F34,F34p,F35,F36,F38)
implicit none
real, intent(in)    :: t0 , ta0 , ts0
real, intent(out)   :: F30, F30p, F31,F32, F32p, F33,F38, F34, F34p, F35, F36
character(len = 17) :: FORMAT6 ="(X, A, 39X, A23)"
character(len = 17) :: FORMAT7 ="(X, A, 37X, A25)"
character(len = 42) :: FORMAT8 ="(X, A, 37X, A, 2X, A4, A, 2X F8.3, 6X, A)"
! ET
F30  =  F(T0) - SIGMA*T0**4 + F(TA0) - SIGMA*TA0**4 - F(TS0) + SIGMA*TS0**4
F30p = -2.*F(TC0) + F(TS0)
F33 =  DF(TA0)-4.*SIGMA*TA0**3                            ! [W M^-2]
F31 =  DF(T0) - 4.*SIGMA*T0**3 + F33
F32  =  4.*SIGMA*TS0**3 - DF(TS0)
F32P =  DF(TS0)
F38  =  F33*TAP2TMP + (F31-F33)*TP2TMP! [W M^-2 K^-1]
! ES
F34  =  SIGMA*TA0**4 - F(TA0)-SIGMA*TS0**4             ! [W M^-2]
F34p =  F(TC0)                                         ! [W M^-2]
F35  = -F33 !* TAP2TMP                                 ! [W M^-2 K^-1]
F36  = -4.*SIGMA*TS0**3                                ! [W M^-2 K^-1]
!    F30 = -114.059723
!    F30p=  -33.0982513
!    F31 =   -4.84247971
!    F32 =    3.16219831
!    F32p=    2.25595331
!
!    F33 =   -3.16219878
!    F34 = -139.039246
!    F34p=   86.0687485
!    F35 =    3.16219878
!    F36 =   -5.41815186
write(*,FORMAT6) '|','Radiative coeficients:|'
write(*,FORMAT7) '|','************************|'
write(*,FORMAT8) '|', '*','F30 ','=',F30, '|'
write(*,FORMAT8) '|','*','F30p' , '=' , F30p, '|'
write(*,FORMAT8) '|','*','F32 ' , '=' , 4*SIGMA*TS0**3 + 0.5*DF(TS0), '|'
write(*,FORMAT8) '|','*','F33 ' , '=' , F33, '|'
write(*,FORMAT8) '|','*','F34 ' , '=' , F34, '|'
write(*,FORMAT8) '|','*','F34p' , '=' , F34p, '|'
write(*,FORMAT8) '|','*','F35 ' , '=' , F35, '|'
write(*,FORMAT8) '|','*','F36 ' , '=' , F36, '|'
write(*,FORMAT7) '|','************************|'
!write(*,*) F30, F30p,F31, F32, F32p,F33, F34, F34p, F35, F36,F38
end subroutine RADIATIVE
!===============================================================================
! === LW ABSORPTION FUNCTIONS ==================================================
!===============================================================================
real function F(T)
implicit none
integer :: YEAR
real, intent(in) :: T
real :: ABS1=0.0, ABS2=0.067, ABS3=0.802, ABS4=0.931, ABS5=0.310
F  = (1-ABS1)*(DFTL(T,12e-6)-DFTL(T,08e-6)) + &
     (1-ABS2)*( FTL(T,13e-6)-FTL(T,12e-6) ) + & ! CO2
     (1-ABS3)*( FTL(T,14e-6)-FTL(T,13e-6) ) + & ! CO2
     (1-ABS4)*( FTL(T,17e-6)-FTL(T,16e-6) ) + &
     (1-ABS5)*( FTL(T,18e-6)-FTL(T,17e-6) )
end function F
!
real function DF(T)
implicit none
real, intent(in) :: T
real :: ABS1=0.0, ABS2=0.067, ABS3=0.802, ABS4=0.931, ABS5=0.310
DF = (1-ABS1)*(DFTL(T,12e-6)-DFTL(T,08e-6)) + &
     (1-ABS2)*(DFTL(T,13e-6)-DFTL(T,12e-6)) + &
     (1-ABS3)*(DFTL(T,14e-6)-DFTL(T,13e-6)) + &
     (1-ABS4)*(DFTL(T,17e-6)-DFTL(T,16e-6)) + &
     (1-ABS5)*(DFTL(T,18e-6)-DFTL(T,17e-6))
end function DF
!
real function FTL(T,LB)
implicit none
real, intent(in) :: T, LB
FTL  = C1 * EXP(-C2/(LB*T)) * ( T/(C2*LB**3) + 3.*(T/(C2*LB))**2 + &
       6.*(T**3/(LB*C2**3) + 6.*(T/C2)**4))
end function FTL
!
real function DFTL(T,LB)
implicit none
real, intent(in) :: T, LB
DFTL = C1 * EXP(-C2/(LB*T)) * ( 1./(T*LB**4) + 4./(C2*LB**3) + &
       12.*T/(LB*C2)**2 + 24.*(T**2)/(LB*C2**3) + (T**3/C2**4))
end function DFTL
!
real function XKADEM(LATITUDE)
real, intent(in) :: LATITUDE
real, dimension(20) :: KSA=(/ 0.35 , 0.35 , 0.34 , 0.34 , 0.33 , 0.33 , 0.32 , &
                              0.32 , 0.32 , 0.33 , 0.34 , 0.36 , 0.38 , 0.40 , &
                              0.45 , 0.50 , 0.55 , 0.55 , 0.55 , 0.55 /)
XKADEM = KSA(int(abs(LATITUDE)/5)+1)
end function XKADEM

real function A2ADEM(LATITUDE,SEASON)
integer :: SEASON
real, intent(in) :: LATITUDE
real, dimension(20) :: A22

if ( SEASON == 1 ) then
A22=(/ 0.134 , 0.134 , 0.123 , 0.123 , 0.122 , 0.122 , 0.105 , &
       0.105 , 0.100 , 0.100 , 0.095 , 0.095 , 0.108 , 0.108 , &
       0.143 , 0.143 , 0.143 , 0.143 , 0.143 , 0.143 /)
endif

if ( SEASON == 2) then
A22=(/ 0.143 , 0.143 , 0.145 , 0.145 , 0.145 , 0.145 , 0.142 , &
       0.142 , 0.135 , 0.135 , 0.132 , 0.132 , 0.147 , 0.147 , &
       0.169 , 0.169 , 0.164 , 0.164 , 0.164 , 0.164 /)
endif

if ( SEASON == 3) then
A22=(/ 0.142 , 0.142 , 0.141 , 0.141 , 0.140 , 0.140 , 0.123 , &
       0.123 , 0.116 , 0.116 , 0.110 , 0.110 , 0.120 , 0.120 , &
       0.149 , 0.149 , 0. , 0. , 0. , 0. /)
endif

if ( SEASON == 4) then
A22=(/ 0.137 , 0.137 , 0.134 , 0.134 , 0.126 , 0.126 , 0.115 , &
       0.115 , 0.104 , 0.104 , 0.1 , 0.1 , 0.114 , 0.114 , &
       0.167 , 0.167 , 0. , 0. , 0. , 0. /)
endif

A2ADEM = A22(int(abs(LATITUDE)/5)+1)
end function A2ADEM

real function A3ADEM(LATITUDE,SEASON)
integer :: SEASON
real, intent(in) :: LATITUDE
real, dimension(20) :: A33
if (SEASON==1) then
    A33=(/ 0.034 , 0.034 , 0.030 , 0.030 , 0.030 , 0.030 , 0.032 , &
        0.032 , 0.033 , 0.033 , 0.035 , 0.035 , 0.042 , 0.042 , &
        0.038 , 0.038 , 0.033 , 0.033 , 0.033 , 0.033 /)
endif

if (SEASON==2) then
    A33=(/ 0.034 , 0.034 , 0.038 , 0.038 , 0.036 , 0.036 , 0.033 , &
            0.033 , 0.033 , 0.033 , 0.035 , 0.035 , 0.037 , 0.037 , &
            0.033 , 0.033 , 0.034 , 0.034 , 0.034 , 0.034 /)
endif

if (SEASON==3) then
    A33=(/ 0.031 , 0.031 , 0.032 , 0.032 , 0.034 , 0.034 , 0.035 , &
            0.035 , 0.042 , 0.042 , 0.037 , 0.037 , 0.050 , 0.050 , &
            0.030 , 0.030 , 0. , 0. , 0. , 0. /)
endif

if (SEASON==4) then
    A33=(/ 0.033 , 0.033 , 0.032 , 0.032 , 0.039 , 0.039 , 0.036 , &
            0.036 , 0.044 , 0.049 , 0.049 , 0.039 , 0.039 , 0. , &
            0. , 0. , 0. , 0. , 0. , 0. /)
endif
A3ADEM = A33(int(abs(LATITUDE)/5)+1)
end function A3ADEM


!===============================================================================
! === INSO_TOA subroutine ======================================================
! ===== Soubrutine to compute the monthly mean of solar irradiance at the top ==
! ===== of the atmosphere (TOA) using milankovitch formula on a given month. ===
!===============================================================================
!******* REVISAR NDAY=30 *****************|
subroutine INSO_TOA ( MONTH , S )
implicit none
integer :: MONTH , M ,  NINI , NFIN
integer,dimension(12) :: LASTDAY, NDAY
!real, dimension(NLAT,NDAY) :: INSOL
!real, dimension(NDAY,NLAT) :: INSOC
real, dimension(:,:), allocatable :: INSOL
real, dimension(:,:), allocatable :: INSOC
real, dimension(NLON,NLAT) :: VAROUT
real, dimension(NLON,NLAT) :: VAR1 , VARSUM
real, dimension(NLON*NLAT),intent(out) :: S
LASTDAY =  (/ 31,59,90,120,151,181,212,243,273,304,334,365 /)
NDAY = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

NINI = LASTDAY(MONTH)-NDAY(MONTH)+1
NFIN = LASTDAY(MONTH)
allocate(INSOL(NLAT,NDAY(MONTH)))
allocate(INSOC(NDAY(MONTH),NLAT))
VARSUM = 0.
!M = 1
do LS = NINI , NFIN
      call MILANK( -LAT, NLAT, S0, EPS, EX, LS , INSOL(:,m))
      do j = 1 , NLAT
            do i = 1 , NLON
                  VARSUM(i,j) = INSOL(j,LS) + VARSUM(i,j)
            enddo
      enddo
!      if ( MOD(m,NDAY) == 0 ) then

!      endif
!      M = M + 1
enddo
      do j = 1 , NLAT
            do i = 1 , NLON
                  INSO(i,j,k) = VARSUM(i,j) / (NFIN-NINI)
            enddo
      enddo
      call POLE(NLON,NLAT,INSO(:,:,N))
m = 1
do j = 1 , NLAT
      do i = 1 , NLON
            S(M) = INSO(i,j,k)
            M = M + 1
      enddo
enddo
end subroutine INSO_TOA
!===============================================================================
! === POLE subroutine ==========================================================
! ===== subroutine to average data in poles from the nearest data. =============
!===============================================================================
subroutine POLE( XLON , XLAT, VAR )
implicit none
integer :: XLON , XLAT
real :: TSUM1 , TSUM2
real, dimension(XLON,XLAT) :: VAR
TSUM1 = 0.
TSUM2 = 0.
do I = 1 , XLON
      TSUM1 = TSUM1 + VAR(I,2)
      TSUM2 = TSUM2 + VAR(I,XLAT-1)
enddo
do I = 1 , XLON
      VAR(I,1) = TSUM1/XLON
      VAR(I,XLAT)= TSUM2/XLON
enddo
end subroutine POLE
!===============================================================================
! === LAMBDA subroutine ========================================================
!===============================================================================
subroutine  LAMBDA  ( TPM , F90A , F89, PLON, PLAT )
implicit none
integer,intent(in) :: PLON, PLAT
real, dimension(PLON,PLAT),intent(in)  :: F90A, F89
real, dimension(PLON,PLAT),intent(out) :: TPM
do J = 1 , PLAT
      do I = 1 , PLON
            TPM(I,J) = - F90A(I,J)/F89(I,J)
      enddo
enddo
call POLE(PLON,PLAT,TPM)
end subroutine LAMBDA
!===============================================================================
! === RELAX subroutine =========================================================
!===============================================================================
subroutine RELAX( ICASE, TPM, F90, F89, F84, F85, COSF, PLON, PLAT)
implicit none
real, dimension(PLON,PLAT) :: TPM, F90, F89, F84, F85, COSF
!real, dimension(:,:,:), allocatable :: RELAXTPM
real :: ML , MF , DML , DMF , TPMN , TPMNB , TPMNR , TPMNL , TPMNT
real :: E1, E2, E3, E4, E5, E6, RES, ERRO, SS, AR
real :: TSUM1, TSUM2, S, W, DELTAF, DELTAL, t1, t2
integer :: PLON, PLAT, ICASE, RCYCLE, RLIMIT
DELTAF  = pi/PLAT
DELTAL  = 2.*pi/PLON
AR      = 0.5
ERRO    = 1.0E-2
RLIMIT  = 9999
!allocate(RELAXTPM(PLON,PLAT,int(RLIMIT)))
!if ( N == 1 ) call MAKEDAT2D('TPMIN.NC','TPM',TPM,'K',LON,LAT,PLON,PLAT)
call cpu_time ( t1 )
do RCYCLE = 1 , RLIMIT
      W = 0.
      do J = 2 , PLAT-1
            do I = 1 , PLON
                  if ( I == 1 ) then
                        TPMNL    =  TPM(PLON,J)
                  else
                        TPMNL    =  TPM(I-1,J)
                  endif
                  if( I == PLON ) then
                        TPMNR    =  TPM(1,J)
                  else
                        TPMNR    =  TPM(I+1,J)
                  endif
                  TPMNT    =  TPM(I,J+1)
                  TPMNB    =  TPM(I,J-1)
                  TPMN     =  TPM(I,J)
                  ML       =  1./DELTAL/COSF(I,J)
                  MF       =  1./DELTAF
                  DML      =  R/ML/2.
                  DMF      =  R/MF/2.
                  E1       =  KAUST + F84(i,j)*DML
                  E2       =  (MF/ML)**2*(KAUST + F85(i,j)*DMF)
                  E3       = -E1 + 2.*KAUST
                  E4       = -E2 + 2.*(MF/ML)**2*KAUST
                  E5       =  2.*KAUST*(1. + (MF/ML)**2) - 4.*F89(I,J)*DML**2
                  E6       =  4.*DML**2*F90(i,j)
                  RES      =  E1*TPMNR+E2*TPMNT+E3*TPMNL+E4*TPMNB-E5*TPMN+E6
                  S        =  RES*AR/E5
                  if ( ABS( S ) > ERRO )  W = 1.
                  TPM(I,J) =  TPM(I,J) + S
            enddo
      enddo
      call POLE( PLON , PLAT, TPM )
!!     For print relax process
!      do j = 1 , PLAT
!            do i = 1 , PLON
!                  RELAXTPM(i,j,RCYCLE) = TPM(i,j)
!            enddo
!      enddo
      if ( W == 0. )  exit
      if ( RCYCLE == RLIMIT ) write(*,"(X,A,3X,A16,42X,A)") '|','Did not converge','|'
enddo
call cpu_time ( t2 )
write(*,'(X,A,3X , F5.0, A19, F5.3,X, A17, I2,9X,A)')  '|',real(RCYCLE) , &
         "Cycles to relax in ", ( T2 - T1 ) , "seconds for case " , ICASE, '|'
!if (N==7) call MAKEDAT3D('RELAX.NC','TPM',RELAXTPM,'none','none','none','none','DegC',LON,LAT,1.,PLON,PLAT,REAL(RCYCLE))
!deallocate(RELAXTPM)
end subroutine RELAX
!
SUBROUTINE PRESS(HGT,TEMP,PRES,P1,GAMMA,H5)
!===============================================================================
! === HGTS subroutine ==========================================================
!====== Computes heights and pressures at any level from observed ==============
!====== heights (HGT) and temperatures (TEMP) at a given level of an isobaric ==
!====== surface. TT and PRESS are temperature and pressure at z=h0. ============
!===============================================================================
!integer :: NLONS,NLAT
real, intent(in) :: P1 , GAMMA , h5
real, dimension(NLON,NLAT), intent(out) :: PRES
real, dimension(NLON,NLAT), intent(in)  :: TEMP , HGT
real :: TT
do j=1,NLAT
    do i=1,NLON
        TT = TEMP(i,j) - GAMMA*( h5-HGT(i,j) )
        PRES(i,j) = P1 * ( TT / TEMP(i,j) )**(GRB)
    enddo
enddo
end subroutine PRESS

SUBROUTINE HGT(TEMP,HOBS,HEIGHT,GAMMA,G,RD,P0,H0,P5,H5)
real, dimension(NLON,NLAT), intent(out) :: HEIGHT
real, dimension(NLON,NLAT), intent(in)  :: TEMP , HOBS
real :: GAMMA , h5 , TT, PP, h0,Hup,Hlow,Tup,Tlow,rholow,Plow,P0,Rd,p5,g
do j = 1 , NLAT
    do i = 1 , NLON
        TT     = TEMP(i,j) - gamma*H0/2.
        Hup    = HOBS(i,j)
        Tup    = TT - gamma*(Hup-h0)
        Tlow   = TT + gamma*(h0-h5)
        PP     = P0*(TT/Tup)**GRB
        Plow   = PP*(Tlow/TT)**GRB
        rholow = Plow/Rd/tlow
        HEIGHT(i,j) = H5 + (Plow-P5)/rholow/g
    enddo
enddo

end subroutine HGT
!===============================================================================
! === CLIM_FILES subroutine ====================================================
!===============================================================================
subroutine CLIM_FILES
implicit none
character (len = 14), parameter :: rutaClim='input/monthly/'
write(*,FORMAT2) 'Loading climatology'
write(*,FORMAT6) 'Binary field for oceans and continents                   '
open(unit=30,FILE=rutaClim//'GRID.DATA',status='OLD')
open(unit=31,FILE=rutaClim//'LAND.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Clear Sky Downward Shortwave Radiation at surface        '
open(unit=33,FILE=rutaClim//'Qq0.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Shortwave Radiation at surface                           '
open(unit=53,FILE=rutaClim//'Qq.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Albedo Surface                                           '
open(unit=36,FILE=rutaClim//'ALBS.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Sensible Heat Net Flux at surface                        '
open(unit=37,FILE=rutaClim//'G2.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Latent Heat Net Flux at surface                          '
open(unit=38,FILE=rutaClim//'G3.DATA',status='UNKNOWN',	&
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Surface Temperature                                      '
open(unit=43,FILE=rutaClim//'TS.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Temperature at 250mb (top of the layer model)            '
open(unit=44,FILE=rutaClim//'T250.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Surface Temperature Anomaly (NCEP 1948-2015)                  '
open(unit=45,FILE=rutaClim//'TSDN.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
!write(*,FORMAT6) 'Mid-Troposphere Temperature Anomaly (NCEP 1948-2015)          '
!open(unit=46,FILE=rutaClim//'TMDN.DATA',status='UNKNOWN', &
!     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Surface Wind Speed                                       '
open(unit=47,FILE=rutaClim//'VAN.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Meridional tropospheric wind                             '
open(unit=48,FILE=rutaClim//'VM.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Zonal tropospheric wind                                  '
open(unit=49,FILE=rutaClim//'UM.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) '700mb Geopotential Height                                '
open(unit=50,FILE=rutaClim//'H700.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) '500mb Geopotential Height                                '
open(unit=51,FILE=rutaClim//'H500.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) '250mb Geopotential Height                                '
open(unit=52,FILE=rutaClim//'H250.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
write(*,FORMAT6) 'Explicit advection                                       '
open(unit=56,FILE=rutaClim//'UDTDX.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
open(unit=57,FILE=rutaClim//'VDTDY.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
open(unit=58,FILE=rutaClim//'WDTDY.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
if ( XK0 == 1 ) then
write(*,FORMAT6) '-- XK Derivated from NCEP & ERA reanalisys            '
open(unit=34,FILE=rutaClim//'XK.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) '-- XK obteined by Adem                                '
endif
if( CLDERA == 0 ) then
    open(unit=32,FILE=rutaClim//'CLD.DATA',status='UNKNOWN', &
         access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
    write(*,FORMAT6) '-- Total cloud cover from NCEP reanalysis.               '
else
    open(unit=32,FILE=rutaClim//'CLDERA.DATA',status='UNKNOWN',	&
         access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
    write(*,FORMAT6) '-- Total cloud cover from ERA reanalysis.                   '
endif
if ( MILAN == 0 ) then
write(*,FORMAT6) "-- Insolation at TOA from NCEP reanalysis                   "
open(unit=35,FILE=rutaClim//'S.DATA',status='UNKNOWN',	&
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) "-- Insolation at TOA calculated with Milankovitch's formula "
endif
if ( LH == 0 ) then
write(*,FORMAT6) "-- Latent heat of condensation from NCEP precipitation      "
open(unit=39,FILE=rutaClim//'G5.DATA',status='UNKNOWN',	&
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) "-- Latent heat of condensation from Kuo's scheme'           "
open(unit=39,FILE=rutaClim//'QC.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
endif
if ( TM01 == 0 ) then
write(*,FORMAT6) '-- Mid-troposphere temperature NCEP (500mb)                 '
open(unit=40,FILE=rutaClim//'T500.DATA',status='UNKNOWN',	&
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) '-- Mid-troposphere temperature derivated from NCEP    '
open(unit=40,FILE=rutaClim//'TM.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
endif
if ( A20 == 1 ) then
write(*,FORMAT6) '-- A2 Derivated from NCEP & ERA reanalisys                  '
open(unit=41,FILE=rutaClim//'AA2.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) '-- A2 obteined by Adem                                      '
endif
if ( A30 == 1 ) then
write(*,FORMAT6) '-- A3 Derivated from NCEP & ERA reanalisys                  '
open(unit=42,FILE=rutaClim//'BB3.DATA',status='UNKNOWN', &
     access='DIRECT',form='UNFORMATTED',recl=NLON*NLAT*8)
else
write(*,FORMAT6) '-- A3 obteined by Adem                                      '
endif
call GETDAT2D(rutaClim//'land.sfc.gauss.nc','land',land,NLON,NLAT)
write(*,FORMAT6) 'Topography                                               '
call GETDAT2D(rutaClim//'topo_gauss.nc','topo',TOPO,NLON,NLAT)
read(30,*) LON
read(30,*) LAT
write(*,FORMAT2) ' - - SUCCESS - - '
!m=1
!do j=1,NLAT
!    do i=1,NLON

        ! write(*,*) j,m
        !if ( j > 4 .and. j < 16 .and. i > 158 .and. i < 182 ) OCEANMASK(M) = 1.
        !if ( j > 80 ) OCEANMASK(M) = 1.
!!        m=m+1
!    enddo
!enddo
!write(*,*) lat
!write(*,*)
!call exit(0)
end subroutine CLIM_FILES
!===============================================================================
! === CLIM_DATA subroutine =====================================================
!===============================================================================
subroutine CLIM_DATA( N , MONTH , MONTHI )
implicit none
integer :: N , MONTH , MONTHI
!write(*,*) month
read( 31 , rec = MONTH ) LANDMASK
read( 32 , rec = MONTH ) CLD
read( 33 , rec = MONTH ) QQ0
!
if ( XK0 == 1 ) then
    read( 34 , rec = MONTH ) XK
else
    do m=1,NLON*NLAT
        XK(M)=XKADEM(int(M/NLON)+1.)
    enddo
endif
!
if ( A20 == 1 ) then
    read( 41 , rec = MONTH )  A2
    do m = 1 , NLON*NLAT
        if (A2(M) > 0.3) A2(M) = 0.3
    enddo
else
    do m=1,NLON*NLAT
        A2(M)=A2ADEM(int(M/NLON)+1.,SEA(CURRENTMES))
    enddo
endif
!
if ( A30 == 1 ) then
    read( 42 , rec = MONTH )  A3
    do m = 1 , NLON*NLAT
        if (A3(M) > 0.25) A3(M) = 0.25
    enddo
else
    do m=1,NLON*NLAT
        A3(M)=A3ADEM(int(M/NLON)+1.,SEA(CURRENTMES))
    enddo
endif
!
if ( MILAN == 1 ) then
    call INSO_TOA ( MONTH , S )
else
    read( 35 , rec = MONTH ) S
endif
!
read( 36 , rec = MONTH )  ALBS
read( 37 , rec = MONTH )  G2N
read( 38 , rec = MONTH )  G3N
read( 39 , rec = MONTH )  G5N
read( 40 , rec = MONTH )  TMN
read( 43 , rec = MONTH )  TSPN               ! TS
read( 44 , rec = MONTH )  TPN	               ! T
read( 45 , rec = REG )  TSPDN	     ! Ts anomaly previa
!write(*,*)
!write(*,*) YEARINI, MESINI
!write(*,*)
M = 1
do j = 1 , NLAT
    do i = 1 , NLON
        !OCEANMASK(M) =
        AUX3 (i,j,N)=  TSPDN(M)*LANDMASK(M)
        m=m+1
    enddo
enddo

! read( 46 , rec = MONTH )  TMDN               ! Tm anomaly previa
read( 47 , rec = MONTH )  VAN
read( 48 , rec = MONTH )  VM
read( 49 , rec = MONTH )  UM
read( 50 , rec = MONTH )  H700
read( 51 , rec = MONTH )  H500
read( 52 , rec = MONTH )  H250
read( 53 , rec = MONTH )  QQ
read( 56 , rec = MONTH )  ADV_U
read( 57 , rec = MONTH )  ADV_V
read( 58 , rec = MONTH )  CONV
TMN  = TMN  + 273.16
TSPN = TSPN + 273.16
CLD  = CLD  / 100.
ALBS = ALBS / 100.

!call PRESS(H500,TMN,P500(:,:,n),P5,GAMMA,H5)
end subroutine CLIM_DATA
!===============================================================================
! == NETCDF 3D-DATA write subroutine ===========================================
!===============================================================================
subroutine MAKEDAT3D( FILE_NAME , VAR_NAME , DATA_WROUT , LN_DATA , &
                      DATA_UNITS , VARDESC , DATASET , LVLDESC , LONS , LATS , &
                      RECS , NLONS , NLATS , NRECS )
implicit none
integer :: lat, lon, ncid
integer :: lon_varid, lat_varid, rec_varid, data_varid
integer :: lon_dimid, lat_dimid, rec_dimid
integer, intent(in) :: NRECS, NLATS, NLONS
character (len = *), intent(in) :: FILE_NAME
character (len = *), intent(in) :: VAR_NAME
character (len = *), intent(in) :: DATA_unitS
character (len = *), intent(in) :: LN_DATA
character (len = *), intent(in) :: VARDESC
character (len = *), intent(in) :: DATASET
character (len = *), intent(in) :: LVLDESC
real 	:: LATS(NLATS), LONS(NLONS), RECS(NRECS)
integer :: dimids(3)
integer :: start(3), cont(3)
character (len = *), parameter :: LAT_NAME  = "lat"
character (len = *), parameter :: LON_NAME  = "lon"
character (len = *), parameter :: rec_NAME  = "time"
character (len = *), parameter :: LAT_unitS = "degree_N"
character (len = *), parameter :: LON_unitS = "degree_E"
character (len = *), parameter :: REC_unitS = "YYYYMM"
real, dimension(:,:,:) :: DATA_WROUT
call check(nf90_create(FILE_NAME, nf90_clobber, ncid) )
call check(nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
call check(nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
call check(nf90_def_dim(ncid, rec_NAME, NRECS, rec_dimid) )
dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
call check(nf90_def_var(ncid, LAT_NAME, NF90_real, lat_dimid, lat_varid) )
call check(nf90_def_var(ncid, LON_NAME, NF90_real, lon_dimid, lon_varid) )
call check(nf90_def_var(ncid, rec_NAME, NF90_real, rec_dimid, rec_varid) )
call check(nf90_def_var(ncid, VAR_NAME, NF90_real, dimids, data_varid) )
call check( nf90_put_att(ncid, lat_varid, 'units', LAT_units) )
call check( nf90_put_att(ncid, lon_varid, 'units', LON_units) )
call check( nf90_put_att(ncid, rec_varid, 'units', REC_units) )
call check( nf90_put_att(ncid, data_varid, 'long_name', LN_DATA) )
call check( nf90_put_att(ncid, data_varid, 'units', DATA_UNITS) )
call check( nf90_put_att(ncid, data_varid, 'add_offset', 0.0 ) )
call check( nf90_put_att(ncid, data_varid, 'scale_factor', 1.0 ) )
call check( nf90_put_att(ncid, data_varid, 'precision', 2 ) )
call check( nf90_put_att(ncid, data_varid, 'least_significant_digit', 1 ) )
call check( nf90_put_att(ncid, data_varid, 'var_desc', VARDESC ) )
call check( nf90_put_att(ncid, data_varid, 'dataset', DATASET ) )
call check( nf90_put_att(ncid, data_varid, 'level_desc', LVLDESC ) )
call check( nf90_enddef(ncid) )
call check( nf90_put_var(ncid, lat_varid, LATS) )
call check( nf90_put_var(ncid, lon_varid, LONS) )
call check( nf90_put_var(ncid, rec_varid, RECS) )
cont = (/ NLONS, NLATS, NRECS /)
start = (/ 1, 1, 1 /)
call check(nf90_put_var(ncid,data_varid,DATA_WROUT,start,cont) )
call check( nf90_close(ncid) )
end subroutine MAKEDAT3D
!===============================================================================
! === NETCDF 2D-DATA write subroutine ==========================================
!===============================================================================
subroutine MAKEDAT2D( FILE_NAME , VAR_NAME , DATA_WROUT , DATA_unitS , &
                      LONS , LATS , NLONS , NLATS )
implicit none
integer :: lat, lon, ncid
integer :: lon_varid, lat_varid, data_varid
integer :: lon_dimid, lat_dimid
integer, intent(in) :: NLATS, NLONS
character (len = *), intent(in) :: FILE_NAME
character (len = *), intent(in) :: VAR_NAME
character (len = *), intent(in) :: DATA_unitS
real :: LATS(NLATS), LONS(NLONS)
integer :: dimids(2), start(2), cont(2)
character (LEN = *), parameter :: LAT_NAME = "lat"
character (LEN = *), parameter :: LON_NAME = "lon"
character (LEN = *), parameter :: unitS = "units"
character (LEN = *), parameter :: LAT_unitS = "degN"
character (LEN = *), parameter :: LON_unitS = "degE"
real, dimension(:,:) :: DATA_WROUT
call check(nf90_create(FILE_NAME, nf90_clobber, ncid))
call check(nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
call check(nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
call check(nf90_def_var(ncid, LAT_NAME, NF90_real, lat_dimid, lat_varid) )
call check(nf90_def_var(ncid, LON_NAME, NF90_real, lon_dimid, lon_varid) )
call check( nf90_put_att(ncid, lat_varid, unitS, LAT_unitS) )
call check( nf90_put_att(ncid, lon_varid, unitS, LON_unitS) )
dimids = (/ lon_dimid, lat_dimid /)
call check( nf90_def_var(ncid, VAR_NAME, NF90_real, dimids, data_varid) )
call check( nf90_put_att(ncid, data_varid, unitS, DATA_unitS) )
call check( nf90_enddef(ncid) )
call check( nf90_put_var(ncid, lat_varid, LATS) )
call check( nf90_put_var(ncid, lon_varid, LONS) )
start = (/ 1, 1 /)
cont = (/ NLONS, NLATS /)
call check(nf90_put_var(ncid,data_varid,DATA_WROUT,start,cont) )
call check( nf90_close(ncid) )
end subroutine MAKEDAT2D
!===============================================================================
! === NETCDF 3D-DATA read subroutine ===========================================
!===============================================================================
subroutine GETDAT3D(FILE_NAME,VAR_NAME,DATA_OUT,NLONS,NLATS,NRECS)
implicit none
integer :: lat, lon, rec, i, ncid
character (len = *), intent(in) :: FILE_NAME
character (len = *), intent(in) :: VAR_NAME
real :: lats(NLATS), lons(NLONS)
integer :: lon_varid, lat_varid, data_varid
integer, intent(in)	:: NRECS, NLATS, NLONS
integer			:: start(3), cont(3)
real, dimension(:,:,:), allocatable :: DATA_OUT
allocate(DATA_OUT(NLONS, NLATS, NRECS))
call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
call check( nf90_inq_varid(ncid, VAR_NAME, data_varid) )
cont = (/ NLONS, NLATS, NRECS /)
start = (/ 1, 1, 1 /)
call check( nf90_get_var(ncid, data_varid,DATA_OUT,start,cont))
call check( nf90_close(ncid) )
end subroutine GETDAT3D
!===============================================================================
! === NETCDF 2D-DATA read subroutine ===========================================
!===============================================================================
subroutine GETDAT2D(FILE_NAME,VAR_NAME,DATA_OUT,NLONS,NLATS)
implicit none
integer :: lat, lon, rec, i, ncid
character (len = *), intent(in) :: FILE_NAME
character (len = *), intent(in) :: VAR_NAME
integer, intent(in) :: NLATS, NLONS
integer :: lon_varid, lat_varid, data_varid
integer	:: start(2), cont(2)
real :: lats(NLATS), lons(NLONS)
real, dimension(:,:), allocatable :: DATA_OUT
allocate(DATA_OUT(NLONS, NLATS))
call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
call check( nf90_inq_varid(ncid, VAR_NAME, data_varid) )
cont = (/ NLONS, NLATS /)
start = (/ 1, 1 /)
call check( nf90_get_var(ncid, data_varid,DATA_OUT,start,cont))
call check( nf90_close(ncid) )
end subroutine GETDAT2D
!===============================================================================
! === NETCDF ERROR-CHECK subroutine ============================================
!===============================================================================
subroutine CHECK(status)
implicit none
integer, intent ( in) :: status
if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 1
endif
end subroutine CHECK
!===============================================================================
! === PROMPT subroutine ========================================================
!===============================================================================
subroutine WELLCOME
implicit none
!call system('clear')
write(*,*) '|*************************************************************|'
write(*,*) '|========   MODELO TERMODINAMICO DEL CLIMA GLOBAL   ==========|'
write(*,*) '|===================   CCA-UNAM 2015  ========================|'
write(*,FORMAT3)
end subroutine
end program mtc_simula
