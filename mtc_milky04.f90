PROGRAM mtc_milky04
USE netcdf
!include "/opt/local/lib/netcdf.inc"
IMPLICIT NONE
!===============================================================================
! GLOBAL THERMODYNAMIC CLIMATE MODEL ===========================================
!===============================================================================
! CCA-UNAM 2016                                                                =
! mtc_global.f90                                                               =
!                                                                              =
!                                                                              =
!                                                                              =
!===============================================================================
! FLAGS SWITCH:
! NAME     DESCRIPTION                                         BOOLEAN
!
! cldera   Cloud cover from ERA-Interim                        (1)
!          Cloud cover from NCEP                               (0)
! milan    Solar Radiation at TOA from Milankovitch            (1)
!          Solar Radiation at TOA from NCEP                    (0)
! g22p     binary field of land/sea ocean                      (1)
!          binary field of land/sea cont                       (0)
! lh       Latent Heat of Condensation from calorMON.f90       (1)
!          Latent Heat of Condensation from NCEP precipitation (0)
! xk0      Cloud Transmisivity from CALCULATION                (1)
!          Cloud Transmisivity from Adem 1962                  (0)
! a20      Atmospheric absorptivity from CALCULATION           (1)
!          Atmospheric absorptivity from Adem 1962             (0)
! a30      Cloud absorptivity from CALCULATION                 (1)
!          Cloud absorptivity from Adem 1962                   (0)
!===============================================================================
INTEGER :: cldera = 0 , milan = 1, g22p = 1
INTEGER :: lh = 1 , xk0 = 1 , a20 = 1 , a30 = 1
INTEGER :: i , j , m , n , ext
INTEGER :: icase, casen, casedn, reg
INTEGER :: currentdate, currentmes, currentyear
INTEGER :: fordate, formes, foryear
CHARACTER( LEN = 6 ) :: labeldate
!============================ Constants parameters =============================
REAL, PARAMETER :: cc      = 0.837e3           ! Spec. Heat continent
REAL, PARAMETER :: cs      = 4.1855e3          ! [J Kg^-1 K^-1] Spec. Heat water
REAL, PARAMETER :: cv      = 1.004e3           ! [J Kg^-1 K^-1]
REAL, PARAMETER :: dt      = 2.68e6            ! timestep
REAL, PARAMETER :: d2      = 2.01              ! cloud - G5 coefficient
!REAL, PARAMETER :: erre    = 4e-5
REAL, PARAMETER :: exentri = 0.06              ! Earth's orbital excentricity
REAL, PARAMETER :: g       = 9.81              ! [M S^-2] Grav. ac.
REAL, PARAMETER :: gamma   = 6.5e-3            ! [DegC M^-1] Lapse rate
REAL, PARAMETER :: h0      = 11e3              ! [M] Height at 266MB
REAL, PARAMETER :: h5      = 5.6e3             ! [M] Height at 500MB
REAL, PARAMETER :: hs0     = 60.               ! [M] Ocean mixed layer deep
REAL, PARAMETER :: hc      = 0.!2.5            ! [M] continent deep
REAL, PARAMETER :: kaust   = 2e6               ! [K S^-1] Austausch coeficient
REAL, PARAMETER :: k2      = 2.6838e-04        !
REAL, PARAMETER :: k3      = 2.6838e-04        ! 26.8380051 hemispheric
REAL, PARAMETER :: k4b     = 5.1844e-04        ! 4.05013077E-02 hemispheric
REAL, PARAMETER :: om      = 7.2921159e-5      ! Earth's angular velocity
REAL, PARAMETER :: p0      = 2.66e4            ! [Pa] Pressure at 266MB
REAL, PARAMETER :: p5      = 5.00e4
REAL, PARAMETER :: pi      = 3.141592653589793
REAL, PARAMETER :: r       = 6.37e6            ! [M] Earth's radius
REAL, PARAMETER :: rd      = 287.053           ! Gas constant for dry air.
REAL, PARAMETER :: rhos    = 1.e3              ! [Kg M^-3] Water density
REAL, PARAMETER :: rhoc    = 1.43e3            ! [Kg M^-3] Soil density
REAL, PARAMETER :: s0      = 1366              ! [W m-2] Solar constant
REAL, PARAMETER :: t0      = 216.5             ! [K] Temperature AT 266MB
REAL, PARAMETER :: tc0     = 261.              ! [K] Cloud Temperature
REAL, PARAMETER :: ts0     = 288.
!=========================== Calculated parameters =============================
REAL, PARAMETER :: a       = 5.*h0/(4.*t0+gamma*h0)
REAL, PARAMETER :: ds      = cs * rhos * hs0
REAL, PARAMETER :: dc      = cc * rhoc * hc
REAL, PARAMETER :: eps     = 23.5*pi/180.
REAL, PARAMETER :: grb     = g / ( rd * gamma )
REAL, PARAMETER :: ta0     = ts0
REAL, PARAMETER :: tm0     = t0 + gamma*h0/2.
REAL, PARAMETER :: f2      = cv*p0/g*((ta0/t0)**grb-1)      ! [7295820.50] Cv*a0
REAL, PARAMETER :: b2      = cv*t0/(gamma*(grb+1))*(1-(ta0/t0)**(grb+1))! FF8,F8
REAL, PARAMETER :: b3      = p0/t0 * (1-grb) * b2 - g / gamma * f2      ! F8P
REAL, PARAMETER :: b1      = b3 - p0/t0*b2                              ! F8PP
!=========================== Time & space parameters ===========================
INTEGER, PARAMETER       :: nlon    = 192                  ! longitude dim
INTEGER, PARAMETER       :: nlat    = 94                   ! latitude dim
INTEGER, PARAMETER       :: nsim    = 12                  ! forecast extension
INTEGER, PARAMETER       :: dateini = 199706              ! first forecast date
INTEGER, PARAMETER       :: dateend = 199706               ! last  forecast date
INTEGER, PARAMETER       :: yearini = int( dateini/100. )
INTEGER, PARAMETER       :: yearend = int( dateend/100. )
INTEGER, PARAMETER       :: mesini  = dateini - yearini*100
INTEGER, PARAMETER       :: mesend  = dateend - yearend*100
INTEGER, PARAMETER       :: ntime   = (yearend-yearini)*12 + (mesend-mesini) + 1
INTEGER, DIMENSION(12)   :: prev    = (/ 12,1,2,3,4,5,6,7,8,9,10,11 /)
INTEGER, DIMENSION(12)   :: season  = (/  4,4,1,1,1,2,2,2,3,3, 3, 4 /)
INTEGER, DIMENSION(nsim) :: mnth
!=========================== Declare data arrays ===============================
REAL :: t1, t2
!REAL :: f16, f22
REAL :: f30, f30p, f31, f32, f32p, f33, f34, f34p, f35, f36, f37
REAL :: ocean, cont
REAL :: f38, f41, f42, f43, f44, f46, f47, f48, f49
REAL :: f68, f68p, f69, f70, f70p
REAL :: f71, f71p, f72, f73, f75, f76, f77, f78
REAL :: f80, f81, f82, f83, f86, xlapn
REAL :: fcor, fg, fr, erre
REAL :: tt, tlower, tupper, plower, rholow, pr
REAL, DIMENSION(nlon)           :: lon                   ! longitude vector
REAL, DIMENSION(nlat)           :: lat, phi              ! latitude vector
REAL, DIMENSION(nlon*nlat)      :: landmask
REAL, DIMENSION(nlon*nlat)      :: qq0, xk, a1, a2, a3, sw, albs
REAL, DIMENSION(nlon*nlat)      :: esn, esdn, etn, etdn
REAL, DIMENSION(nlon*nlat)      :: cld, clddn, hum
REAL, DIMENSION(nlon*nlat)      :: g2n, g3n, g5n, g2dn, g3dn, g5dn
REAL, DIMENSION(nlon*nlat)      :: d7, d8
REAL, DIMENSION(nlon*nlat)      :: ap, bp, cp, dp, ep        ! vars clap's coef.
REAL, DIMENSION(nlon*nlat)      :: tmpn, tspn, tspdn     ! vars from input files
REAL, DIMENSION(nlon*nlat)      :: tspp, tmpp            ! vars primas
REAL, DIMENSION(nlon*nlat)      :: h250, h500, t500
REAL, DIMENSION(nlon*nlat)      :: tpn, tmn, tmdn, tsn, tsdn
REAL, DIMENSION(nlon*nlat)      :: van, vmn, umn!, umdn, vmdn, um, vm
REAL, DIMENSION(nlon*nlat)      :: f72p, f73p, f75p, f76p, f77p, f79
REAL, DIMENSION(nlon*nlat)      :: f89, f90
REAL, DIMENSION(nlon,nlat)      :: f84, f85, f87, f88, f87a, f88a, f89a, f90a
REAL, DIMENSION(nlon,nlat)      :: tpm, senp, cosp
REAL, DIMENSION(nlon,nlat)      :: pp, pn, pc, pcn, pcdn, tn, hup, h500a
REAL, DIMENSION(nlon,nlat)      :: f84n, f85n
REAL, DIMENSION(nlon,nlat)      :: dpdx, dpdy, dtdx, dtdy, dtndx, dtndy
REAL, DIMENSION(nlon,nlat)      :: jac1, jac2
REAL, DIMENSION(nlon,nlat)      :: tm1, g5a, bp1, cp1, dp1, ep1, g5dnm
REAL, DIMENSION(nlon,nlat,nsim) :: tpmn,  tpmdn!, tpsn, tpsdn
REAL, DIMENSION(nlon,nlat,nsim) :: aux1, aux2
!=============================== Begin Model ===================================
CALL system('clear')                        ! clean terminal window ------------
CALL wellcome                               ! prompt wellcome ------------------
CALL cpu_time( t1 )                         ! begin time count -----------------
CALL load_clima                             ! load climatology -----------------
CALL lw_func(t0,ta0,ts0,f30,f30p,f31,f32,f32p,f33,f34,f34p,f35,f36,f37)
CALL expParameters                          ! prompt experiment parameters -----

DO n = 1 , ntime
    WRITE(labeldate,"(I6)") currentdate
    reg = (yearini-1948)*12+mesini+n-2
    IF ( reg < 1 ) THEN
        WRITE(*,*) '** Error: Invalid Record Number. Check for forecast date **'
        STOP
    ENDIF
    CALL clim_data( reg , currentmes , prev(currentmes) )
    WRITE(*,"(30('*'))")
    WRITE(*,"(4X,A16,X,I3,5X,'*')")  'BEGIN EXPERIMENT',n
    WRITE(*,"(30('*'))")
    WRITE(*,"(A21,X,I6,X,'*')") '       current date =',currentdate
    WRITE(*,"(A21,X,I6,X,'*')") '      preview month =',prev(currentmes)
    WRITE(*,"(A21,X,I6,X,'*')") '             season =',season(currentmes)
    WRITE(*,"(A21,X,I6,X,'*')") 'rec (from Jan 1948) =',reg
    WRITE(*,"(30('*'))")
    !mnth(n)=currentdate
    ! set forecast PARAMETERs
    fordate=currentdate
    foryear=int(fordate/100)
    formes=currentmes
    ! Experiment loop ----------------------------------------------------------
    DO ext = 1 , nsim
        CALL clim_data( reg , formes , prev(formes) )
        WRITE(*,"(33('_'))")
        WRITE(*,"(A25,X,I3,5X,('|'))") '        forecast number =',ext
        WRITE(*,"(A25,X,I6,X,('|'))")  '      forecast for date =',fordate
        WRITE(*,"(A25,X,I2,5X,('|'))") 'Initial condition month =',prev(formes)
        WRITE(*,"(A25,X,I3,4X,('|'))") '                    rec =',reg
        WRITE(*,"(33('_'),('|'),17('_'))")
        mnth(ext)=fordate
        ! --- Computes climatic values of pressure -----------------------------
        m = 1
        DO j = 1 , nlat
            DO i = 1 , nlon
                tt      = t500(m) - gamma*(h0 - h500(m))
                pn(i,j) = 500*100*(tt/t500(m))**grb
                tn(i,j) = tpn(m) - t0
                m       = m + 1
            ENDDO
        ENDDO

        ! Case normal and departure to normal loop -----------------------------
        DO icase = 1 , 2
            casen  = 2 - icase              ! Normal case flag
            casedn = icase - 1              ! Departure to normal case flag
            IF ( ext == 1 ) THEN
                WRITE(*,"(A45,6X,('|'))") &
                                 "Calculating T'm & T's from initial conditions"
                DO m = 1 , nlon*nlat
                    tmpp(m) = tmpn(m) - tm0
                    tspp(m) = (tspn(m) + tspdn(m)*casedn) - ts0
                ENDDO
            ELSE
                WRITE(*,"(A38,13X,('|'))") &
                                        "Calculating T'm & T's from MTCG output"
                m=1
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        tmpp(m) = tpmn(i,j,ext-1) + tpmdn(i,j,ext-1)*casedn
                        tspp(m) = (tspn(m) + tspdn(m)*casedn) - ts0
                        m       = m + 1
                    ENDDO
                ENDDO
            ENDIF
            ! Model Coeficients computed below ---------------------------------
            DO m = 1 , nlon*nlat
                cont   = landmask(m)
                ocean  = 1.-landmask(m)
                a1(m)  = qq0(m)*( 1.-( 1.-xk(m) )*cld(m) )*( 1.-albs(m) )
                !f16    = -t0/p0/gamma/(1.-grb)*(1-(ta0/t0)**(grb+1))
                !f22    = -rd*h0/p0
                f38    = ocean*g3n(m) + cont*( g3n(m)-casedn*(1.-d7(m))*esn(m))
                f41    = ocean*k4b*van(m)*0.981*casedn
                f42    = - hum(m)/0.981*f41 +                                  &
                         casedn*cont*(1.-d7(m))*k3*van(m)*(1.-hum(m))
                f43    = cont*casedn *( 1 - d7(m) )
                f44    = g2n(m) + cont*(1.-g22p)*d8(m)
                f46    = van(m)*(ocean*k2+cont*k3*g22p)*casedn
                f47    = cont*(1.-g22p)*d7(m)
                f48    = f34 + f34p*cld(m) + a1(m)
                f49    = - qq0(m)*(1.-xk(m))*(1.-albs(m))
                f68    = bp(m)*casedn/1e3
                !f68    = bp(m)*casedn
                f68p   = 0.
                f71    = (ds*ocean+dc*cont)/dt - f36*(1.-f47-f43) + f41 + f46
                f72    = ( ( ds*ocean + dc*cont )*tspp(m)/dt +                 &
                         f48*( 1. - f47 - f43 ) - f44 - f38 ) / f71
                f73    = ( f35*( 1. - f47 - f43 ) + f46 - f42 ) / f71
                f75    = ( 1. - f47 - f43 ) * ( f49 + f34p ) / f71
                f76    = ( f41 + f46 ) / f71
                f77    = ( f42 - f46 ) / f71
                f71p   = 1. - casedn*f75*d2*f68p
                f72p(m)= f72/f71p
                f73p(m)= f73/f71p
                f75p(m)= f75/f71p
                f76p(m)= (f76- f75*casedn*d2*f68p)/f71p
                f77p(m)= f77/f71p
                f78    = f32 + cld(m)*f32p + f46 + f47*f36 + f68p +            &
                         d2*casedn*f68p*( f30p + a3(m)*sw(m) )
                f79(m) = 1. + d2*casedn*( (f49+f34p)*f47 + f30p + a3(m)*sw(m)+ &
                         f78*f75p(m) )
                f80    = f37 - f46 + f47*f35 + f78*f73p(m) + f79(m)*f68
                f81    = f30 + f44 + f47*f48 + f72p(m)*f78 +                   &
                         cld(m)*( f30p + a3(m)*sw(m) ) + a2(m)*sw(m) +         &
                         g5n(m) + ( f78*f76p(m) - f46 -f68p -                  &
                         d2*casedn*f68p*( f30p + a3(m)*sw(m) ) )*tsn(m)*casedn &
                         + (f78*f77p(m) + f46 -f79(m)*f68 )*tmn(m)*casedn
                f82    = f2  !- casedn*b2*f16*p0/f22/t0
                f83    = -f2 !+ casedn*b2*f16*p0/f22
                f86    = 1.  !- rd*f16*b2/cv/h0/f22
                f89(m) = f82/f83/dt -f86*f80/f83
                f90(m) = -tmpp(m)*f82/f83/dt - f86*f81/f83
            ENDDO
            IF ( icase==1 ) THEN
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        pp(i,j) = pn(i,j)
                    ENDDO
                ENDDO
            ELSE
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        pp(i,j) = pn(i,j) + pcdn(i,j)
                    ENDDO
                ENDDO
            ENDIF
            ! Convert some vectors into matrices--------------------------------
            m = 1
            DO j = 1 , nlat
                DO i = 1 , nlon
                    tm1(i,j) = tmn(m)
                    g5a(i,j) = g5n(m)
                    bp1(i,j) = bp(m)
                    cp1(i,j) = cp(m)
                    dp1(i,j) = dp(m)
                    m = m + 1
                ENDDO
            ENDDO
            ! Calculate space derivates of P and T with a centered finite ------
            ! differences ------------------------------------------------------
            CALL deriv_sph_ct(nlon,nlat,lon,lat,pp,dpdx,dpdy)
            CALL deriv_sph_ct(nlon,nlat,lon,lat,tn,dtdx,dtdy)
            CALL deriv_sph_ct(nlon,nlat,lon,lat,tm1,dtndx,dtndy)
            CALL jacobin(tm1,tn+t0,jac1)
            CALL jacobin(tn+t0,pn,jac2)
            ! ------------------------------------------------------------------
            m = 1
            DO j = 1 , nlat
                DO i = 1 , nlon
                    phi(j)    = -lat(j)*pi/180.
                    senp(i,j) = SIN(phi(j))
                    cosp(i,j) = COS(phi(j))
                    fcor      = 2*om*SIN(phi(j))
                    erre      = 5e-4*EXP(-ABS(phi(j)/(10.*pi/180.)))
                    fg        = fcor / (fcor**2+erre**2)
                    fr        = erre / (fcor**2+erre**2)
                    xlapn     = tm1(i+1,j) + tm1(i-1,j) + tm1(i,j+1) +         &
                                tm1(i,j-1) - 4.*tm1(i,j)
                    f69       = cp(m)*casedn
                    f70       = dp(m)*casedn
                    f70p      = 0.
                    f84(i,j)  = b2*( fg*dpdy(i,j) + fr*dpdx(i,j) ) / f83  +    &
                                b1*fg*dtdy(i,j)/f83 + b1*fr*dtdy(i,j)/f83 +    &
                                b3*fr*dtdy(i,j)/f83 + b3*fg*dtdy(i,j)/f83 +    &
                                b3*fg*dtndy(i,j)/f83 + b3*fr*dtndy(i,j)/f83 !+  &
                                !b1*fg*dtndy(i,j)/f83 + b1*fr*dtndy(i,j)/f83
                    IF ( icase==1 ) f84n(i,j) = f84(i,j)
                    f85(i,j)  = b2*( fg*dpdx(i,j) - fr*dpdy(i,j) ) / f83  +    &
                                b1*fg*dtdx(i,j)/f83 + b1*fr*dtdx(i,j)/f83 +    &
                                b3*fr*dtdx(i,j)/f83 + b3*fg*dtdx(i,j)/f83 +    &
                                b3*fg*dtndx(i,j)/f83 + b3*fr*dtndx(i,j)/f83 !+  &
                                !b1*fg*dtndx(i,j)/f83 + b1*fr*dtndx(i,j)/f83
                    IF ( icase==1 ) f85n(i,j) = f85(i,j)
                    f87(i,j)  = -f79(m)*f69/f83
                    f88(i,j)  = -f79(m)*f70/f83
                    f90(m)    = f90(m) + f86*f79(m)*( f69*dtndx(i,j) +         &
                                f70*dtndy(i,j) + f70p*xlapn ) / f83
                    f87a(i,j) = -umn(m) + ( (f84(i,j)-f84n(i,j)) + f87(i,j) )*casedn
                    f88a(i,j) = -vmn(m) - kaust*senp(i,j)/cosp(i,j)/r +        &
                                ( (f85(i,j)-f85n(i,j)) + f88(i,j) )*casedn
                    !f87a(i,j) = f84(i,j) + f87(i,j)
                    !f88a(i,j) = f85(i,j) + f88(i,j) - kaust*senp(i,j)/cosp(i,j)/r
                    f89a(i,j) = f89(m)
                    f90a(i,j) = f90(m) + 1./(cosp(i,j)*r**2)*( b3*jac1(i,j) + b2*jac2(i,j) ) /f83
                    m         = m + 1
                ENDDO
            ENDDO
            ! === CALL LAMBDA ==================================================
            CALL lambda( tpm , F90a , F89a , nlon , nlat )
            ! === CALL RELAX ===================================================
            CALL relax( icase, tpm, f90a, f89a, f87a, f88a, cosp, nlon, nlat )

            IF ( icase==1 ) THEN
                m = 1
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        hup(i,j) = h250(m)
                        m = m + 1
                    ENDDO
                ENDDO
            ELSE
                m = 1
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        hup(i,j) = h250(m) + a * ( tpm(i,j) - tmn(m) )
                        m = m + 1
                    ENDDO
                ENDDO
            ENDIF
            m = 1
            DO j = 1 , nlat
                DO i = 1 , nlon
                    tt         = tpm(i,j) + tm0 - gamma*h0/2
                    tupper     = tt - gamma * ( hup(i,j) - h0 )
                    tlower     = tt - gamma * ( h0- h5)
                    pr         = p0*(tt/tupper)**grb
                    plower     = pr*(tlower/tt)**grb
                    rholow     = plower/rd/tlower
                    h500a(i,j) = h5 + (plower -p5)/rholow/g
                    m = m + 1
                ENDDO
            ENDDO
            m = 1
            DO j = 1 , nlat
                DO i = 1 , nlon
                    tt        = tpm(i,j) + t0
                    tupper    = tt - gamma * ( hup(i,j) - h0 )
                    tlower    = tt + gamma * h0
                    pc(i,j)   = p0*(tt/tupper)**grb
                    m = m + 1
                ENDDO
            ENDDO
            IF ( icase == 1 ) THEN
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        pcn(i,j) = pc(i,j)
                        pcdn(i,j) = 0.
                    ENDDO
                ENDDO
            ELSE
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        pcdn(i,j) = pc(i,j) - pcn(i,j)
                    ENDDO
                ENDDO
            ENDIF
            DO j = 2 , nlat-1
                DO i = 2 , nlon-1
                    g5dnm(i,j) = g5a(i,j) + &
                                 casedn*( bp1(i,j)*( tpm(i,j) - tm1(i,j) ) +   &
                                 cp1(i,j)*( tpm(i-1,j) - tpm(i+1,j) -          &
                                 tm1(i-1,j) + tm1(i+1,j) ) +                   &
                                 dp1(i,j)*( tpm(i,j-1) - tpm(i,j+1) -          &
                                 tm1(i,j-1) + tm1(i,j+1) ) +                   &
                                 ep1(i,j)*(tpm(i-1,j) - tpm(i+1,j) -           &
                                 tpm(i,j-1) - tpm(i,j+1) - 4*tpm(i,j) -        &
                                 tm1(i-1,j) - tm1(i+1,j) - tm1(i,j-1) -        &
                                 tm1(i,j+1) + 4*tm1(i,j)) )
                ENDDO
            ENDDO
            m = 1
            DO j = 1 , nlat
                DO i = 1 , nlon
                    IF ( g5dnm(i,j) < 0. ) g5dnm = 0.
                    cont     = landmask(m)
                    ocean    = 1. - landmask(m)
                    tmdn(m)  = tpm(i,j)
                    tsdn(m)  = f72p(m) + f73p(m)*tmdn(m) +                     &
                               f76p(m)*tsn(m) + f77p(m)*tmn(m) +               &
                               f75p(m)*d2*casedn*( g5dnm(i,j) - g5n(m) )

                    g5dn(m)  = g5dnm(i,j) + casedn*ap(m)*( tsdn(m) - tsn(m) )
                    clddn(m) = cld(m) + casedn*d2*( g5dn(m) - g5n(m) )
                    a1(m)    = qq0(m)*( 1.-( 1.-xk(m) )*clddn(m) )*(1.-albs(m))
                    esdn(m)  = a1(m) + f34 + f34p*clddn(m) + f35*tmdn(m)+      &
                               f36*tsdn(m)
                    etdn(m)  = f30 + f30p*clddn(m) + f37*tmdn(m) +             &
                               (f32 + f32p*clddn(m))*tsdn(m) +                 &
                               sw(m)*(clddn(m)*a3(m)+a2(m))
                    g2dn(m)  = ( ocean + cont*g22p )*g2n(m) + casedn*van(m)*   &
                               ( ocean*k2 + cont*k3*g22p )*                    &
                               ( tsdn(m) - tsn(m) - tmdn(m) + tmn(m) )
                    g3dn(m)  = ocean*( g3n(m) + casedn*k4b*van(m)*(            &
                               0.981*(tsdn(m)-tsn(m))-hum(m)*(tmdn(m)-tmn(m))))&
                               + cont*( g3n(m) + casedn*( 1. - d7(m) ) *       &
                               ( esdn(m) - esn(m) +                            &
                               k3*van(m)*( 1.-hum(m) )*( tmdn(m)-tmn(m) )))
                    m        = m + 1
                ENDDO
            ENDDO
            IF ( icase == 1 ) THEN
                m = 1
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        tmn(m) = tmdn(m)
                        tsn(m) = tsdn(m)
                        etn(m) = etdn(m)
                        esn(m) = esdn(m)
                        g2n(m) = g2dn(m)
                        g3n(m) = g3dn(m)
                        g5n(n) = g5dn(m)
                        tpmn(i,j,ext) = tmn(m)
                        !tpsn(i,j,ext) = tsn(m)
                        !aux1(i,j,ext) = f84(i,j)   ! PRINT OUT NORMAL VARIABLE
                        m      = m + 1
                    ENDDO
                ENDDO
            ELSE
                m = 1
                DO j = 1 , nlat
                    DO i = 1 , nlon
                        tpmdn(i,j,ext) = tmdn(m) - tmn(m)
                        !tpsdn(i,j,ext) = tsdn(m) - tsn(m)
                        aux2(i,j,ext)  = (f84(i,j)-f84n(i,j))   ! PRINT OUT ANOMAL VARIABLE
                        m = m + 1
                    ENDDO
                ENDDO
            ENDIF
        ENDDO ! ******************************************* END OF icase DO LOOP
        ! Preparing date PARAMETERs for new simulation -------------------------
        reg=reg+1
        formes=formes+1
        IF ( mod(formes,13)==0 ) THEN
            foryear=foryear+1
            formes=1
        ENDIF
        fordate=foryear*100+formes
    ENDDO
    ! preparing date PARAMETERs for new initial conditions ---------------------
    currentmes=currentmes+1
    IF (mod(currentmes,13)==0) THEN
        currentyear=currentyear+1
        currentmes=1
    ENDIF
    currentdate=currentyear*100+currentmes
ENDDO! ****************************************************** END OF ext DO LOOP
WRITE(*,"(51('_'),'|')")
WRITE(*,"(A35,X,I3)") 'Writing NetCDF files for experiment',n-1
! writing Mid-troposphere temperature
WRITE(*,*) 'Writing Tm climatology [dec C]'
CALL makedat3d('output/TPMN'//labeldate//'.NC' , 'TPMN' ,  tpmn+tm0- 273.15 ,  &
               'Monthly Long-Term Mean of Mid-troposphere temperature' ,       &
               'degC' , 'Mid-troposphere temperature' , 'MTCG output' ,        &
               'Troposphere' , lon , lat , REAL(mnth) , nlon , nlat , nsim )
WRITE(*,*) 'Writing Tm anomaly [dec C]'
CALL makedat3d('output/TPMDN'//labeldate//'.NC' , 'TPMDN' ,  tpmdn ,           &
               'Monthly Anomaly of Mid-troposphere temperature' ,              &
               'degC' , 'Mid-troposphere temperature' , 'MTCG output' ,        &
               'Troposphere' , lon , lat , REAL(mnth) , nlon , nlat , nsim )
!! writing surface temperature
!CALL makedat3d('output/TPSN'//labeldate//'.NC','TPSN',tpsn+ts0-273.15,         &
!               'Monthly Long-Term Mean of Surface temperature' ,               &
!               'degC' , 'Surface temperature' , 'MTCG output' ,                &
!               'Surface' , lon , lat , REAL(mnth) , nlon , nlat , nsim )
!WRITE(*,*) 'Writing TS anomaly [dec C]'
!CALL makedat3d('output/TPSDN'//labeldate//'.NC','TPSDN', tpsdn ,               &
!               'Monthly Anomaly of Surface temperature' ,                      &
!               'degC' , 'Surface temperature' , 'MTCG output' ,                &
!               'Surface' , lon , lat , REAL(mnth) , nlon , nlat , nsim )
CALL makedat3d('output/VARN'//labeldate//'.NC','AUX', aux1 ,                   &
               'n/a' , 'n/a' , 'n/a' , 'n/a' , 'n/a',                          &
                lon , lat , REAL(mnth) , nlon , nlat , nsim )
CALL makedat3d('output/VARDN'//labeldate//'.NC','AUX', aux2 ,                  &
               'n/a' , 'n/a' , 'n/a' , 'n/a' , 'n/a',                          &
                lon , lat , REAL(mnth) , nlon , nlat , nsim )

CALL cpu_time( t2 )
WRITE(*,"(38('='))")
WRITE(*,"(A21,X,F7.3,X,A8)") 'End of experiments in',t2-t1,'seconds.'
WRITE(*,"(38('='))")
CONTAINS
!===============================================================================
! == DERIVATIVE SUBROUTINE =====================================================
! == Centred finite differences in spherical coordinates =======================
!===============================================================================
SUBROUTINE deriv_sph_ct(nlon,nlat,lon,lat,f,dfdx,dfdy)
REAL :: phi, dphi, dlam, northpole, southpole
INTEGER :: i , j, nlon, nlat
REAL, PARAMETER :: a  = 6.37e6
REAL, PARAMETER :: pi = 3.141592
REAL, DIMENSION(nlon) :: lon
REAL, DIMENSION(nlat) :: lat
REAL, DIMENSION(nlon,nlat) :: f, dfdx, dfdy
DO j = 1 , nlat
    phi = lat(j)*pi/180.
    DO i = 2 , nlon-1
        dlam = ( lon(i+1) - lon(i-1) ) * pi/180.
        dfdx(i,j) = 1/(a*cos(phi)) * ( f(i+1,j) - f(i-1,j)  ) / (2*dlam)
    ENDDO
    dfdx(nlon,j) = 1/(a*cos(phi)) * ( f(1,j) - f(nlon-1,j)  ) / (2*dlam)
    dfdx(1,j)    = 1/(a*cos(phi)) * ( f(2,j) - f(nlon,j)  ) / (2*dlam)
ENDDO
DO j = 2 , nlat-1
    dphi = ( lat(j+1) - lat(j-1) ) * pi/180.
    DO i = 1 , nlon
        dfdy(i,j) = 1/a * ( f(i,j+1) - f(i,j-1) ) / (dphi)
    ENDDO
ENDDO
northpole=SUM(dfdy(:,2))/nlon
southpole=SUM(dfdy(:,nlat-1))/nlon
DO i = 1 , nlon
    dfdy(i,1)    = northpole
    dfdy(i,nlat) = southpole
ENDDO
END SUBROUTINE deriv_sph_ct
!===============================================================================
! === LOAD_CLIMA SUBROUTINE ====================================================
!===============================================================================
SUBROUTINE load_clima
IMPLICIT NONE
CHARACTER (LEN = 14), PARAMETER :: rutaclim='input/monthly/'
WRITE(*,*) 'Loading climatology'
WRITE(*,*) 'Binary field for oceans and continents'
OPEN(UNIT=30,FILE=rutaclim//'GRID.DATA',STATUS='OLD')
OPEN(UNIT=31,FILE=rutaclim//'LAND.DATA',STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Clear Sky Downward Shortwave Radiation at surface'
OPEN(UNIT=33,FILE=rutaclim//'Qq0.DATA', STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Albedo Surface'
OPEN(UNIT=36,FILE=rutaclim//'ALBS.DATA',STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Sensible Heat Net Flux at surface'
OPEN(UNIT=37,FILE=rutaclim//'G2.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Latent Heat Net Flux at surface'
OPEN(UNIT=38,FILE=rutaclim//'G3.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) '-- Mid-troposphere temperature derivated from NCEP'
OPEN(UNIT=40,FILE=rutaclim//'TM.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Surface Temperature'
OPEN(UNIT=43,FILE=rutaclim//'TS.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Top of the Layer Temperature'
OPEN(UNIT=44,FILE=rutaclim//'TPN.DATA', STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Surface Temperature Anomaly (NCEP 1948-2015)'
OPEN(UNIT=45,FILE=rutaclim//'TSDN.DATA',STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Surface Wind Speed'
OPEN(UNIT=47,FILE=rutaclim//'VAN.DATA', STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Meridional tropospheric wind'
OPEN(UNIT=48,FILE=rutaclim//'VM.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Zonal tropospheric wind'
OPEN(UNIT=49,FILE=rutaclim//'UM.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'Relative Humidity'
OPEN(UNIT=59,FILE=rutaclim//'HUM.DATA', STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) 'D7 & D8 coefficients'
OPEN(UNIT=60,FILE=rutaclim//'D7.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(UNIT=61,FILE=rutaclim//'D8.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) "Clapp's coefficients"
OPEN(UNIT=62,FILE=rutaclim//'AP.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(UNIT=63,FILE=rutaclim//'BP.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(UNIT=64,FILE=rutaclim//'CP.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(UNIT=65,FILE=rutaclim//'DP.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(UNIT=66,FILE=rutaclim//'EP.DATA',  STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
WRITE(*,*) '--  temperature & height ERA (500mb)'
OPEN(unit=67,FILE=rutaclim//'TEMP500ERA.DATA',STATUS='UNKNOWN',            &
     ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(unit=68,FILE=rutaclim//'HGT500ERA.DATA',STATUS='UNKNOWN',             &
     ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
OPEN(unit=69,FILE=rutaclim//'H250.DATA',STATUS='UNKNOWN', ACCESS='DIRECT', &
     FORM='UNFORMATTED',RECL=nlon*nlat*8)
IF ( xk0 == 1 ) THEN
    WRITE(*,*) '-- XK Derivated from NCEP & ERA reanalisys'
    OPEN(UNIT=34,FILE=rutaclim//'XK.DATA',STATUS='UNKNOWN',     &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ELSE
    WRITE(*,*) '-- XK obteined by Adem'
ENDIF
IF ( cldera == 0 ) THEN
    OPEN(UNIT=32,FILE=rutaclim//'CLD.DATA',STATUS='UNKNOWN',    &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
    WRITE(*,*) '-- Total cloud cover from NCEP reanalysis'
ELSE
    OPEN(UNIT=32,FILE=rutaclim//'CLDERA.DATA',STATUS='UNKNOWN', &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
    WRITE(*,*) '-- Total cloud cover from ERA reanalysis'
ENDIF
IF ( milan == 0 ) THEN
    WRITE(*,*) "-- Insolation at TOA from NCEP reanalysis"
    OPEN(UNIT=35,FILE=rutaclim//'SW.DATA',STATUS='UNKNOWN',     &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ELSE
    WRITE(*,*) "-- Insolation at TOA calculated with Milankovitch's formula"
ENDIF
IF ( lh == 0 ) THEN
    WRITE(*,*) "-- Latent heat of condensation from NCEP precipitation"
    OPEN(UNIT=39,FILE=rutaclim//'G5.DATA',STATUS='UNKNOWN',     &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ELSE
    WRITE(*,*) "-- Latent heat of condensation from Kuo's scheme'"
    OPEN(UNIT=39,FILE=rutaclim//'QC.DATA',STATUS='UNKNOWN',     &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ENDIF
IF ( a20 == 1 ) THEN
    !WRITE(*,*) '-- A2 Derivated from NCEP & ERA reanalisys'
    WRITE(*,*) '-- A2 Derivated from Precipitable water'
    OPEN(UNIT=41,FILE=rutaclim//'A2.DATA',STATUS='UNKNOWN',     &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ELSE
    WRITE(*,*) '-- A2 obteined by Adem'
ENDIF
IF ( a30 == 1 ) THEN
    WRITE(*,*) '-- A3 Derivated from NCEP & ERA reanalisys'
    OPEN(UNIT=42,FILE=rutaclim//'BB3.DATA',STATUS='UNKNOWN',    &
         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=nlon*nlat*8)
ELSE
    WRITE(*,*) '-- A3 obteined by Adem'
ENDIF
READ(30,*) lon
READ(30,*) lat
WRITE(*,*) ' - - SUCCESS - - '
END SUBROUTINE load_clima
!===============================================================================
! === CLIM_DATA SUBROUTINE =====================================================
!===============================================================================
SUBROUTINE clim_data( n , month , monthi )
IMPLICIT NONE
INTEGER :: n , month , monthi
READ( 31 , REC = month ) landmask
READ( 32 , REC = month ) cld
READ( 33 , REC = month ) qq0
IF ( xk0 == 1 ) THEN
     READ( 34 , REC = month ) xk
ELSE
    DO m=1,nlon*nlat
        xk(m)=xkadem(int(m/nlon)+1.)
    ENDDO
ENDIF
IF ( a20 == 1 ) THEN
    READ( 41 , REC = month )  a2
ELSE
    DO m=1,nlon*nlat
        a2(m)=a2adem(int(m/nlon)+1.,season(currentmes))
    ENDDO
ENDIF
IF ( a30 == 1 ) THEN
     READ( 42 , REC = month )  a3
ELSE
    DO m=1,nlon*nlat
        a3(m)=a3adem(int(m/nlon)+1.,season(currentmes))
    ENDDO
ENDIF
IF ( milan == 1 ) THEN
    CALL inso_toa( month , s0 , eps, exentri,  sw )
ELSE
    READ( 35 , REC = month ) sw
ENDIF
READ( 36 , REC = month  )  albs
READ( 37 , REC = month  )  g2n
READ( 38 , REC = month  )  g3n
READ( 39 , REC = month  )  g5n
READ( 40 , REC = monthi )  tmpn    ! Tm normal
READ( 43 , REC = month  )  tspn    ! Ts normal prescirpcion oceano en mes actual
READ( 45 , REC = n      )  tspdn   ! Ts anomaly previa
READ( 44 , REC = month  )  tpn     ! T  normal
!READ( 46 , REC = MONTH )  TMDN    ! Tm anomaly previa
READ( 47 , REC = month  )  van
READ( 48 , REC = month  )  vmn
READ( 49 , REC = month  )  umn
READ( 67 , REC = month  )  t500
READ( 68 , REC = month  )  h500
READ( 69 , REC = month  )  h250
READ( 59 , REC = month  )  hum
READ( 60 , REC = month  )  d7
READ( 61 , REC = month  )  d8
READ( 62 , REC = month  )  ap
READ( 63 , REC = month  )  bp
READ( 64 , REC = month  )  cp
READ( 65 , REC = month  )  dp
READ( 66 , REC = month  )  ep
tmpn = tmpn + 273.16
tspn = tspn + 273.16
tpn  = tpn  + 273.16
cld  = cld  / 100.
END SUBROUTINE clim_data
!===============================================================================
! === LW_functions SUBROUTINE ==================================================
! ===== SUBROUTINE to display the coeficients from Long-wave radiation =========
! ===== balance in MTC. ( F30, F30p, F32, F34, F34p, F35, F36, f37 ).  =========
!===============================================================================
SUBROUTINE LW_func(t0,ta0,ts0,f30,f30p,f31,f32,f32p,f33,f34,f34p,f35,f36,f37)
IMPLICIT NONE
REAL, INTENT(IN)  :: t0 , ta0 , ts0
REAL, INTENT(OUT) :: f30, f30p, f31, f32, f32p, f33, f34, f34p, f35, f36, f37
REAL, PARAMETER   :: sig     = 5.6704e-8       ! [W M^-2 K^-4] Steffan-Boltzmann
REAL, PARAMETER   :: tap2tmp = 1. + a*gamma/2. ! [1.0599] T'A to T'M
REAL, PARAMETER   :: tp2tmp  = 1. - a*gamma/2. ! [0.94] T' to T'M
! ET
f30  =  f(t0) - sig*t0**4 + f(ta0) - sig*ta0**4 - f(ts0) + sig*ts0**4
f30p = -2.*f(tc0) + f(ts0)                                       ! [W M^-2]
f33  =  df(ta0)- 4.*sig*ta0**3                                   ! [W M^-2 K^-1]
f31  =  df(t0) - 4.*sig*t0**3 + f33                              ! [W M^-2 K^-1]
f32  =  4.*sig*ts0**3 - df(ts0)                                  ! [W M^-2 K^-1]
f32p =  df(ts0)                                                  ! [W M^-2 K^-1]
f37  =  f33*tap2tmp + (f31-f33)*tp2tmp         ! f33 +f31        ! [W M^-2 K^-1]
! ES
f34  =  sig*ta0**4 - f(ta0)-sig*ts0**4                           ! [W M^-2]
f34p =  f(tc0)                                                   ! [W M^-2]
f35  = -f33 * tap2tmp                                            ! [W M^-2 K^-1]
f36  = -4.*sig*ts0**3                                            ! [W M^-2 K^-1]
END SUBROUTINE LW_func
!===============================================================================
! === LW ABSORPTION FUNCTIONS ==================================================
!===============================================================================
REAL FUNCTION f(t)
IMPLICIT NONE
REAL, INTENT(IN) :: t
REAL :: abs1=0.0, abs2=0.067, abs3=0.802, abs4=0.931, abs5=0.310
f  = (1-abs1)*( ftl(t,12e-6)-ftl(t,08e-6))  + &
     (1-abs2)*( ftl(t,13e-6)-ftl(t,12e-6) ) + & ! CO2
     (1-abs3)*( ftl(t,14e-6)-ftl(t,13e-6) ) + & ! CO2
     (1-abs4)*( ftl(t,17e-6)-ftl(t,16e-6) ) + &
     (1-abs5)*( ftl(t,18e-6)-ftl(t,17e-6) )
END FUNCTION f
REAL FUNCTION df(t)
IMPLICIT NONE
REAL, INTENT(IN) :: t
REAL :: abs1=0.0, abs2=0.067, abs3=0.802, abs4=0.931, abs5=0.310
df = (1-abs1)*(dftl(t,12e-6)-dftl(t,08e-6)) + &
     (1-abs2)*(dftl(t,13e-6)-dftl(t,12e-6)) + &
     (1-abs3)*(dftl(t,14e-6)-dftl(t,13e-6)) + &
     (1-abs4)*(dftl(t,17e-6)-dftl(t,16e-6)) + &
     (1-abs5)*(dftl(t,18e-6)-dftl(t,17e-6))
END FUNCTION df
REAL FUNCTION ftl(t,lb)
IMPLICIT NONE
REAL, INTENT(IN) :: t, lb
REAL, PARAMETER  :: c       = 299792458.           ! [M S^-1] Light speed
REAL, PARAMETER  :: hpl     = 6.62606896e-34       ! [J S] Planck constant
REAL, PARAMETER  :: kb      = 1.3806503e-23        ! [J K^-1] Boltzmann constant
REAL, PARAMETER  :: c1      = 2*pi*hpl*c**2        ! [3.742E-16 W M^2]
REAL, PARAMETER  :: c2      = (hpl*c/kb)
ftl  = c1 * exp(-c2/(lb*t)) * ( t/(c2*lb**3) + 3.*(t/(c2*lb))**2 + &
       6.*(t**3/(lb*c2**3) + 6.*(t/c2)**4))
END FUNCTION ftl
REAL FUNCTION dftl(t,lb)
IMPLICIT NONE
REAL, INTENT(IN) :: t, lb
REAL, PARAMETER  :: c       = 299792458.           ! [M S^-1] Light speed
REAL, PARAMETER  :: hpl     = 6.62606896e-34       ! [J S] Planck constant
REAL, PARAMETER  :: kb      = 1.3806503e-23        ! [J K^-1] Boltzmann constant
REAL, PARAMETER  :: c1      = 2*pi*hpl*c**2        ! [3.742E-16 W M^2]
REAL, PARAMETER  :: c2      = (hpl*c/kb)
dftl = c1 * exp(-c2/(lb*t)) * ( 1./(t*lb**4) + 4./(c2*lb**3) + &
       12.*t/(lb*c2)**2 + 24.*(t**2)/(lb*c2**3) + (t**3/c2**4))
END FUNCTION dftl
!===============================================================================
! == SHORTWAVE RADIATION functions =============================================
!===============================================================================
REAL FUNCTION xkadem(latitude)
REAL, INTENT(IN) :: latitude
REAL, DIMENSION(20) :: ksa=(/ 0.35 , 0.35 , 0.34 , 0.34 , 0.33 , 0.33 , 0.32 , &
                              0.32 , 0.32 , 0.33 , 0.34 , 0.36 , 0.38 , 0.40 , &
                              0.45 , 0.50 , 0.55 , 0.55 , 0.55 , 0.55 /)
xkadem = ksa(int(abs(latitude)/5)+1)
END FUNCTION xkadem
REAL FUNCTION a2adem(latitude,season)
INTEGER :: season
REAL, INTENT(IN) :: latitude
REAL, DIMENSION(20) :: A22
IF ( season == 1 ) THEN
    a22=(/ 0.134 , 0.134 , 0.123 , 0.123 , 0.122 , 0.122 , 0.105 , &
           0.105 , 0.100 , 0.100 , 0.095 , 0.095 , 0.108 , 0.108 , &
           0.143 , 0.143 , 0.143 , 0.143 , 0.143 , 0.143 /)
ENDIF
IF ( SEASON == 2) THEN
    a22=(/ 0.143 , 0.143 , 0.145 , 0.145 , 0.145 , 0.145 , 0.142 , &
           0.142 , 0.135 , 0.135 , 0.132 , 0.132 , 0.147 , 0.147 , &
           0.169 , 0.169 , 0.164 , 0.164 , 0.164 , 0.164 /)
ENDIF
IF ( SEASON == 3) THEN
    a22=(/ 0.142 , 0.142 , 0.141 , 0.141 , 0.140 , 0.140 , 0.123 , &
           0.123 , 0.116 , 0.116 , 0.110 , 0.110 , 0.120 , 0.120 , &
           0.149 , 0.149 , 0. , 0. , 0. , 0. /)
ENDIF
IF ( SEASON == 4) THEN
    a22=(/ 0.137 , 0.137 , 0.134 , 0.134 , 0.126 , 0.126 , 0.115 , &
           0.115 , 0.104 , 0.104 , 0.1 , 0.1 , 0.114 , 0.114 ,     &
           0.167 , 0.167 , 0. , 0. , 0. , 0. /)
ENDIF
a2adem = a22(int(abs(latitude)/5)+1)
END FUNCTION a2adem
REAL FUNCTION A3ADEM(latitude,season)
INTEGER :: season
REAL, INTENT(IN) :: latitude
REAL, DIMENSION(20) :: a33
IF (season==1) THEN
    a33=(/ 0.034 , 0.034 , 0.030 , 0.030 , 0.030 , 0.030 , 0.032 , &
           0.032 , 0.033 , 0.033 , 0.035 , 0.035 , 0.042 , 0.042 , &
           0.038 , 0.038 , 0.033 , 0.033 , 0.033 , 0.033 /)
ENDIF

IF (season==2) THEN
    a33=(/ 0.034 , 0.034 , 0.038 , 0.038 , 0.036 , 0.036 , 0.033 , &
           0.033 , 0.033 , 0.033 , 0.035 , 0.035 , 0.037 , 0.037 , &
           0.033 , 0.033 , 0.034 , 0.034 , 0.034 , 0.034 /)
ENDIF

IF (season==3) THEN
    a33=(/ 0.031 , 0.031 , 0.032 , 0.032 , 0.034 , 0.034 , 0.035 , &
           0.035 , 0.042 , 0.042 , 0.037 , 0.037 , 0.050 , 0.050 , &
           0.030 , 0.030 , 0. , 0. , 0. , 0. /)
ENDIF

IF (season==4) THEN
    a33=(/ 0.033 , 0.033 , 0.032 , 0.032 , 0.039 , 0.039 , 0.036 , &
           0.036 , 0.044 , 0.049 , 0.049 , 0.039 , 0.039 , 0. ,    &
           0. , 0. , 0. , 0. , 0. , 0. /)
ENDIF
a3adem = a33(int(abs(latitude)/5)+1)
END FUNCTION a3adem
!===============================================================================
FUNCTION accday(meses,year)
INTEGER, DIMENSION(12) :: acc
INTEGER :: year, meses, mes, i, accday
acc=(/31,59,90,120,151,181,212,243,273,304,334,365/)
IF (mod(year,4)==0) THEN
    DO i = 2 , 12
        acc(i)=acc(i)+1
    ENDDO
ENDIF
mes=mod(meses,12)
IF (mes==0) mes=12
accday=acc(mes)
RETURN
END FUNCTION accday
!===============================================================================
SUBROUTINE jacobin(da,db,jacob)
REAL, DIMENSION(nlon,nlat) :: da, dadx, dady, db, dbdx, dbdy, jacob
CALL deriv_sph_ct(nlon,nlat,lon,lat,da,dadx,dady)
CALL deriv_sph_ct(nlon,nlat,lon,lat,db,dbdx,dbdy)
DO j=1,nlat
     DO i=1,nlon
         jacob(i,j) = dadx(i,j)*dbdy(i,j) - dady(i,j)*dbdx(i,j)
     ENDDO
ENDDO
END SUBROUTINE jacobin
!===============================================================================
! === MILANK subroutine ========================================================
! ===== Soubrutine to compute the latitudinal solar irradiance at the top of ===
! ===== the atmosphere (TOA) using milankovitch formula on a given day. ========
! ==============================================================================
!    IT = s0 (a0/r)^2 [ w1 sin(phi) sin(delta) + sin(w1) cos(phi) cos(delta) ]
!===============================================================================
SUBROUTINE milank( lats, nlat, s0, epsi, exen, ls , it)
IMPLICIT NONE
INTEGER :: j
REAL :: a0r, delta, w1, cosz, yr
REAL, PARAMETER  :: lsp = 90.*pi/180.
REAL, INTENT(IN) :: s0, epsi, exen
INTEGER, INTENT(IN) :: ls
INTEGER, INTENT(IN) :: nlat
REAL, DIMENSION(nlat) :: phi
REAL, DIMENSION(nlat), INTENT(IN)  :: lats
REAL, DIMENSION(nlat), INTENT(OUT) :: it
yr    = (ls*360./365.)*pi/180.
a0r   = ( 1. + exen*COS(yr-lsp) ) / ( 1. - exen**2 )
delta = ASIN( SIN(epsi)*SIN(yr-lsp) )
DO j = 1 , nlat
    phi(j) = lats(j)*pi/180.
    w1    = ACOS( - TAN(phi(j)) * TAN(delta) )
    cosz  = ( w1*SIN(phi(j))*SIN(delta) + SIN(w1)*COS(phi(j))*COS(delta) )
    it(j) = ( s0/pi ) * cosz*(a0r)**2
    IF ( ABS(phi(j)) .GT. pi/2.-ABS(delta) ) THEN
        it(j) = 0.02
        IF ( phi(j)*delta .GT. 0. ) THEN
            it(j) = s0*(a0r)**2*SIN(phi(j))*SIN(delta)
        ENDIF
    ENDIF
ENDDO
END SUBROUTINE milank
!===============================================================================
! === INSO_TOA subroutine ======================================================
! ===== Soubrutine to compute the monthly mean of solar irradiance at the top ==
! ===== of the atmosphere (TOA) using milankovitch formula on a given month. ===
!===============================================================================
SUBROUTINE inso_toa( month , s0 , eps, ex,  solar )
IMPLICIT NONE
REAL :: s0, eps, ex
INTEGER :: month , m ,  nini , nfin, nday, year, ls
REAL, DIMENSION(:,:), ALLOCATABLE  :: insol
REAL, DIMENSION(nlon,nlat) :: varsum
REAL, DIMENSION(nlon*nlat),INTENT(OUT) :: solar
IF( month == 1 ) THEN
    nini = 1
ELSE
    nini = accday(month-1,year)+1
ENDIF
nfin = accday(month,year)
nday = nfin - nini
ALLOCATE(insol(nlat,nday))
varsum = 0.
m = 1
DO ls = nini , nfin
    CALL milank( lat, nlat, s0, eps, ex, ls, insol(:,ls))
    DO j = 1 , nlat
        DO i = 1 , nlon
            varsum(i,j) = insol(j,ls) + varsum(i,j)
        ENDDO
    ENDDO
    m = m + 1
ENDDO
!write(*,*) nday
DO j = 1 , nlat
    DO i = 1 , nlon
        varsum(i,j) = varsum(i,j) / nday
    ENDDO
ENDDO
!CALL POLE(NLON,NLAT,varsum)
m = 1
DO j = 1 , nlat
    DO i = 1 , nlon
        solar(m) = varsum(i,j)
        m = m + 1
    ENDDO
ENDDO
DEALLOCATE(insol)
RETURN
END SUBROUTINE INSO_TOA
!===============================================================================
! === LAMBDA SUBROUTINE ========================================================
!===============================================================================
SUBROUTINE  lambda  ( tpm , f90 , f89, plon, plat )
IMPLICIT NONE
INTEGER,INTENT(IN) :: plon, plat
REAL, DIMENSION(plon,plat),INTENT(IN)  :: f90, f89
REAL, DIMENSION(plon,plat),INTENT(OUT) :: tpm
DO j = 1 , plat
    DO i = 1 , plon
        tpm(i,j) = - f90(i,j)/f89(i,j)
    ENDDO
ENDDO
CALL pole(plon,plat,tpm)
END SUBROUTINE lambda
!===============================================================================
! === POLE SUBROUTINE ==========================================================
! ===== Subroutine to average data in poles from the nearest data. =============
!===============================================================================
SUBROUTINE pole( xlon , xlat, var )
IMPLICIT NONE
INTEGER :: xlon , xlat
REAL :: tsum1 , tsum2
REAL, DIMENSION(xlon,xlat) :: var
tsum1 = 0.
tsum2 = 0.
DO i = 1 , xlon
    tsum1 = tsum1 + var(i,2)
    tsum2 = tsum2 + var(i,xlat-1)
ENDDO
DO i = 1 , xlon
    var(i,1) = tsum1/xlon
    var(i,xlat)= tsum2/xlon
ENDDO
END SUBROUTINE pole
!===============================================================================
! === RELAX SUBROUTINE =========================================================
!===============================================================================
SUBROUTINE relax( icase, tpm, f90, f89, f84, f85, cosf, plon, plat)
IMPLICIT NONE
REAL, DIMENSION(plon,plat) :: tpm, f90, f89, f84, f85, cosf
REAL :: ml , mf , dml , dmf , tpmn , tpmnb , tpmnr , tpmnl , tpmnt
REAL :: e1, e2, e3, e4, e5, e6, res, erro, ar
REAL :: s, w, deltaf, deltal, t1, t2
INTEGER :: plon, plat, icase, rcycle, rlimit
deltaf  = pi/plat
deltal  = 2.*pi/plon
ar      = 0.5
erro    = 5.0e-3
rlimit  = 9999
CALL cpu_time ( t1 )
DO rcycle = 1 , rlimit
    w = 0.
    DO j = 2 , plat-1
        DO i = 1 , plon
            IF ( i == 1 ) THEN
                tpmnl    =  tpm(plon,j)
            ELSE
                tpmnl    =  tpm(i-1,j)
            ENDIF
            IF( i == plon ) THEN
                tpmnr    =  tpm(1,j)
            ELSE
                tpmnr    =  tpm(i+1,j)
            ENDIF
            tpmnt    =  tpm(i,j+1)
            tpmnb    =  tpm(i,j-1)
            tpmn     =  tpm(i,j)
            ml       =  1./deltal/cosf(i,j)
            mf       =  1./deltaf
            dml      =  r/ml/2.
            dmf      =  r/mf/2.
            e1       =  kaust + f84(i,j)*dml
            e2       =  (mf/ml)**2*(kaust + f85(i,j)*dmf)
            e3       = -e1 + 2.*kaust
            e4       = -e2 + 2.*(mf/ml)**2*kaust
            e5       =  2.*kaust*(1. + (mf/ml)**2) - 4.*f89(i,j)*dml**2
            e6       =  4.*dml**2*f90(i,j)
            res      =  e1*tpmnr+e2*tpmnt+e3*tpmnl+e4*tpmnb-e5*tpmn+e6
            s        =  res*ar/e5
            IF ( abs( s ) > erro )  w = 1.
            tpm(i,j) =  tpm(i,j) + s
        ENDDO
    ENDDO
    CALL pole( plon , plat, tpm )
    IF ( w == 0. )  EXIT
    IF ( rcycle == rlimit ) WRITE(*,*) 'Did not converge'
ENDDO
CALL cpu_time ( t2 )
WRITE(*,"(I4,X,A18,X,F5.3,X,A16,X,I1,3X,('|'))") rcycle , &
                                        "Cycles to relax in", ( T2 - T1 ) ,    &
                                        "seconds for case" , icase
END SUBROUTINE relax
!===============================================================================
! == NETCDF 3D-DATA write SUBROUTINE ===========================================
!===============================================================================
SUBROUTINE makedat3d( file_name , var_name , data_wrout , ln_data , &
                      data_units , vardesc , dataset , lvldesc , lons , lats , &
                      recs , nlons , nlats , nrecs )
IMPLICIT NONE
INTEGER :: ncid
INTEGER :: lon_varid, lat_varid, rec_varid, data_varid
INTEGER :: lon_dimid, lat_dimid, rec_dimid
INTEGER, INTENT(IN) :: nrecs, nlats, nlons
CHARACTER (LEN = *), INTENT(IN) :: file_name
CHARACTER (LEN = *), INTENT(IN) :: var_name
CHARACTER (LEN = *), INTENT(IN) :: data_units
CHARACTER (LEN = *), INTENT(IN) :: ln_data
CHARACTER (LEN = *), INTENT(IN) :: vardesc
CHARACTER (LEN = *), INTENT(IN) :: dataset
CHARACTER (LEN = *), INTENT(IN) :: lvldesc
REAL    :: lats(nlats), lons(nlons), recs(nrecs)
INTEGER :: dimids(3)
INTEGER :: start(3), cont(3)
CHARACTER (LEN = *), PARAMETER :: lat_name  = "lat"
CHARACTER (LEN = *), PARAMETER :: lon_name  = "lon"
CHARACTER (LEN = *), PARAMETER :: rec_name  = "time"
CHARACTER (LEN = *), PARAMETER :: lat_units = "degree_N"
CHARACTER (LEN = *), PARAMETER :: lon_units = "degree_E"
CHARACTER (LEN = *), PARAMETER :: rec_units = "YYYYMM"
REAL, DIMENSION(:,:,:) :: data_wrout
CALL check(nf90_create(file_name, nf90_clobber, ncid) )
CALL check(nf90_def_dim(ncid, lat_name, nlatS, lat_dimid) )
CALL check(nf90_def_dim(ncid, lon_name, nlons, lon_dimid) )
CALL check(nf90_def_dim(ncid, rec_name, nrecs, rec_dimid) )
dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
CALL check(nf90_def_var(ncid, lat_name, nf90_real, lat_dimid, lat_varid) )
CALL check(nf90_def_var(ncid, lon_name, nf90_real, lon_dimid, lon_varid) )
CALL check(nf90_def_var(ncid, rec_name, nf90_real, rec_dimid, rec_varid) )
CALL check(nf90_def_var(ncid, var_name, nf90_real, dimids, data_varid) )
CALL check( nf90_put_att(ncid, lat_varid, 'units', lat_units) )
CALL check( nf90_put_att(ncid, lon_varid, 'units', lon_units) )
CALL check( nf90_put_att(ncid, rec_varid, 'units', rec_units) )
CALL check( nf90_put_att(ncid, data_varid, 'long_name', ln_data) )
CALL check( nf90_put_att(ncid, data_varid, 'units', data_units) )
CALL check( nf90_put_att(ncid, data_varid, 'add_offset', 0.0 ) )
CALL check( nf90_put_att(ncid, data_varid, 'scale_factor', 1.0 ) )
CALL check( nf90_put_att(ncid, data_varid, 'precision', 2 ) )
CALL check( nf90_put_att(ncid, data_varid, 'least_significant_digit', 1 ) )
CALL check( nf90_put_att(ncid, data_varid, 'var_desc', vardesc ) )
CALL check( nf90_put_att(ncid, data_varid, 'dataset', dataset ) )
CALL check( nf90_put_att(ncid, data_varid, 'level_desc', lvldesc ) )
CALL check( nf90_enddef(ncid) )
CALL check( nf90_put_var(ncid, lat_varid, lats) )
CALL check( nf90_put_var(ncid, lon_varid, lons) )
CALL check( nf90_put_var(ncid, rec_varid, recs) )
cont = (/ nlons, nlats, nrecs /)
start = (/ 1, 1, 1 /)
CALL check(nf90_put_var(ncid,data_varid,data_wrout,start,cont) )
CALL check( nf90_close(ncid) )
END SUBROUTINE makedat3d
!===============================================================================
! === NETCDF ERROR-CHECK SUBROUTINE ============================================
!===============================================================================
SUBROUTINE check(status)
IMPLICIT NONE
INTEGER, INTENT (IN) :: status
IF (status /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(status))
      STOP 1
ENDIF
END SUBROUTINE check
!===============================================================================
! === PROMPT subroutine ========================================================
!===============================================================================
SUBROUTINE wellcome
IMPLICIT NONE
WRITE(*,*) '|*************************************************************|'
WRITE(*,*) '|====   GLOBAL THERMODYNAMIC CLIMATE MODEL   =================|'
WRITE(*,*) '|===================   CCA-UNAM 2015  ========================|'
WRITE(*,'(X,"|",61X,"|")')
END SUBROUTINE wellcome
!===============================================================================
! === EXPERIMENT subroutine ====================================================
!===============================================================================
SUBROUTINE expParameters
IMPLICIT NONE
WRITE(*,*)
WRITE(*,"(A24,X,I6)") 'Initial date experiment:', dateini
WRITE(*,"(A24,X,I6)") '  Final date experiment:', dateend
WRITE(*,"(A24,X,I6)") '     Forecast extension:', nsim
WRITE(*,"(A24,X,I6)") '        No. experiments:', ntime
WRITE(*,*)
currentdate = yearini*100 + mesini
currentmes  = mesini
currentyear = yearini
END SUBROUTINE expParameters
END PROGRAM mtc_milky04
