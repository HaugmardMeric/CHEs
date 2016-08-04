! novembre 2013
! gère les problèmes liés au temps :
! calendrier Julien, base 60 ...
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE time

    use modparam
    use typetemps, only : date_min, date_sec, date_secPS

    implicit none

    private

    public  :: tempszero
    public  :: basetime
    public  :: difftime
    public  :: JDATE,GDATE


    ! le temps est défini ici, soit en :
    ! - annee, mois, jour, heure, minutes et secondes
    ! - annee, mois, jour, heure, minutes, secondes pour ondes P, secondes pour ondes S
    ! les ondes S sont donc TOJOURS associées à des ondes P

    ! -------                                                  --------    .
    interface basetime
      module procedure basetimeP, basetimePS   ! basetime pour differents types
    end interface basetime
    ! -------                                                  --------    .
    interface difftime
      module procedure difftimeP, difftimePS  ! difftime pour differents types
    end interface difftime
    ! -------                                                  --------    .


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine JDATE (JD,YEAR,MONTH,DAY)
    ! -------                                                  --------    .mh
    ! COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
    ! Reference: Fliegel, H. F. & van Flandern, T. C. 1968, Communications of the ACM, 11, 657.
    ! source : http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    ! for years AD 1801–2099 :
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (out) :: JD
    integer(KIND=wi), intent (in) :: YEAR,MONTH,DAY
    integer(KIND=wi) :: I,J,K
    ! -----------------------------------------------------------------    .
    I= YEAR
    J= MONTH
    K= DAY
    JD= K-32075+1461*(I+4800+(J-14)/12)/4+ &
        367*(J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4
    ! -----------------------------------------------------------------    .
  end subroutine

! ---------------------------------------------------------------------    .

  subroutine GDATE (JD,YEAR,MONTH,DAY)
    ! -------                                                  --------    .mh
    !COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY) GIVEN THE JULIAN DATE (JD)
    ! Reference: Fliegel, H. F. & van Flandern, T. C. 1968, Communications of the ACM, 11, 657.
    ! source : http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    ! for years AD 1801–2099 :
    ! -------                                                  --------    .
    implicit none
    integer(KIND=wi), intent (in) :: JD
    integer(KIND=wi), intent (out) :: YEAR,MONTH,DAY
    integer(KIND=wi) :: I,J,K,L,N
    ! -----------------------------------------------------------------    .
    L= JD+68569
    N= 4*L/146097
    L= L-(146097*N+3)/4
    I= 4000*(L+1)/1461001
    L= L-1461*I/4+31
    J= 80*L/2447
    K= L-2447*J/80
    L= J/11
    J= J+2-12*L
    I= 100*(N-49)+I+L
    YEAR= I
    MONTH= J
    DAY= K
    ! -----------------------------------------------------------------    .
  end subroutine

! ---------------------------------------------------------------------    .

  subroutine basetimePS(thedate)
    ! -------                                                  --------    .mh
    ! respect des bases 60 pour les sec et min, et 24 pour les heures
    ! ansi que du découpage des années en mois et jours avec prise en compte des années bisextiles
    ! les ondes P et S : type(date) identique -> secS peut-être >60 !!!
    ! -------                                                  --------    .
    implicit none
    type(date_secPS), intent (inout) :: thedate
    ! -----------------------------------------------------------------    .
    do while(thedate%secP.ge.60.0_wr)                                      ! pas plus de 60 sec dans une minute
      thedate%secP=thedate%secP-60.0_wr
      thedate%secS=thedate%secS-60.0_wr
      thedate%date%min=thedate%date%min+1
    enddo
    ! -------                                                  --------    .
    do while(thedate%secP.lt.0.0_wr)                                       ! pas moins de 0 sec dans une minute
      thedate%secP=thedate%secP+60.0_wr
      thedate%secS=thedate%secS+60.0_wr
      thedate%date%min=thedate%date%min-1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%min.ge.60)                                       ! pas plus de 60 min dans une heure
      thedate%date%min=thedate%date%min-60
      thedate%date%hour=thedate%date%hour+1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%min.lt.0)                                        ! pas moins de 0 min dans une heure
      thedate%date%min=thedate%date%min+60
    thedate%date%hour=thedate%date%hour-1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%hour.ge.24)                                      ! pas plus de 24 heures dans une journée
      thedate%date%hour=thedate%date%hour-24
      thedate%date%day=thedate%date%day+1
      call JDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
      call GDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%hour.lt.0)                                       ! pas moins de 0 heure dans une journée
      thedate%date%hour=thedate%date%hour+24
      thedate%date%day=thedate%date%day-1
      call JDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
      call GDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine basetimePS

! ---------------------------------------------------------------------    .

  subroutine basetimeP(thedate)
    ! -------                                                  --------    .mh
    ! respect des bases 60 pour les sec et min, et 24 pour les heures
    ! ansi que du decoupage des années en mois et jours avec prise en compte des années bisextiles
    ! -------                                                  --------    .
    implicit none
    type(date_sec), intent (inout) :: thedate
    ! -----------------------------------------------------------------    .
    do while(thedate%sec.ge.60.0_wr)                                       ! pas plus de 60 sec dans une minute
      thedate%sec=thedate%sec-60.0_wr
      thedate%date%min=thedate%date%min+1
    enddo
    ! -------                                                  --------    .
    do while(thedate%sec.lt.0.0_wr)                                        ! pas moins de 0 sec dans une minute
      thedate%sec=thedate%sec+60.0_wr
      thedate%date%min=thedate%date%min-1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%min.ge.60)                                       ! pas plus de 60 min dans une heure
      thedate%date%min=thedate%date%min-60
      thedate%date%hour=thedate%date%hour+1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%min.lt.0)                                        ! pas moins de 0 min dans une heure
      thedate%date%min=thedate%date%min+60
      thedate%date%hour=thedate%date%hour-1
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%hour.ge.24)                                      ! pas plus de 24 heure dans une journée
      thedate%date%hour=thedate%date%hour-24
      thedate%date%day=thedate%date%day+1
      call JDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
      call GDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
    enddo
    ! -------                                                  --------    .
    do while(thedate%date%hour.lt.0)                                       ! pas moins de 0 heure dans une journée
      thedate%date%hour=thedate%date%hour+24
      thedate%date%day=thedate%date%day-1
      call JDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
      call GDATE (thedate%date%Jday,thedate%date%year,thedate%date%month,thedate%date%day)
    enddo
    ! -----------------------------------------------------------------    .
  end subroutine basetimeP

! ---------------------------------------------------------------------    .

  subroutine difftimePS(deltaP,deltaS,thedate1,thedate2)
    ! -------                                                  --------    .mh
    ! difference relative entre deux temps absolus, résultat en sec
    ! ondes P et S
    ! -------                                                  --------    .
    implicit none
    type(date_secPS), intent (in) :: thedate1,thedate2
    real(KIND=wr), intent (out) :: deltaP,deltaS
    real(KIND=wr) :: delta
    ! -----------------------------------------------------------------    .
    delta = real(thedate1%date%Jday-thedate2%date%Jday,wr)*24.0_wr*60.0_wr*60.0_wr + &
    real(thedate1%date%hour-thedate2%date%hour,wr)*60.0_wr*60.0_wr + &
    real(thedate1%date%min-thedate2%date%min,wr)*60.0_wr
    ! -------                                                  --------    .
    deltaP = delta + thedate1%secP - thedate2%secP
    deltaS = delta + thedate1%secS - thedate2%secS
    ! -----------------------------------------------------------------    .
  end subroutine difftimePS

! ---------------------------------------------------------------------    .

  subroutine difftimeP(deltatps,thedate1,thedate2)
    ! -------                                                  --------    .mh
    ! difference relative entre deux temps absolus, résultat en sec
    ! -------                                                  --------    .
    implicit none
    type(date_sec), intent (in) :: thedate1,thedate2
    real(KIND=wr), intent (out) :: deltatps
    real(KIND=wr) :: delta
    ! -----------------------------------------------------------------    .
    delta = real(thedate1%date%Jday-thedate2%date%Jday,wr)*24.0_wr*60.0_wr*60.0_wr + &
    real(thedate1%date%hour-thedate2%date%hour,wr)*60.0_wr*60.0_wr + &
    real(thedate1%date%min-thedate2%date%min,wr)*60.0_wr
    ! -------                                                  --------    .
    deltatps = delta + thedate1%sec - thedate2%sec
    ! -----------------------------------------------------------------    .
  end subroutine difftimeP

! ---------------------------------------------------------------------    .

   subroutine tempszero(atime)
    ! -------                                                  --------    .mh
    ! initialise le temps à zéro, pratique parfois ...
    ! -------                                                  --------    .
    implicit none
    type(date_min), intent(out) :: atime
    ! -----------------------------------------------------------------    .
    atime%Jday=0
    atime%year=0
    atime%month=0
    atime%day=0
    atime%hour=0
    atime%min=0
    ! -----------------------------------------------------------------    .
  end subroutine tempszero

END MODULE time



! *********************************************************************    .
! *********************************************************************    .


