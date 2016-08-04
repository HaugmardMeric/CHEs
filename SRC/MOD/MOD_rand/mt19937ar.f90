! ---------------------------------------------------------------------    .
! -------  modification début 27/03/2014                       --------    .1_mh
! ---------------------------------------------------------------------    .
!   compilation :
!       mt19937ar.o : MOD/MOD_rand/mt19937ar.f90 modparam.o
!       gfortran -c -fno-range-check $(OPTIONC) $(FFLAGS) $<
!                                                                          
! ---------------------------------------------------------------------    .
! Le statisticien William Youden (en) écrit69 en 1962 une explication      !
! du but et de la position de la loi normale dans les sciences.            !
! Il la présente sous forme de courbe en cloche :                          !
!                                                                          !
!                                   THE                                    !
!                                  NORMAL                                  !
!                               LAW OF ERROR                               !
!                            STANDS OUT IN THE                             !
!                          EXPERIENCE OF MANKIND                           !
!                         AS ONE OF THE BROADEST                           !
!                        GENERALIZATIONS OF NATURAL                        !
!                      PHILOSOPHY . IT SERVES AS THE                       !
!                     GUIDING INSTRUMENT IN RESEARCHES                     !
!                 IN THE PHYSICAL AND SOCIAL SCIENCES AND                  !
!                IN MEDICINE AGRICULTURE AND ENGINEERING .                 !
!           IT IS AN INDISPENSABLE TOOL FOR THE ANALYSIS AND THE           !
! INTERPRETATION OF THE BASIC DATA OBTAINED BY OBSERVATION AND EXPERIMENT  !
!                                                                          !
! ---------------------------------------------------------------------    .
! -------  modification fin                                    --------    .1_mh
! ---------------------------------------------------------------------    .

! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.

! Code converted to Fortran 95 by Josi Rui Faustino de Sousa
! Date: 2002-02-01

! Before using, initialize the state by using init_genrand(seed)
! or init_by_array(init_key, key_length).

! This library is free software.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
! Any feedback is very welcome.
! http://www.math.keio.ac.jp/matumoto/emt.html
! email: matumoto@math.keio.ac.jp

      MODULE mt19937

    ! ---------------------------------------------------------------------    .
    ! -------  modification début 27/03/2014                       --------    . 1_mh
    ! ---------------------------------------------------------------------    .
        use modparam
    ! ---------------------------------------------------------------------    .
    ! -------  modification fin                                    --------    . 1_mh
    ! ---------------------------------------------------------------------    .

        implicit none

        intrinsic :: bit_size

        private

    ! -------------------------------------------------------------------    .
    ! -------  modification début 27/03/2014                     --------    . 2_mh
    ! -------------------------------------------------------------------    .
        public  :: initseed
    ! -------------------------------------------------------------------    .
    ! -------  modification fin                                  --------    . 2_mh
    ! -------------------------------------------------------------------    .

        public  :: init_genrand, init_by_array
        public  :: genrand_int32, genrand_int31
        public  :: genrand_real1, genrand_real2, genrand_real3, genrand_res53
        public  :: normal

    ! -------------------------------------------------------------------    .
    ! -------  modification début 27/03/2014                     --------    . 3_mh
    ! -------------------------------------------------------------------    .
    !   integer,  parameter  :: intg = selected_int_kind( 9 )
    !   integer,  parameter  :: long = selected_int_kind( 18 )
    !   integer,  parameter  :: flot = selected_real_kind( 6, 37 )
    !   integer,  parameter  :: dobl = selected_real_kind( 15, 307 )
    !   integer,  public, parameter :: wi = intg
    !   integer,  public, parameter :: wl = long
    !   integer,  public, parameter :: wr = dobl
    ! -------------------------------------------------------------------    .
    ! -------  modification fin                                  --------    . 3_mh
    ! -------------------------------------------------------------------    .

    ! Period parameters

        integer( kind = wi ), parameter :: n = 624_wi
        integer( kind = wi ), parameter :: m = 397_wi
        integer( kind = wi ), parameter :: hbs = bit_size( n ) / 2_wi
        integer( kind = wi ), parameter :: qbs = hbs / 2_wi
        integer( kind = wi ), parameter :: tbs = 3_wi * qbs

        integer( kind = wi )  :: mt(n)                                       ! the array for the state vector
        logical( kind = wi )  :: mtinit = .false._wi                         ! means mt[N] is not initialized
        integer( kind = wi )  :: mti = n + 1_wi                              ! mti==N+1 means mt[N] is not initialized

      contains

  ! ---------------------------------------------------------------------    .
  ! -------  modification début 27/03/2014                       --------    . 4_mh
  ! ---------------------------------------------------------------------    .
  subroutine initseed(lib)
    ! -------                                                    --------    . mh
    ! initialisation de la graine du générateur de nb aleatoire
    ! fixe (lib=0) ou non (lib=1)
    ! -------                                                    --------    .
    implicit none
    ! -------                                                    --------    .
    logical( kind = wi ), intent (in)  :: lib
    integer( kind = wi )  :: i,seed,countn,count_rate,count_max
    ! -------------------------------------------------------------------    .
    if(lib) then
      call system_clock(countn,count_rate,count_max)
      do while (countn > 10000000)
        i = int(10000123*genrand_real1())
        countn = countn - i
      end do
      if (countn < 0) then
        write(*,*)'probleme dans initseed : mauvaise génération de graine'
        stop
      endif
      seed = countn
    else
      seed = 12345678
    endif
    call init_genrand (seed)
    write(*,*)'générateur de nombre aléatoire :       ',seed
    ! -------------------------------------------------------------------    .
  end subroutine initseed
  ! ---------------------------------------------------------------------    .
  ! -------  modification fin                                    --------    . 4_mh
  ! ---------------------------------------------------------------------    .

        elemental function uiadd( a, b ) result( c )

          implicit none
          intrinsic :: ibits, ior, ishft
          integer( kind = wi ), intent( in )  :: a, b
          integer( kind = wi )  :: c
          integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

          a1 = ibits( a, 0, hbs )
          a2 = ibits( a, hbs, hbs )
          b1 = ibits( b, 0, hbs )
          b2 = ibits( b, hbs, hbs )
          s1 = a1 + b1
          s2 = a2 + b2 + ibits( s1, hbs, hbs )
          c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
          return
        end function uiadd

        elemental function uisub( a, b ) result( c )

          implicit none
          intrinsic :: ibits, ior, ishft
          integer( kind = wi ), intent( in )  :: a, b
          integer( kind = wi )  :: c
          integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

          a1 = ibits( a, 0, hbs )
          a2 = ibits( a, hbs, hbs )
          b1 = ibits( b, 0, hbs )
          b2 = ibits( b, hbs, hbs )
          s1 = a1 - b1
          s2 = a2 - b2 + ibits( s1, hbs, hbs )
          c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
          return
        end function uisub

        elemental function uimlt( a, b ) result( c )
    
          implicit none
          intrinsic :: ibits, ior, ishft
          integer( kind = wi ), intent( in )  :: a, b
          integer( kind = wi )  :: c
          integer( kind = wi )  :: a0, a1, a2, a3
          integer( kind = wi )  :: b0, b1, b2, b3
          integer( kind = wi )  :: p0, p1, p2, p3

          a0 = ibits( a, 0, qbs )
          a1 = ibits( a, qbs, qbs )
          a2 = ibits( a, hbs, qbs )
          a3 = ibits( a, tbs, qbs )
          b0 = ibits( b, 0, qbs )
          b1 = ibits( b, qbs, qbs )
          b2 = ibits( b, hbs, qbs )
          b3 = ibits( b, tbs, qbs )
          p0 = a0 * b0
          p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
          p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
          p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
          c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
          c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
          c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )
          return
        end function uimlt

  ! initializes mt[N] with a seed
        subroutine init_genrand( s )

          implicit none
          intrinsic :: iand, ishft, ieor, ibits
          integer( kind = wi ), intent( in )  :: s
          integer( kind = wi )  :: i, mult_a

          data mult_a /z'6C078965'/

          mtinit = .true._wi
          mt(1) = ibits( s, 0, 32 )
          do i = 2, n, 1
             mt(i) = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
             mt(i) = uimlt( mt(i), mult_a )
             mt(i) = uiadd( mt(i), uisub( i, 1_wi ) )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
      ! In the previous versions, MSBs of the seed affect
      ! only MSBs of the array mt[].
      ! 2002/01/09 modified by Makoto Matsumoto
             mt(i) = ibits( mt(i), 0, 32 )
      ! for >32 bit machines
          end do
          return
        end subroutine init_genrand

  ! initialize by an array with array-length
  ! init_key is the array for initializing keys
  ! key_length is its length
        subroutine init_by_array( init_key )

          implicit none
          intrinsic :: iand, ishft, ieor
          integer( kind = wi ), intent( in )  :: init_key(:)
          integer( kind = wi )  :: i, j, k, tp, key_length
          integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

          data seed_d /z'12BD6AA'/
          data mult_a /z'19660D'/
          data mult_b /z'5D588B65'/
          data msb1_d /z'80000000'/

          key_length = size( init_key, dim = 1 )
          call init_genrand( seed_d )
          i = 2_wi
          j = 1_wi
          do k = max( n, key_length ), 1, -1
             tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
             tp = uimlt( tp, mult_a )
             mt(i) = ieor( mt(i), tp )
             mt(i) = uiadd( mt(i), uiadd( init_key(j), uisub( j, 1_wi ) ) ) ! non linear
             mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
             i = i + 1_wi
             j = j + 1_wi
             if ( i > n ) then
                mt(1) = mt(n)
                i = 2_wi
             end if
             if ( j > key_length) j = 1_wi
          end do
          do k = n-1, 1, -1
             tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
             tp = uimlt( tp, mult_b )
             mt(i) = ieor( mt(i), tp )
             mt(i) = uisub( mt(i), uisub( i, 1_wi ) ) ! non linear
             mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
             i = i + 1_wi
             if ( i > n ) then
                mt(1) = mt(n)
                i = 2_wi
             end if
          end do
          mt(1) = msb1_d ! MSB is 1; assuring non-zero initial array
          return
        end subroutine init_by_array

  ! generates a random number on [0,0xffffffff]-interval
        function genrand_int32( ) result( y )

          implicit none
          intrinsic :: iand, ishft, ior, ieor, btest, ibset, mvbits
          integer( kind = wi )  :: y
          integer( kind = wi )  :: kk
          integer( kind = wi )  :: seed_d, matrix_a, matrix_b, temper_a,temper_b

          data seed_d   /z'5489'/
          data matrix_a /z'9908B0DF'/
          data matrix_b /z'0'/
          data temper_a /z'9D2C5680'/
          data temper_b /z'EFC60000'/

          if ( mti > n ) then ! generate N words at one time
             if ( .not. mtinit ) call init_genrand( seed_d ) ! if init_genrand() has not been called, a default initial seed is used
             do kk = 1, n-m, 1
                y = ibits( mt(kk+1), 0, 31 )
                call mvbits( mt(kk), 31, 1, y, 31 )
                if ( btest( y, 0 ) ) then
                   mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_a )
                else
                   mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_b )
                end if
             end do
             do kk = n-m+1, n-1, 1
                y = ibits( mt(kk+1), 0, 31 )
                call mvbits( mt(kk), 31, 1, y, 31 )
                if ( btest( y, 0 ) ) then
                   mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_a )
                else
                   mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_b )
                end if
             end do
             y = ibits( mt(1), 0, 31 )
             call mvbits( mt(n), 31, 1, y, 31 )
             if ( btest( y, 0 ) ) then
                mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_a )
             else
                mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_b )
             end if
             mti = 1_wi
          end if
          y = mt(mti)
          mti = mti + 1_wi
          ! Tempering
          y = ieor( y, ishft( y, -11) )
          y = ieor( y, iand( ishft( y, 7 ), temper_a ) )
          y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
          y = ieor( y, ishft( y, -18 ) )
          return
        end function genrand_int32

  ! generates a random number on [0,0x7fffffff]-interval
        function genrand_int31( ) result( i )

          implicit none
          intrinsic :: ishft
          integer( kind = wi )  :: i
          i = ishft( genrand_int32( ), -1 )
          return
        end function genrand_int31

  ! generates a random number on [0,1]-real-interval
        function genrand_real1( ) result( r )

          implicit none
          real( kind = wr )  :: r
          integer( kind = wi )  :: a, a1, a0

          a = genrand_int32( )
          a0 = ibits( a, 0, hbs )
          a1 = ibits( a, hbs, hbs )
           r = real( a0, kind = wr ) * 2.3283064370807973754315e-10_wr
!          r = real( a0, kind = wr ) / 4294967295.0_wr
           r = real( a1, kind = wr ) * 1.525878906605271367963e-5_wr    + r
!          r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
    ! divided by 2^32-1
          return
        end function genrand_real1

  ! generates a random number on [0,1)-real-interval
        function genrand_real2( ) result( r )

          implicit none
          intrinsic :: ibits
          real( kind = wr )  :: r
          integer( kind = wi )  :: a, a1, a0

          a = genrand_int32( )
          a0 = ibits( a, 0, hbs )
          a1 = ibits( a, hbs, hbs )
          r = real( a0, kind = wr ) * 2.3283064365386962890625e-10_wr
!          r = real( a0, kind = wr ) / 4294967296.0_wr
          r = real( a1, kind = wr ) * 1.52587890625e-5_wr + r
!          r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32
          return
        end function genrand_real2

  ! generates a random number on (0,1)-real-interval
        function genrand_real3( ) result( r )

          implicit none
          real( kind = wr )  :: r
          integer( kind = wi )  :: a, a1, a0

          a = genrand_int32( )
          a0 = ibits( a, 0, hbs )
          a1 = ibits( a, hbs, hbs )
          r = ( real( a0, kind = wr ) + 0.5_wr ) * 2.3283064365386962890625e-10_wr
!          r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
          r = real( a1, kind = wr ) * 1.52587890625e-5_wr + r
!          r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32
          return
        end function genrand_real3

  ! generates a random number on [0,1) with 53-bit resolution
        function genrand_res53( )  result( r )

          implicit none
          intrinsic :: ishft
          real( kind = wr )  :: r
          integer( kind = wi )  :: a, a0, a1
          integer( kind = wi )  :: b, b0, b1

          a = ishft( genrand_int32( ), -5 )
          a0 = ibits( a, 0, hbs )
          a1 = ibits( a, hbs, hbs )
          b = ishft( genrand_int32( ), -6 )
          b0 = ibits( b, 0, hbs )
          b1 = ibits( b, hbs, hbs )
          r = real( a1, kind = wr ) * 4.8828125e-4_wr
!          r = real( a1, kind = wr ) / 2048.0_wr
          r = real( a0, kind = wr ) * 7.450580596923828125e-9_wr + r
!          r = real( a0, kind = wr ) / 134217728.0_wr + r
          r = real( b1, kind = wr ) * 7.27595761418342590332e-12_wr + r
!          r = real( b1, kind = wr ) / 137438953472.0_wr + r
          r = real( b0, kind = wr ) * 1.1102230246251565404236e-16_wr + r
!          r = real( b0, kind = wr ) / 9007199254740992.0_wr + r
          return
        end function genrand_res53
  ! These real versions are due to Isaku Wada, 2002/01/09 added

  ! tirage au sort selon une loi normale

       function normal0()

          implicit none
          real(kind=wr)::r1,r2,rsq,normal0
          real(wr),save::g
          logical,save::gaus_stored=.false.
              if (gaus_stored) then
                      normal0=g
                      gaus_stored=.false.
              else
              do
                      r1=2._wr*genrand_real1()-1._wr
                      r2=2._wr*genrand_real1()-1._wr
                      rsq=r1**2+r2**2
                      if (rsq >0.0_wr .and. rsq < 1.0) exit
                      end do
                      rsq=sqrt(-2.0_wr*log(rsq)/rsq)
                      normal0=r1*rsq
                      g=r2*rsq
                      gaus_stored=.true.
              end if
          return          
        end function normal0

       function normal(m,s) 

          implicit none
          real(kind=wr):: normal
          real(wr),intent(in)::m,s
              normal=s*normal0()+m
          return  
        end function normal

      END MODULE mt19937