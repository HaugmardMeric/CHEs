! module indiquant l'avancement d'une boucle avec une barre de progression
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE cpt_temps

    use modparam

    implicit none

    private

    public  :: progress


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine progress(n,ntotal,one_string)
    ! -------                                                  --------    .mh
    implicit none
    integer(KIND=wi), intent(in) :: n                                      ! itération en cours
    integer(KIND=wi), intent(in) :: ntotal                                 ! itérations totales
    character(LEN=*), intent(in) :: one_string                             ! chaîne à afficher avant la barre [##### 100 % #####]
    ! -------                                                  --------    .
    character(LEN=255) :: prog,oldprog
    real(KIND=wr) :: tl
    integer(KIND=wi) :: i,oldtime,newtime,ratetime,cptmax
    ! -------                                                  --------    .
    save :: oldprog                                                        ! garde en mémoire l'ancienne chaîne et le temps initial
    save :: oldtime
    ! -----------------------------------------------------------------    .
    if (n.eq.1) then
      write(oldprog,'(a255)')"x"                                           ! chaîne vide, première itération
      call system_clock(oldtime)                                           ! temps à la première itération
    endif
    call system_clock(newtime,ratetime,cptmax)                             ! temps à l'itération courante
    ! -------                                                  --------    . différence de temps
    tl = real(newtime-oldtime,wr)
    if(tl.lt.0.0_wr) tl = real(newtime-oldtime+cptmax,wr)
    tl = tl/real(ratetime,wr)
    ! -------                                                  --------    .
    if ((n.gt.0).and.(n.lt.ntotal)) then
      tl=(1.0_wr*real(ntotal,wr)/real(n,wr))*tl-tl                         ! temps restant en sec
    else
      tl=tl                                                                ! temps total en sec
    endif
    if (tl.le.0.0_wr) tl=0.0_wr                                            ! première(s) itération(s)
    ! -------                                                  --------    .
    write(prog,'(a30,1x,''['')')one_string//"   "                          ! construction de la chaîne de caractère
    do i=1,40
      prog(32+i:32+i)=' '
    enddo
    write(prog(48:56),'(f7.1,''%'')') 100.0_wr*real(n,wr)/real(ntotal,wr)  ! pourcentage restant
    do i=1,40
      if ((1.0_wr*real(n,wr)/real(ntotal,wr)).gt.(1.0_wr*real(i,wr)/40.0_wr)) then
        if (prog(32+i:32+i).eq.' ') prog(32+i:32+i)='#'                    ! barre de progression ##
      endif
    enddo
    prog(72:72)=']'
    write(prog(75:77),'(i2.2,'':'')')int(tl/3600.0_wr)                     ! temps restant en heure

    write(prog(78:80),'(i2.2,'':'')')int((tl-real(int(tl/3600.0_wr,wi),wr)*3600.0_wr)/60.0_wr,wi) ! temps restant en min

    write(prog(81:82),'(i2.2)')int((tl-real(int(tl/60.0_wr,wi),wr)*60.0_wr),wi) ! temps restant en sec
    write(prog(83:86),'(a4)')" sec"
    if (prog.ne.oldprog) write(*,'(a,a,$)') prog(1:86),char(13)            ! permet de rester sur la ligne courante et n'écrit que si nécessaire
    oldprog=prog
    if (n.eq.ntotal) write(*,*)                                            ! fin,char(13)
    ! -----------------------------------------------------------------    .
  end subroutine progress

END MODULE cpt_temps



! *********************************************************************    .
! *********************************************************************    .


