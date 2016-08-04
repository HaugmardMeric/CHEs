! Librairie de subroutine concernant le McMC : algorithme de Métropolis-Hastings 
! septembre 2013
! *********************************************************************    .
! ------- Méric Haugmard meric.haugmard@univ-nantes.fr         --------    .
! *********************************************************************    .
! ---------------------------------------------------------------------    .

MODULE algo_metropolis

    use modparam
    use typetemps, only : parametresinv,parametres,fcout,accept
    use mt19937, only : genrand_real1

    implicit none

    private

    public  :: metropolis


CONTAINS

! ---------------------------------------------------------------------    .

  subroutine metropolis(p,param_best,misfit,acceptance,mod_save,newkeep)
    ! -------                                                  --------    .mh
    ! permet l'acceptation et le rejet des modèles selon un critère d'acceptation
    ! -------                                                  --------    .
    implicit none
    type(parametresinv),intent(inout) :: p
    type(parametres),intent(inout) :: param_best
    type(fcout),intent(inout) :: misfit
    type(accept),intent(inout) :: acceptance
    logical,intent(in) :: mod_save
    logical, intent(out) :: newkeep
    ! -------                                                  --------    .
    real(KIND=wr) :: aleatoire,h,delta
    ! -----------------------------------------------------------------    .
    h=-700._wr                                                             ! exp(h)< 10^307 (_wr)
    ! -------                                                  --------    .
    newkeep = .true.
    if (misfit%new.gt.misfit%old) then                                     ! pire
      aleatoire = genrand_real1()
      delta=max(misfit%old-misfit%new,h)                                   ! pb exec sinon car déjà : exp(-700)=9.9e-305
      if (exp(delta).le.aleatoire) then                                    ! repêchage ?
        newkeep = .false.                                                  ! repêchage : non
        misfit%new = misfit%old
        p%valNew=p%valOld
        acceptance%N = acceptance%N+int(1,wl)
      else
        acceptance%NO = acceptance%NO+int(1,wl)                            ! repêchage : oui
      endif
    else
      acceptance%O = acceptance%O+int(1,wl)                                ! meilleur
    endif
    !------- if newkeep = .true. -> new parameters             --------    .
    if (newkeep) then
      p%valOld = p%valNew
      misfit%old = misfit%new
    end if
    !------- the lowest misfit is updated                      --------    .
    if ((misfit%new.le.misfit%best).and.(mod_save)) then
      misfit%best=misfit%new
      param_best=p%valNew
    end if
    ! -----------------------------------------------------------------    .
  end subroutine metropolis

END MODULE algo_metropolis



! *********************************************************************    .
! *********************************************************************    .


