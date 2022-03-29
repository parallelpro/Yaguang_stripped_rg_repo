! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none

      ! these routines are called by the standard run_star check_model
      contains
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going   

        
        !  ! protostar to rgb numax=20
        !  if ((s% center_h1 < 0.0001) .and. (s%nu_max <= 20)) then
        !     !(s%delta_Pg <= (0.13 * s%delta_nu / (s%nu_max**2.0) * 1.e6)) .and. (s%delta_Pg .ne. 0.) ) then
        !     ! stop when star hydrogen mass drops to specified level
        !     extras_check_model = terminate
        !     termination_code_str(t_xtra1) = 'have reached RGB numax = 2 muHz.'
        !     return
        !  end if

        !  ! rgb numax=20 to HeZAMS
        !  if ((s%power_he_burn/s%power_nuc_burn)>0.9) .and. ((s%center_he4+s%center_he3) <0.95) ) then !(phase_helium_burning >=5)
        !     extras_check_model = terminate
        !     termination_code_str(t_xtra1) = 'have reached cHeB.'
        !     return
        !  end if

        !  ! HeZAMS to AGB
        !  if ((s%power_he_burn/s%power_nuc_burn)>0.9) .and. (s%center_he4+s%center_he3) < 0.001) then
        !     extras_check_model = terminate
        !     termination_code_str(t_xtra1) = 'have reached AGB.'
        !     return
        !  end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         
         names(1) = 'helium_ignition'
         vals(1) = merge(1.d0, 0.d0, s%helium_ignition) 
         names(2) = 'phase_of_evolution'
         vals(2) = s%phase_of_evolution

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! record variables that is going to need for saving profiles
         s%xtra(1) = s%delta_nu
         s%xtra(2) = s%Teff
         s%xtra(3) = log10(s%L_phot)

        !  ! ZAMS to RGB tip (Dnu = 2 muHz)
        !  if (((s%center_h1 < 0.0001) .or. ((s%L_nuc_burn_total/s%L_surf)>=0.99) ) &
        !         .and. (s%nu_max>=20.) .and. (.not.(s%power_he_burn/s%power_nuc_burn)>0.9))) then
        !     ! save profile models
        !     s% write_profiles_flag = .true. 

        !     ! use a custom frequency, save 3 profiles per 1muHz in Dnu, 
        !     ! or 0.2 profiles per 1K in Teff, or 100 profiles per 1dex in log10(L/Lsun)
        !     ! d(N_m)/d(Dnu) = d(N_m)/d(N_p) * d(N_p)/d(Dnu) = profile_interval * 2/1muHz
        !     ! d(N_m)/d(Teff) = d(N_m)/d(N_p) * d(N_p)/d(Teff) = profile_interval * 0.1/1K
        !     ! d(N_m)/d(log10(L)) = d(N_m)/d(N_p) * d(N_p)/d(log10(L)) = profile_interval * 100/1dex
        !     s% profile_interval = max( min(int(1./abs(s%xtra(1)-s%xtra_old(1))/2.), &
        !                                    int(1./abs(s%xtra(2)-s%xtra_old(2))/0.1), &
        !                                    int(1./abs(s%xtra(3)-s%xtra_old(3))/100)) , 1)

        !     ! change timestep control
        !     s% delta_lgTeff_hard_limit = 0.0005
        !     s% delta_lgL_limit_L_min = 0.0
        !     s% delta_lgL_hard_limit = 0.004
        !  endif

        !  ! RGB tip (Dnu = 2 muHz) to He-ZAMS
        !  ! to speed things up
        !  if ((s% center_h1 < 0.0001) .and. (s%nu_max<20.) .and. (.not.(s%power_he_burn/s%power_nuc_burn)>0.9))) then
        !     ! save profile models
        !     s% write_profiles_flag = .true.

        !     ! use a custom frequency, save 3 profiles per 1muHz in Dnu, 
        !     ! or 0.2 profiles per 1K in Teff, or 100 profiles per 1dex in log10(L/Lsun)
        !     ! d(N_m)/d(Dnu) = d(N_m)/d(N_p) * d(N_p)/d(Dnu) = profile_interval * 2/1muHz
        !     ! d(N_m)/d(Teff) = d(N_m)/d(N_p) * d(N_p)/d(Teff) = profile_interval * 0.1/1K
        !     ! d(N_m)/d(log10(L)) = d(N_m)/d(N_p) * d(N_p)/d(log10(L)) = profile_interval * 100/1dex
        !     s% profile_interval = max( min(int(1./abs(s%xtra(1)-s%xtra_old(1))/2.), &
        !                                    int(1./abs(s%xtra(2)-s%xtra_old(2))/0.1), &
        !                                    int(1./abs(s%xtra(3)-s%xtra_old(3))/100)) , 1)

        !     ! change timestep control
        !     s% delta_lgTeff_hard_limit = -1
        !     s% delta_lgL_limit_L_min = -100
        !     s% delta_lgL_hard_limit = -1
        !     !s% varcontrol_target = 0.0002
        !     ! mesh controls 
        !     ! s% mesh_delta_coeff = 1.0
        !  endif       

         ! He-ZAMS to AGB
         if (((s%L_by_category(3)/s%L_nuc_burn_total)>0.1) .and. ((s%center_he4+s%center_he3) <0.92) &
                .and.((s%center_he4+s%center_he3) >0.01) .and. (log10(s%star_mdot)<-30)) then
            ! save profile models
            s% write_profiles_flag = .true. 

            ! use a custom frequency, save 3 profiles per 1muHz in Dnu, 
            ! or 0.2 profiles per 1K in Teff, or 100 profiles per 1dex in log10(L/Lsun)
            ! d(N_m)/d(Dnu) = d(N_m)/d(N_p) * d(N_p)/d(Dnu) = profile_interval * 2/1muHz
            ! d(N_m)/d(Teff) = d(N_m)/d(N_p) * d(N_p)/d(Teff) = profile_interval * 0.1/1K
            ! d(N_m)/d(log10(L)) = d(N_m)/d(N_p) * d(N_p)/d(log10(L)) = profile_interval * 100/1dex
            s% profile_interval = max( min(int(1./abs(s%xtra(1)-s%xtra_old(1))/2.), &
                                           int(1./abs(s%xtra(2)-s%xtra_old(2))/0.1), &
                                           int(1./abs(s%xtra(3)-s%xtra_old(3))/100)) , 1)

            ! change timestep control
            s% delta_lgTeff_hard_limit = 0.0005
            s% delta_lgL_limit_L_min = 0.0
            s% delta_lgL_hard_limit = 0.004
            !s% varcontrol_target = 0.0001
            ! s% delta_XHe_cntr_limit = 0.01

            ! mesh controls 
            ! s% mesh_delta_coeff = 1.0
         endif

         !see extras_check_model for information about custom termination codes
         !by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step

      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      
      end module run_star_extras
