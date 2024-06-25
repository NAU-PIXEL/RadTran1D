module module_diag_common

        implicit none

        private

        public :: test_for_bad, test_for_bad_detailed, test_for_bad_detailed_2d

contains
        
subroutine test_for_bad(in_value,field_name, &
                        do_fix,fix_with)

   use, INTRINSIC :: IEEE_ARITHMETIC
   use module_wrf_error

   implicit none

   real               ::           in_value
   real,intent(in   ) ::           fix_with
   character(len=*), intent(in) :: field_name
   logical, intent(in) ::          do_fix
   logical :: no_error
   
   no_error = ieee_is_normal(in_value)

   if(no_error) then
     return
   else
     write(wrf_err_message,*) "NaN or inf in variable " &
                              //trim(field_name)
     call wrf_message(trim(wrf_err_message))

     if(do_fix) in_value = fix_with
   endif

end subroutine test_for_bad


subroutine test_for_bad_detailed(array,field_name, &
                          ims,ime,jms,jme,kms,kme, &
                          its,ite,jts,jte,kts,kte, &
                          do_fix,fix_with,         &
                          orig_int, corrected,     &
                          kdm_mode                 )

   use, INTRINSIC :: IEEE_ARITHMETIC
   use module_wrf_error

   implicit none

   integer, intent(in) :: ims,ime,jms,jme,kms,kme, &
                          its,ite,jts,jte,kts,kte

   character(len=*), intent(in) :: field_name

   real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: array

   logical, intent(in) :: do_fix
   real, intent(in)    :: fix_with

   logical, optional, intent(in)  :: orig_int
   integer, optional, intent(in)  :: kdm_mode
   integer, optional, intent(out), &
                        dimension(ims:ime,kms:kme,jms:jme) :: corrected

   logical :: is_error

   integer :: i,j,k

   if(present(corrected) .and. (.not.present(orig_int))) &
     call wrf_error_fatal("test_for_bad: got corrected but no orig_int flag")

   if(present(orig_int)) then  ! load integer array
     if(present(corrected)) then
       if(orig_int) corrected(its:ite,kts:kte,jts:jte)=  &
                    int(array(its:ite,kts:kte,jts:jte))
     else
       call wrf_error_fatal("test_for_bad orig_int=.t., but no corrected array passed")
     endif
   endif

   do j=jts,jte
    do k=kts,kte
     do i=its,ite

      is_error = .not.ieee_is_normal(array(i,k,j))

      if(is_error) then

        if(present(kdm_mode)) call kdm_mode_message(kdm_mode)

        write(wrf_err_message,*) "NaN or inf in field " &
                                 //trim(field_name)//   &
                                 " at i,j,k=",i,j,k
        call wrf_message(trim(wrf_err_message))

         if(do_fix) call fix_the_nan(array(i,k,j),fix_with,corrected(i,k,j),orig_int)

      endif

     enddo
    enddo
   enddo

end subroutine test_for_bad_detailed

subroutine test_for_bad_detailed_2d(array,field_name,  &
                                    ims,ime,jms,jme,   &
                                    its,ite,jts,jte,   &
                                    do_fix, fix_with,  &
                                    orig_int,corrected,&
                                    kdm_mode )

   use, INTRINSIC :: IEEE_ARITHMETIC
   use module_wrf_error

   implicit none

   integer, intent(in) :: ims,ime,jms,jme, &
                          its,ite,jts,jte

   character(len=*), intent(in) :: field_name

   real, dimension(ims:ime,jms:jme), intent(inout) :: array

   logical, intent(in) :: do_fix
   real, intent(in)    :: fix_with

   logical, optional, intent(in)  :: orig_int
   integer, optional, intent(in)  :: kdm_mode
   integer, optional, intent(inout), &
                        dimension(ims:ime,jms:jme) :: corrected

   logical :: no_error

   integer :: i,j

   if(present(corrected) .and. (.not.present(orig_int))) &
     call wrf_error_fatal("test_for_bad: got corrected but no orig_int flag")

   if(present(orig_int) .and. (.not.present(corrected))) &
     call wrf_error_fatal("test_for_bad orig_int=.t., but no corrected array passed")

   do j=jts,jte
     do i=its,ite

      no_error = ieee_is_normal(array(i,j))

      if(.not.no_error) then

        if(present(kdm_mode)) call kdm_mode_message(kdm_mode)

        write(wrf_err_message,*) "NaN or inf in field " &
                                 //trim(field_name)//   &
                                 " at i,j=",i,j
        call wrf_message(trim(wrf_err_message))

        if(do_fix) call fix_the_nan(array(i,j),fix_with,corrected(i,j),orig_int)

      endif

     enddo
   enddo

end subroutine test_for_bad_detailed_2d

subroutine fix_the_nan(array,fix_with,corrected,orig_int)

          implicit none

          real, intent(inout) :: array
          real, intent(in   ) :: fix_with
          integer, intent(inout), optional :: corrected
          logical, intent(in   ), optional :: orig_int

            if(present(orig_int)) then
              if(orig_int) then 
                corrected=int(fix_with)
              else
                array = fix_with
              endif
            else
              array = fix_with  ! replace a NaN with value
            endif

end subroutine fix_the_nan

subroutine kdm_mode_message(kdm_mode)

          use module_wrf_error

          implicit none

          integer, intent(in) :: kdm_mode

          if(kdm_mode==0) then      ! kdm module variable kdm_do_vis
             call wrf_message("KDM visible mode")
          elseif(kdm_mode==1) then  ! kdm module variable kdm_do_ir
             call wrf_message("KDM infrared mode")
          else
             call wrf_error_fatal("kdm_mode specified but not valid")
          endif

end subroutine kdm_mode_message

end module module_diag_common

