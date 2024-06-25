MODULE module_read_soundings

#ifdef mpas
      use mpas_constants , only : g => gravity, &
                                r_d => rgas
#else
      use module_model_constants, only: g, r_d
#endif
      use module_m_time, only: d2j
      use module_mars24, only: mars24_Mars_Solar_Date, &
                               mars24_j2000_offset_tt, &
                               mars24_Mars_year
      use module_wrf_error
      use module_nrutils, only : sort2            ! share

      implicit none

      real, parameter :: R_air = r_d

      PRIVATE

      PUBLIC :: READ_SOUNDINGS, READ_CARMA_PROFILES

      character(len=25) :: carma_data_path='./Data/titan/CARMA/Titan/'

contains

      subroutine read_soundings( rs_file, info_type,                &
             t_sounding,                                            &
             sound_lat, sound_lon, sound_elev, sound_psfc,          &
             sound_ls,  sound_time, sound_doy, sound_msd,           &
             z_profile_out, t_profile_out, n_profile )

      implicit none

      character(len=50), intent(in) :: rs_file
      integer, intent(in)           :: info_type, t_sounding

      ! t_sounding is passed because some profiles setup here are T(p) and
      ! some are t(z).  Specifically:
      !  t_sounding = 42  =>  T(z)
      !  t_sounding = 43  =>  T(p)

      REAL, INTENT(OUT), OPTIONAL   :: sound_lat, sound_lon
      REAL, INTENT(OUT), OPTIONAL   :: sound_elev, sound_psfc
      REAL, INTENT(OUT), OPTIONAL   :: sound_ls,  sound_time 
      REAL, INTENT(OUT), OPTIONAL   :: sound_doy, sound_msd

      logical :: end_of_file

      real :: c, p_lowest, t_lowest 
      real :: psfc, elev_lowest, ground_radius
      real, dimension(1000) :: dist, temp
      real, intent(out), optional, dimension(1000) :: z_profile_out, &
                                                     t_profile_out
      integer, intent(out), optional :: n_profile

      real(kind(0.d0)) :: calc_jdate, j_o, msd, mars_year, day_of_year
      integer, dimension(8) :: dat

      character(len=110) :: a,b,e,date_string

      character(len=100) :: header_file, data_file, path

      integer :: k,nl, i_year, i_month, i_day, calc_err
      integer :: i_hour, i_min, i_sec

      integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9
      real    :: r1, r2, r3, r4, r5, r6, r7, r8, r9
      real    :: r10, r11, r12, r13, r14, r15, r16

      logical :: passed_thru_but_found_no_file = .true.

      ! MOLA Radius = radius of surface
      ! MOLA Areoid = radius of the 0km reference surface 
      ! radius      = radius of the lowest measurement point
      !  Note: - the MOLA elevation of the surface can be found 
      !          from the "MOLA Radius" minus "MOLA Areoid"
      !        - the MOLA elevation of any measurement point
      !          can be found from the measurement point radius
      !          minus the "MOLA Areoid"


      ! Each of the following major blocks is for a different style of input
      ! profile.  The test what what type of profile is on the format of the
      ! filename, so this cannot be used with arbitrary filenames - but since
      ! the assuptions about what is contained in any given style of file is
      ! specific, this routine is not generic anyways.


! -- Profile from GCM output built with python script WRFV3/python_code/get_profile.py --

      if(rs_file(1:14) .eq. 'init_t_profile') then

          if(t_sounding .ne. 43) then
               write(wrf_err_message,*) 'Inconsistent t_sounding type. ', &
               'Set t_sounding=43 if you really want to use a GCM profile to init.'
               CALL wrf_error_fatal(TRIM( wrf_err_message ))
          endif
 
          passed_thru_but_found_no_file = .false.

          data_file   = rs_file
          write(0,*) 'Reading initial souding from file on the local path: ',data_file

          end_of_file = .false.
          k = 0

          open(unit=10,file=data_file,form='formatted',status='old')
          rewind(10)

          do while (.not. end_of_file)

            k = k+1
            if(k.gt.100) then
                    write(wrf_err_message,*) &
                    'woops, GCM file has more points than expected!'
                    CALL wrf_error_fatal(TRIM( wrf_err_message ))
            endif
            if(k .le. 2) then  ! just skip the header
                    read(10,*) a
            else
                    ! don't read anything into layer (1) as we're
                    ! going to fake that layer, plus layer nl+1
                    read(10,*,end=20) dist(k-1),temp(k-1) 
            endif
            go to 21
 20         end_of_file = .true.
 21         continue
          enddo

          close(unit=10,status = 'keep')

          nl = k-1  ! add one fake layer atop and one below

          write(0,*) "profile read in is: T(k),  p(k)"
          do k=1,nl
             write(0,*) k,temp(k),dist(k)
          enddo

          ! figure out if profile is surf -> toa or toa -> surf
          if(dist(2) .gt. dist(nl)) then ! p at bottom > p at top
                  dist(1)  = dist(2) * 100.  ! just go to a ridiculously high P
                  temp(1)  = temp(2)         ! and go down isothermally
                  dist(nl) = 1.e-8           ! this is about 300km altitude
                  temp(nl) = temp(nl-1)      ! also isothermal
          else
                  dist(1)  = 1.e-8
                  temp(1)  = temp(2)
                  dist(nl) = dist(nl-1)*100.
                  temp(nl) = temp(nl-1)
          endif

          write(0,*) "profile read in is: T(k),  p(k)"
          do k=1,nl
             write(0,*) k,temp(k),dist(k)
          enddo

          if(present(z_profile_out)) z_profile_out = dist
          if(present(t_profile_out)) t_profile_out = temp
          if(present(n_profile))     n_profile     = nl

      endif 

! -- Mars Express Radio Science - mex-m-mrs-5-occ-v2.0 ---

      if((rs_file(12:13) .eq. 'ai') .or. (rs_file(12:13) .eq. 'AI')) then

        if(t_sounding .ne. 42) then
               write(wrf_err_message,*) 'Inconsistent t_sounding type. ', &
               'Set t_sounding=42 if you really want to use a RS profile to init.'
               CALL wrf_error_fatal(TRIM( wrf_err_message ))
        endif

        passed_thru_but_found_no_file = .false.
        
        path = "./Data/soundings/mex/"

        header_file = rs_file(:13)//"o"//rs_file(15:27)//".txt"
        data_file   = rs_file(:13)//"x"//rs_file(15:27)//".tab" 

        header_file = trim(path)//header_file
        data_file   = trim(path)//data_file

        if(info_type .lt. 3) then
          write(wrf_err_message,*) '...processing a Mars Express RS profile'
          CALL wrf_message ( TRIM( wrf_err_message ) )
          write(wrf_err_message,*) header_file
          CALL wrf_message ( TRIM( wrf_err_message ) )
          write(wrf_err_message,*) data_file
          CALL wrf_message ( TRIM( wrf_err_message ) )
        endif

        ! info_type:   1 = read just the header and write into also to wrf_err
        !              2 = read just the profile and write info also to wrf_err
        !              3 = read both but no chatter to wrf_err
 
        if((info_type .eq. 1) .or. (info_type .eq. 3)) then  ! grab information on location
          open(unit=10,file=header_file,form='formatted',status='old')
          rewind(10)

          end_of_file = .false.
          k = 0

          if(info_type .eq. 1) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) ' - RS profile HEADER info -'
             CALL wrf_message ( TRIM( wrf_err_message ) )
          endif

          do while (.not. end_of_file)

            c = 0.
            read(10,'(a53,a35)',end=100) a,b
            k = k+1
            if(k.eq.4) date_string = b
            if(info_type.eq.1.and.k.ge.9.and.k.le.27) then
                    write(wrf_err_message,*) trim(a),trim(b)
                    CALL wrf_message ( TRIM( wrf_err_message ) )
            endif
            if((k.ge.5).and.(k.le.30)) read(b,*) c
            if(present(sound_lat) .and. (k.eq.9))  sound_lat = c
            if(present(sound_lon) .and. (k.eq.10)) sound_lon = c
            if(present(sound_ls)  .and. (k.eq.13)) sound_ls = c
            if(k.eq.14) ground_radius = c
            if(present(sound_elev).and. (k.eq.15)) sound_elev = (ground_radius -c)*1000.
            if(k.eq.15) elev_lowest = c
            if(k.eq.16) elev_lowest = (c - elev_lowest)*1000.
            if(k.eq.18) p_lowest = c
            if(k.eq.20) t_lowest = c
            if(present(sound_time) .and. (k.eq.27)) sound_time = c
            go to 110
 100        end_of_file = .true.
 110        continue
          enddo

          if (info_type .eq. 1) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )

             read(date_string(3:6),*) i_year
             read(date_string(8:9),*) i_month
             read(date_string(11:12),*) i_day
             read(date_string(14:15),*) i_hour
             read(date_string(17:18),*) i_min
             read(date_string(20:21),*) i_sec

             write(wrf_err_message,*) date_string
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) i_year,i_month,i_day,i_hour,i_min,i_sec
             CALL wrf_message ( TRIM( wrf_err_message ) )

             dat(1) = i_year
             dat(2) = i_month
             dat(3) = i_day
             dat(4) = 0.     ! this is the offset from UTC
             dat(5) = i_hour
             dat(6) = i_min
             dat(7) = i_sec
             dat(8) = 0.     ! this is the milliseconds - don't care

             call d2j(dat,calc_jdate,calc_err)
             call mars24_j2000_offset_tt(calc_jdate,j_o)
             call mars24_Mars_Solar_Date(j_o,msd)
             call mars24_Mars_Year(j_o,mars_year,day_of_year)
             if(present(sound_doy)) sound_doy=real(day_of_year)
             if(present(sound_msd)) sound_msd=real(msd)

             write(wrf_err_message,*) 'returned jdate=',real(calc_jdate),real(msd)
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) 'mars year and doy=',real(mars_year),real(day_of_year)
             CALL wrf_message ( TRIM( wrf_err_message ) )
          endif

          if(present(sound_psfc)) sound_psfc = p_lowest * exp((elev_lowest-sound_elev)*g/(r_d*t_lowest))

          close(unit=10,status = 'keep')

        endif

        if((info_type .eq. 2) .or. (info_type .eq. 3)) then ! grab the actual thermal profile

          end_of_file = .false.
          k = 0

          open(unit=10,file=data_file,form='formatted',status='old')
          rewind(10)

          do while (.not. end_of_file)


            k = k+1
            if(k.gt.100) then
                    write(wrf_err_message,*) &
                    'woops, RS file has more points than expected!'
                    CALL wrf_error_fatal(TRIM( wrf_err_message ))
            endif
            read(10,'(a50,f8.3,a103,f7.3)',end=200) a,dist(k),b,temp(k)
            dist(k) = 1000.*(dist(k)-ground_radius)
            go to 210
 200        end_of_file = .true.
 210        continue
          enddo

          close(unit=10,status = 'keep')

          nl = k
          dist(nl) = 0.
          temp(nl) = temp(nl-1)

          if(present(z_profile_out)) z_profile_out = dist
          if(present(t_profile_out)) t_profile_out = temp
          if(present(n_profile))     n_profile     = nl

          if(info_type .eq. 2) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) 'Radio Occ profile: (lev, T, z(m))'
             CALL wrf_message ( TRIM( wrf_err_message ) )
             do k = 1, nl
               write(wrf_err_message,*) k,temp(k), dist(k)
               CALL wrf_message ( TRIM( wrf_err_message ) )
             enddo
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
          endif

        endif

      endif

! -- Mars Express Radio Science directly from Dave Hinson ---

      if(rs_file(1:10) .eq. 'hinson_mex') then

        if(t_sounding .ne. 42) then
           write(wrf_err_message,*) 'Inconsistent t_sounding type. ', &
           'Set t_sounding=42 if you really want to use a RS profile to init.'
           CALL wrf_error_fatal(TRIM( wrf_err_message ))
        endif

        passed_thru_but_found_no_file = .false.
        
        path = "./Data/soundings/"

        data_file   = trim(path)//rs_file

        write(wrf_err_message,*) '...processing a Mars Express RS profile'
        CALL wrf_message ( TRIM( wrf_err_message ) )
        write(wrf_err_message,*) data_file
        CALL wrf_message ( TRIM( wrf_err_message ) )

        ! info_type:   1 = read just the header and write into also to wrf_err
        !              2 = read just the profile and write info also to wrf_err
        !              3 = read both but no chatter to wrf_err
 
          open(unit=10,file=data_file,form='formatted',status='old')
          rewind(10)

          end_of_file = .false.
          k = 0

          if(info_type .eq. 1) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) ' - NOTE: These files should not be used ',&
                  'to set location information as many variables are missing!'
             CALL wrf_message ( TRIM( wrf_err_message ) )
          endif

          read(10,*) a,b, date_string, i1, i2, r1, r2, r3, r4, r5, r6, r7, r8, r9, &
                  r10, r11, r12, r13, r14, r15, r16

          if(present(sound_lat))  sound_lat = r3
          if(present(sound_lon))  sound_lon = r5
          if(present(sound_ls))   sound_ls  = r9
          if(present(sound_psfc)) sound_psfc = r12
          if(present(sound_time)) sound_time = r16

          if(info_type .eq. 1 .or. info_type .eq. 3) then
            write(wrf_err_message,*) ' '
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' Profile Lat = ', r3
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' Profile Lon = ', r5
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' Profile Ls = ', r9
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' Profile Psfc = ', r12
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' Profile Time = ', r16
            CALL wrf_message ( TRIM( wrf_err_message ) )
            write(wrf_err_message,*) ' '
            CALL wrf_message ( TRIM( wrf_err_message ) )
          endif

          if (info_type .eq. 1) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )

             read(date_string(1:4),*) i_year
             read(date_string(6:7),*) i_month
             read(date_string(9:10),*) i_day
             read(date_string(12:13),*) i_hour
             read(date_string(15:16),*) i_min
             read(date_string(18:19),*) i_sec

             write(wrf_err_message,*) date_string
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) i_year,i_month,i_day,i_hour,i_min,i_sec
             CALL wrf_message ( TRIM( wrf_err_message ) )

             dat(1) = i_year
             dat(2) = i_month
             dat(3) = i_day
             dat(4) = 0.     ! this is the offset from UTC
             dat(5) = i_hour
             dat(6) = i_min
             dat(7) = i_sec
             dat(8) = 0.     ! this is the milliseconds - don't care

             call d2j(dat,calc_jdate,calc_err)
             call mars24_j2000_offset_tt(calc_jdate,j_o)
             call mars24_Mars_Solar_Date(j_o,msd)
             call mars24_Mars_Year(j_o,mars_year,day_of_year)
             if(present(sound_doy)) sound_doy=real(day_of_year)
             if(present(sound_msd)) sound_msd=real(msd)

             write(wrf_err_message,*) 'returned jdate=',real(calc_jdate),real(msd)
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) 'mars year and doy=',real(mars_year),real(day_of_year)
             CALL wrf_message ( TRIM( wrf_err_message ) )
          endif
          if((info_type .eq. 2) .or. (info_type .eq. 3)) then ! grab the actual thermal profile

            k=0
            do while (.not. end_of_file)

              k = k+1
              if(k.gt.100) then
                    write(wrf_err_message,*) &
                    'woops, RS file has more points than expected!'
                    CALL wrf_error_fatal(TRIM( wrf_err_message ))
              endif
              read(10,*,end=500) r1,r2,r3,r4,r5,r6,r7
              if(k == 1) elev_lowest = (r4/g) - 1000.    ! lowest elevation = geop/grav and 
                                                        ! then surface is usually 1km lower
              dist(k) = r4/g - elev_lowest
              temp(k) = r7
              go to 510
 500          end_of_file = .true.
 510          continue
            enddo

            nl = k

            if(present(z_profile_out)) z_profile_out(nl) = 0.
            if(present(t_profile_out)) t_profile_out(nl) = temp(1)
            if(present(n_profile))     n_profile     = nl

            do k = 1, nl-1
               if(present(z_profile_out)) z_profile_out(k) = dist(nl-k)
               if(present(t_profile_out)) t_profile_out(k) = temp(nl-k)
            enddo

            if(info_type .eq. 2) then
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
             write(wrf_err_message,*) 'Radio Occ profile: (lev, T, z(m))'
             CALL wrf_message ( TRIM( wrf_err_message ) )
             do k = nl-1, 1, -1
               write(wrf_err_message,*) (nl-k),temp(k), dist(k)
               CALL wrf_message ( TRIM( wrf_err_message ) )
             enddo
             write(wrf_err_message,*) ' '
             CALL wrf_message ( TRIM( wrf_err_message ) )
            endif

          endif

          close(unit=10,status = 'keep')

        endif
! ---

        if(passed_thru_but_found_no_file) then

           write(wrf_err_message,*) 'ERROR: read_soundings passed though ',&
                'all file options and nothing matched your request.'
           CALL wrf_error_fatal(TRIM( wrf_err_message ))

        endif

      end subroutine read_soundings

! ------------------------------------------------------- CARMA Profiles ------------

      subroutine read_carma_profiles(load_atmos,p,t,n &
#ifdef USE_CARMA
                                ,ingest_carma_outfile  &
#endif 
                      )

! This is a re-written version of Erika Barth's column_titan.f90 but simplified 
! just to extract the input T vs. P profiles (rather than bundling this with
! interpolation onto a fixed height grid, which WRF does separately).

      implicit none

      integer, intent(in) :: load_atmos
#ifdef USE_CARMA
      character(len=*), intent(in) :: ingest_carma_outfile
#endif 

      real, intent(out), dimension(1000) :: p,t
      integer, intent(out) :: n

      real, dimension(1000) :: p_buffer, t_buffer

#ifdef USE_CARMA
    if(ingest_carma_outfile == "NONE") then
#endif
    select case( load_atmos )
     case( 1 )
       write(wrf_err_message,*) 'Set up Lellouch'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_lellouch(T,p,n)
     case( 2 ) 
       write(wrf_err_message,*) 'Set up HASI profiles'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_hasi(T,p,n)
     case( 3 )
       write(wrf_err_message,*) 'Use Schinder et al. 2011 RSS profile - T12 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder11('12i',T,p,n)
     case( 4 )
       write(wrf_err_message,*) 'Use Schinder et al. 2011 RSS profile - T12 egress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder11('12e',T,p,n)
     case( 5 )
       write(wrf_err_message,*) 'Use Schinder et al. 2011 RSS profile - T14 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder11('14i',T,p,n)
     case( 6 )
       write(wrf_err_message,*) 'Use Schinder et al. 2011 RSS profile - T14 egress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder11('14e',T,p,n)
     case( 7 )
       write(wrf_err_message,*) 'Use Anderson CIRS profile for 85 N'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_cirs(85,T,p,n)
     case( 8 )
       write(wrf_err_message,*) 'Use Anderson CIRS profile for 10 N'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_cirs(10,T,p,n)
     case( 9 )
       write(wrf_err_message,*) 'Use Anderson CIRS profile for 15 S'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_cirs(15,T,p,n)
     case( 10 )
       write(wrf_err_message,*) 'Use Anderson CIRS profile for 55 S'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_cirs(55,T,p,n)
     case( 11 )
       write(wrf_err_message,*) 'Use Vinatier CIRS profile for 87 S (old profile)'
       !-ELB-Keeping this for comparison; call for newer profile below
       !(use -1 for latitude value here to distinguish from new 87 S profile)
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_vinatier(-1,T,p,n)
     case( 12 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T27 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('27i',T,p,n)
     case( 13 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T27 egress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('27e',T,p,n)
     case( 14 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T31 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('31i',T,p,n)
     case( 15 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T31 egress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('31e',T,p,n)
     case( 16 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T46 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('46i',T,p,n)
     case( 17 )
       write(wrf_err_message,*) 'Use Schinder et al. 2012 RSS profile - T57 ingress'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_schinder12('57i',T,p,n)
     case( 18 )
       write(wrf_err_message,*) 'Use Vinatier CIRS profile for 87 S'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_vinatier(87,T,p,n)
     case( 19 )
       write(wrf_err_message,*) 'Use Vinatier CIRS profile for 78 S'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_vinatier(78,T,p,n)
     case( 20 )
       write(wrf_err_message,*) 'Use Vinatier CIRS profile for 73 S'
       CALL wrf_message ( TRIM( wrf_err_message ) )
       call get_atmos_profiles_vinatier(73,T,p,n)
     case default
      write(wrf_err_message,*) 'Defaulting to Lindal Voyager ingress profile'
      CALL wrf_message ( TRIM( wrf_err_message ) )
      write(wrf_err_message,*) 'Start at surface using 4 km layers' 
      CALL wrf_message ( TRIM( wrf_err_message ) )
      call get_atmos_profiles_lindal(T,p,n)
    end select

#ifdef USE_CARMA
    else ! in this case ingest_carma_outfile is not NONE:
       write(wrf_err_message,*) 'Getting t-p profiles from CARMA outfile: ',trim(ingest_carma_outfile)
       call wrf_message ( TRIM( wrf_err_message ) )
       call read_carma_tp_profile(T,p,n,ingest_carma_outfile)
    endif
#endif

    ! want to braket the profile such that the "data" always exceeds all possible
    ! model domain extents - isothermally extend up and down:

    ! sort so that we have p decreasing with increasing index
    p=-p
    CALL SORT2(p(1:n),t(1:n))
    p=-p

    p_buffer(2:n+1)=p(1:n)
    t_buffer(2:n+1)=t(1:n)

    p_buffer(1) = p_buffer(2)*100.
    t_buffer(1) = t_buffer(2)

    p_buffer(n+2) = 1.e-8
    t_buffer(n+2) = t_buffer(n+1)

    p = p_buffer
    t = t_buffer
    n=n+2

      end subroutine read_carma_profiles

      !```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_lindal(T,p,n_profile) 

!+ Titan pressure, temperature, air density from Lindal et al. 1983
!  Voyager 1 radio occultation measurements

implicit none

!== Arguments
real, dimension(:), intent(OUT) :: t,p
integer, intent(OUT) :: n_profile

!== Local declarations
integer :: j, k, nlayers
logical :: use_ingress_data

integer, parameter :: nz_data = 106
real, dimension(nz_data) :: t_ingress,t_egress,p_ingress, &
   p_egress,n_ingress,n_egress

!Data runs from 200 km to surface, mainly 2 km increments down to
!6 km, then 1 or 0.5 km increments (see altitude array set below)

![K]
data t_ingress / 169.4, 169.1, 168.8, 168.6, 168.3, 168.1, 167.9, &
   167.6, 167.4, 167.3, 167.2, 167.2, 167.2, 167.1, 167.0, 166.9, &
   166.9, 166.4, 165.7, 165.4, 165.3, 165.0, 164.5, 163.9, 163.1, &
   162.6, 162.3, 162.0, 161.5, 160.9, 160.1, 159.4, 159.0, 158.4, &
   157.8, 157.0, 156.1, 155.2, 154.0, 153.0, 152.1, 151.3, 150.7, & 
   150.0, 149.2, 148.5, 147.8, 147.0, 145.6, 144.3, 143.1, 141.9, &
   140.6, 139.0, 137.1, 135.9, 135.1, 133.2, 130.4, 127.4, 124.4, & 
   121.4, 118.4, 115.4, 110.8, 104.8,  99.1,  92.1,  85.4,  80.5, &
    77.5,  75.1,  73.6,  72.7,  72.2,  71.9,  71.6,  71.5,  71.4, &
    71.2,  71.2,  71.4,  71.7,  71.9,  72.3,  72.9,  73.5,  74.2, &
    75.1,  76.2,  77.0,  78.2,  79.5,  80.8,  82.2,  83.6,  85.3, &
    87.1,  88.1,  88.9,  89.9,  91.2,  91.9,  92.6,  93.3,  94.0 /

data t_egress  / 174.2, 172.0, 168.4, 164.2, 163.9, 165.4, 166.3, &
   167.0, 167.4, 167.6, 166.9, 165.7, 165.7, 166.7, 167.5, 168.2, &
   169.1, 170.4, 170.5, 169.2, 168.4, 168.1, 168.1, 168.2, 167.8, &
   166.8, 166.1, 165.6, 164.8, 163.5, 162.4, 161.7, 161.3, 161.4, &
   161.2, 160.8, 160.1, 158.9, 157.6, 156.2, 155.0, 154.1, 153.4, &
   152.8, 152.0, 151.1, 150.2, 149.5, 148.3, 147.0, 145.3, 143.6, &
   142.3, 140.9, 139.2, 137.6, 136.2, 134.9, 132.4, 129.7, 126.8, & 
   123.6, 120.4, 117.5, 112.6, 107.4, 101.8,  93.1,  86.3,  81.3, &
    77.8,  75.8,  74.3,  73.2,  72.5,  72.0,  71.8,  71.8,  71.4, &
    71.6,  71.7,  71.7,  71.9,  72.1,  72.4,  73.0,  73.6,  74.4, &
    75.2,  76.0,  77.2,  78.4,  79.5,  80.8,  82.4,  83.7,  85.4, &
    87.0,  88.0,  89.0,  90.2,  91.4,  92.0,  92.7,  93.3,  93.9 /

![mbar]
data p_ingress  /    0.75, 0.79, 0.83, 0.87, 0.91, 0.95, 1.00, &
   1.04, 1.10, 1.15, 1.20, 1.26, 1.33, 1.39, 1.46, 1.53, 1.61, &
   1.69, 1.77, 1.86, 1.95, 2.05, 2.15, 2.26, 2.38, 2.50, 2.63, &
   2.76, 2.91, 3.06, 3.22, 3.39, 3.57, 3.76, 3.96, 4.17, 4.40, &
   4.64, 4.89, 5.17, 5.46, 5.76, 6.09, 6.44, 6.81, 7.21, 7.63, &
   8.08, 8.56, 9.07, 9.62, 10.21, 10.84, 11.52, 12.25, 13.03, &
   13.88, 14.80, 15.79, 16.88, 18.06, 19.38, 20.81, 22.41, 24.20, &
   26.21, 28.54, 31.25, 34.48, 38.31, 42.78, 47.97, 53.97, 60.84, &
   68.67, 77.60, 87.75, 99.29, 112.40, 127.30, 144.26, 163.48, &
   185.22, 209.84, 237.68, 269.08, 304.45, 344.16, 388.68, 438.44, &
   494.04, 556.03, 624.87, 701.09, 785.43, 878.54, 980.98, &
   1093.37, 1153.44, 1216.28, 1282.01, 1350.40, 1385.60, 1421.48, &
   1458.02, 1495.26 /

data p_egress   /    0.80, 0.84, 0.88, 0.92, 0.97, 1.02, 1.07, &
   1.12, 1.17, 1.23, 1.29, 1.35, 1.42, 1.49, 1.56, 1.64, 1.72, &
   1.80, 1.89, 1.98, 2.08, 2.18, 2.29, 2.41, 2.52, 2.65, 2.78, &
   2.92, 3.07, 3.23, 3.40, 3.57, 3.76, 3.96, 4.16, 4.38, 4.62, &
   4.86, 5.12, 5.40, 5.70, 6.01, 6.35, 6.71, 7.09, 7.49, 7.92, &
   8.37, 8.86, 9.38, 9.94, 10.54, 11.19, 11.88, 12.62, 13.42, &
   14.29, 15.21, 16.22, 17.32, 18.52, 19.83, 21.28, 22.88, 24.66, &
   26.69, 28.98, 31.68, 34.92, 38.73, 43.24, 48.44, 54.44, 61.31, &
   69.18, 78.15, 88.36, 99.94, 113.10, 128.06, 144.97, 164.17, &
   185.96, 210.62, 238.51, 269.95, 305.35, 345.15, 389.75, &
   439.69, 495.42, 557.38, 626.23, 702.64, 787.11, 880.38, &
   983.00, 1095.51, 1155.78, 1218.73, 1284.34, 1352.68, 1387.87, &
   1423.77, 1460.34, 1497.59 /

![ x 10^-17 cm^-3]
data n_ingress / 0.32, 0.34, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, &
   0.47, 0.50, 0.52, 0.55, 0.57, 0.60, 0.63, 0.66, 0.70, 0.73, &
   0.77, 0.81, 0.85, 0.90, 0.95, 1.00, 1.06, 1.11, 1.17, 1.24, 1.30, &
   1.38, 1.46, 1.54, 1.63, 1.72, 1.82, 1.92, 2.04, 2.17, 2.30, &
   2.45, 2.60, 2.76, 2.93, 3.11, 3.31, 3.52, 3.74, 3.98, 4.26, & 
   4.55, 4.87, 5.21, 5.58, 6.00, 6.47, 6.95, 7.44, 8.05, 8.77, &
   9.60, 10.51, 11.57, 12.74, 14.07, 15.83, 18.14, 20.87, 24.61, &
   29.29, 34.51, 40.07, 46.34, 53.23, 60.73, 69.07, 78.48, &
   89.07, 101.09, 114.64, 130.33, 147.77, 167.00, 188.70, 213.27, &
   240.49, 270.56, 304.10, 340.67, 380.78, 423.99, 473.54, &
   525.79, 582.13, 643.58, 710.23, 781.93, 856.89, 937.22, &
   979.00, 1023.97, 1067.09, 1108.70, 1129.37, 1150.16, 1171.26, &
   1192.38 /

data n_egress / 0.33, 0.35, 0.38, 0.41, 0.43, 0.45, 0.46, 0.49, &
   0.51, 0.53, 0.56, 0.59, 0.62, 0.65, 0.68, 0.71, 0.74, 0.77, &
   0.80, 0.85, 0.90, 0.94, 0.99, 1.04, 1.09, 1.15, 1.21, 1.28, &
   1.35, 1.43, 1.51, 1.60, 1.69, 1.78, 1.87, 1.97, 2.09, 2.22, &
   2.35, 2.51, 2.66, 2.83, 3.00, 3.18, 3.38, 3.59, 3.82, 4.06, &
   4.33, 4.63, 4.96, 5.32, 5.70, 6.10, 6.57, 7.07, 7.60, 8.17, &
   8.88, 9.67, 10.58, 11.63, 12.81, 14.11, 15.87, 18.02, 20.64, &
   24.67, 29.35, 34.54, 40.31, 46.40, 53.19, 60.85, 69.32, 78.92, &
   89.53, 101.29, 115.39, 130.31, 147.36, 167.19, 189.02, 213.71, &
   241.00, 271.02, 304.26, 340.84, 381.27, 426.53, 473.44, 525.40, &
   583.47, 644.66, 710.04, 782.98, 857.63, 940.30, 982.00, &
   1024.11, 1065.72, 1108.49, 1130.12, 1151.27, 1172.23, 1192.92 /

  n_profile=nz_data

!== Executable statements:

 !Set switch for ingress or egress data
  use_ingress_data = .true.

  write(wrf_err_message,*) 'Lindal profile:'
  CALL wrf_message ( TRIM( wrf_err_message ) )
  if( use_ingress_data ) then
      T(1:n_profile) = t_ingress(:)
      p(1:n_profile) = p_ingress(:) * 100. !mb_to_Pa
  else !use_egress_data
      T(1:n_profile) = t_egress(:)
      p(1:n_profile) = p_egress(:) * 100.  !mb_to_Pa
  endif
     
end subroutine get_atmos_profiles_lindal

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_lellouch(T,p,n_profile) 

!+ Titan atmosphere pressure, temperature, air density from the
!  Lellouch and Hunten Engineering model (1987)

!  Use these profiles to match up with runs done for thesis

!~ 14.07.2011 (last modified; TIM)

implicit none

!== Arguments
real, dimension(:), intent(OUT) :: t, p ![nz_ext]
integer, intent(OUT) :: n_profile

!== Local declarations
integer :: k, i, nz_ext, luni, nz_data

integer, parameter :: nz_data10 = 60, &
                      nz_data2  = 50
real :: dz_data, airdens 
real, dimension(nz_data10) :: t_data, rhoa_data, zc_data

logical :: deep_profile

deep_profile = .true.

!== Executable statements:

! Data grid runs from surface to TOA

  if( deep_profile ) then

   !Data grid is 10 km layers; layer bottom beginning at the surface
   !Top layer ends at 600 km
    nz_data = nz_data10
    dz_data = 10.e3
    zc_data(1) = 5.e3 ! 5 km
    do k=2,nz_data10
      zc_data(k) = zc_data(k-1) + dz_data
    enddo
    luni = 10
    open(unit=luni,file=carma_data_path//'LH87_10km_midpt.txt', &
           status='old')
    do k=1,nz_data10
      read(luni,*) T_data(k), rhoa_data(k)
    enddo
    close(luni)

  else

  !Data grid is 2 km layers; layer bottom beginning at the surface
  !Top layer ends at 100 km
    nz_data = nz_data2
    dz_data = 2.e3
    zc_data(1) = 1.e3 ! 1 km
    do k=2,nz_data2
      zc_data(k) = zc_data(k-1) + dz_data
    enddo
    luni = 10
    open(unit=luni,file=carma_data_path//'LH87_2km_midpt.txt', &
           status='old')
    do k=1,nz_data2
      read(luni,*) T_data(k), rhoa_data(k)
    enddo
    close(luni)

  endif

  do k=1,nz_data
    t(k) = t_data(k)
    airdens = rhoa_data(k) /1000.  ! -> g_cm3_to_kg_m3
    p(k) = airdens * R_air * T(k)
  enddo

  n_profile = nz_data

end subroutine get_atmos_profiles_lellouch

subroutine get_atmos_profiles_schinder11(pid,T,p,n_profile) 

!+ Titan atmosphere pressure, temperature, air density from 
!  Schinder et al. 2011
!
!  T12 (Mar 19, 2006) - ingress at 31.4*S, egress at 52.8*S
!  T14 (May 20, 2006) - ingress at 32.7*S, egress at 34.3*S

!~ ELB 10/05/2012

implicit none

!== Arguments
character(len=3), intent(IN) :: pid !profile id

real, dimension(:), intent(OUT) :: t, p
integer, intent(out) :: n_profile

!== Local declarations
integer :: k, l, nl, luni

integer, parameter :: nz_data = 90

real, dimension(nz_data) :: t_data, p_data, zc_data, &
  Ti_data, Te_data, pi_data, pe_data

real :: refi, refe

character(len=2) :: flyby
character(len=1) :: sounding

!== Executable statements:

! Data grid runs from TOA to surface

 ! Deconstruct input profile id string to get flyby number 
 ! and ingress/egress sounding
  flyby = pid(1:2)
  sounding = pid(3:3)
  
  nl = nz_data/2 !Number of table lines in file

  luni = 10
  if( flyby == '12' ) then
    open(unit=luni,file=carma_data_path//'schinder_2011_T12.txt', &
           status='old')
  elseif( flyby == '14' ) then
    open(unit=luni,file=carma_data_path//'schinder_2011_T14.txt', &
           status='old')
  else
    write(wrf_err_message,*) 'Schinder 2011 RSS: Incorrect profile id ',pid 
    CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
  endif

  do l=1,5
    read(luni,*) !Skip file header
  enddo
  do k=1,nl
    read(luni,*) zc_data(k), &
       Ti_data(k), pi_data(k), refi, Te_data(k), pe_data(k), refe, &
       zc_data(k+nl), Ti_data(k+nl), pi_data(k+nl), refi, &
       Te_data(k+nl), pe_data(k+nl), refe 
  enddo
  close(luni)

 !*** Set to ingress or egress data ***
 if( sounding == 'i') then
   p_data(:) = pi_data(:)  
   T_data(:) = Ti_data(:)
 elseif( sounding == 'e') then
   p_data(:) = pe_data(:)  
   T_data(:) = Te_data(:)
 else
   write(wrf_err_message,*) 'Schinder 2011 RSS: Incorrect profile id ',pid 
   CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
 endif

 !Convert pressures mbar->Pa
  p_data(:) = 100. * p_data(:)

 !Convert altitudes km->m
  zc_data(:) = 1000. * zc_data(:)

  do k=1,nl
    p(k) = p_data(k)
    t(k) = t_data(k)
  enddo
  n_profile = nl

end subroutine get_atmos_profiles_schinder11

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_schinder12(pid,T,p,n_profile) 

!+ Titan atmosphere pressure, temperature, air density from 
!  Schinder et al. 2012
!
!  T27 (Mar 26, 2007) - ingress at 69.0*S, egress at 52.9*N
!  T31 (May 28, 2007) - ingress at 74.3*S, egress at 74.1*S
!  T46 (Nov 3, 2008)  - ingress at 32.4*S
!  T57 (June 22, 2009)- ingress at 79.8*N

!~ ELB 07/05/2016

implicit none

!== Arguments
character(len=3), intent(IN) :: pid !profile id

real, dimension(:), intent(OUT) :: t, p
integer, intent(out) :: n_profile

!== Local declarations
integer :: k, l, nl, na, nzd, luni

integer, parameter :: nz_data = 90

real, dimension(nz_data) :: t_data, p_data, zc_data, &
  Ti_data, Te_data, pi_data, pe_data

real :: refi, refe

character(len=2) :: flyby
character(len=1) :: sounding

!== Executable statements:

! Data grid runs from TOA to surface

 ! Deconstruct input profile id string to get flyby number 
 ! and ingress/egress sounding
  flyby = pid(1:2)
  sounding = pid(3:3)
  

  luni = 10
  if( flyby == '27' ) then
    open(unit=luni,file=carma_data_path//'schinder_2012_T27.txt', &
           status='old')
    nl = 44 !Number of full table lines in file
    na = 2  !Number of additional (partial) lines
  elseif( flyby == '31' ) then
    open(unit=luni,file=carma_data_path//'schinder_2012_T31.txt', &
           status='old')
    nl = 44 !Number of full table lines in file
    na = 2  !Number of additional (partial) lines
  elseif( flyby == '46' .or. flyby == '57' ) then
    open(unit=luni,file=carma_data_path//'schinder_2012_T46_T57.txt', &
           status='old')
   !Both soundings are ingress, but set 'i' and 'e' to use correct columns
    if( flyby == '46' ) sounding = 'i'
    if( flyby == '57' ) sounding = 'e'
    nl = 42 !Number of full table lines in file
    na = 3  !Number of additional (partial) lines
  else
    write(wrf_err_message,*) 'Schinder 2012 RSS: Incorrect profile id ',pid 
    CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
  endif

  do l=1,3
    read(luni,*) !Skip file header
  enddo
  do k=1,nl
    read(luni,*) zc_data(k), &
       Ti_data(k), pi_data(k), refi, Te_data(k), pe_data(k), refe, &
       zc_data(k+nl+na), Ti_data(k+nl+na), pi_data(k+nl+na), refi, &
       Te_data(k+nl+na), pe_data(k+nl+na), refe 
  enddo
  do k=nl+1,nl+na
    read(luni,*) zc_data(k), &
       Ti_data(k), pi_data(k), refi, Te_data(k), pe_data(k), refe
  enddo
  close(luni)

 nzd = 2*nl+na
 p_data(:) = 1500. !Set some default values for cases with less 
 T_data(:) = 0.    ! 90 data points
 !*** Set to ingress or egress data ***
 if( sounding == 'i') then
   p_data(1:nzd) = pi_data(1:nzd)  
   T_data(1:nzd) = Ti_data(1:nzd)
 elseif( sounding == 'e') then
   p_data(1:nzd) = pe_data(1:nzd)  
   T_data(1:nzd) = Te_data(1:nzd)
 else
   write(wrf_err_message,*) 'Schinder 2012 RSS: Incorrect profile id ',pid 
   CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
 endif

   !print *,'pressure:',p_data(:)
   !print *,'temperature:',T_data(:)

 !Convert pressures mbar->Pa
  p_data(:) = 100. * p_data(:)

 !Convert altitudes km->m
  zc_data(:) = 1000. * zc_data(:)

  do k=1,nzd
    p(k) = p_data(k)
    t(k) = t_data(k)
  enddo
  n_profile = nzd

end subroutine get_atmos_profiles_schinder12

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_niemann(T,p,n_profile)

!+ Temperature/pressure profiles from Niemann et al. 2010
!  (Composition of Titan’s lower atmosphere and simple surface
!   volatiles as measured by the Cassini‐Huygens probe gas
!   chromatograph mass spectrometer experiment; JGR v. 115, E12006)
!
!~ ELB 08/20/2012

implicit none

!== Arguments:
real, dimension(:), intent(OUT) :: T, p
integer, intent(out) :: n_profile

!== Local declarations:
integer, parameter :: nz_tab1 = 35, &
                      nz_tab2 = 30, &
                      nz_data = nz_tab1 + nz_tab2 - 3
integer :: k, l, luni, ich4, ih2, iz, nz_mod
real :: stdev  !dummy variable for table entries
integer, dimension(nz_tab1) :: time1 !use time entries for
integer, dimension(nz_tab2) :: time2 !combining altitudes
real, dimension(nz_tab1) :: xCH4, z_ch4, p_ch4, T_ch4
real, dimension(nz_tab2) :: xH2, z_h2, p_h2, T_h2
real, dimension(nz_data) :: z_data, p_data, T_data

!== Executable statements:

  write(wrf_err_message,*) "Don't use this for p, T profiles, use HASI data &
  &instead. This has data gap 45-75 km and bad data point in H2" 
  CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
  !!! Rewrite this subroutine to only use for CH4, H2 profiles !!!

 !Read in p, T data from CH4 table (Table 1)
  luni = 10
  open(unit=luni,file=carma_data_path//'niemann_2010_table1.txt', &
         status='old')
  do l=1,9
    read(luni,*) !Skip file header
  enddo
  do k=1,nz_tab1
    read(luni,*) time1(k), xCH4(k), stdev, z_ch4(k), p_ch4(k), T_ch4(k)
  enddo
  close(luni)

 !Read in p, T data from H2 table (Table 2)
  luni = 10
  open(unit=luni,file=carma_data_path//'niemann_2010_table2.txt', &
         status='old')
  do l=1,9
    read(luni,*) !Skip file header
  enddo
  do k=1,nz_tab2
    read(luni,*) time2(k), xH2(k), stdev, z_h2(k), p_h2(k), T_h2(k)
  enddo
  close(luni)
  
 !Table values run from altitude to surface; Combine independent z 
 !values from both tables (note: both contain z=0.0, 2.0, 128.6 points 
 !so total number of altitudes is nz_tab1 + nz_tab2 - 3
  ich4=1 ; ih2=1
  do k=1,nz_data
    if( ih2 <= nz_tab2 .and. time2(ih2) < time1(ich4) ) then 
      z_data(k) = z_h2(ih2)
      T_data(k) = T_h2(ih2)
      p_data(k) = p_h2(ih2)
      if( abs(time2(ih2)-time1(ich4)) <= 1) ich4 = ich4+1 !skip duplicate altitude
      ih2 = ih2+1
    else
      if( ich4 <= nz_tab1 ) then
        z_data(k) = z_ch4(ich4)
        T_data(k) = T_ch4(ich4)
        p_data(k) = p_ch4(ich4)
      if( abs(time2(ih2)-time1(ich4)) <= 1) ih2 = ih2+1  !skip duplicate altitude
        ich4 = ich4+1
      endif
    endif
  !print *,k,z_data(k),T_data(k),p_data(k)
  enddo

 !Convert pressures hPa->Pa
  p_data(:) = 100. * p_data(:)

 !Convert altitudes km->m
  z_data(:) = 1000. * z_data(:)

  do k=1,nz_data
    p(k) = p_data(k)
    t(k) = t_data(k)
  enddo
  n_profile = nz_data

end subroutine get_atmos_profiles_niemann

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_hasi(T,p,n_profile)

!+ Temperature/density profiles from HASI data below 150 km

!~ ELB 02/14/2012

implicit none

!== Arguments:
real, dimension(:), intent(OUT) :: T, p
integer, intent(out) :: n_profile

!== Local declarations:
real, parameter :: dz = 2.e3
integer :: luni, l, k, iz
integer :: nz_ext
real :: Tskip, rhoskip

integer, parameter :: nhlo = 1810, &  !number of data lines in low altitude HASI file
                      nhhi = 978,  &  !number of data lines in high altitude HASI file 
                      nhasi = nhlo+nhhi 
real, dimension(nhasi) :: T_hasi, z_hasi, rho_hasi

  ! Note: Units listed in files for rho are incorrect; rho values are in mks

!== Executable statements:

   !Read in HASI T, rho data (0-150 km)
    luni = 10
    open(unit=luni,file=carma_data_path//'HASI-TEM_temperature_9C580.txt', &
           status='old')
    do l=1,7
      read(luni,*) !Skip file header
    enddo
    do k=1,nhlo
      read(luni,*) T_hasi(k), z_hasi(k)  ![K], [km]
    enddo
    close(luni)

    luni = 10
    open(unit=luni,file=carma_data_path//'HASI_density_profile_9C585.txt', &
           status='old')
    do l=1,7
      read(luni,*) !Skip file header
    enddo
    do k=1,nhlo
      read(luni,*) rho_hasi(k), z_hasi(k)  ![kg/m3], [km]
    enddo
    close(luni)

   !Read in HASI T, rho data (160 - 1500 km)
    luni = 10
    open(unit=luni,file=carma_data_path//'HASI-TEM_temperature_9C57F.txt', &
           status='old')
    do l=1,7
      read(luni,*) !Skip file header
    enddo
     ! The high altitude HASI file contains three derived temperature 
     ! profiles.  Values don't differ significantly below ~1200 km so
     ! just picking one for now since CARMA models usually run at
     ! lower altitudes.
    do k=nhasi,nhlo+1,-1
      read(luni,*) z_hasi(k), T_hasi(k), Tskip, Tskip  ![km], [K]
    enddo
    close(luni)

    luni = 10
    open(unit=luni,file=carma_data_path//'HASI_density_profile_9C584.txt', &
           status='old')
    do l=1,7
      read(luni,*) !Skip file header
    enddo
     ! See comment above for T-data
    do k=nhasi,nhlo+1,-1
      read(luni,*) z_hasi(k), rho_hasi(k), rhoskip, rhoskip  ![km], [kg/m3] 
    enddo
    close(luni)


  z_hasi(:) = z_hasi*1000. !convert to m

  do k=1,nhasi
    t(k) = t_hasi(k)
    p(k) = rho_hasi(k) * R_air * T(k)
  enddo
  n_profile = nhasi

end subroutine get_atmos_profiles_hasi

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_cirs(lat,T,p,n_profile)

!+ Temperature/density profiles from CIRS data (supplied by C. Anderson)

!~ ELB 12/18/2012

implicit none

!== Arguments:
integer, intent(IN) :: lat
real, dimension(:), intent(OUT) :: T, p
integer, intent(out) :: n_profile

!== Local declarations:
integer :: luni, l, k, iz
integer :: nz_ext
real :: qskip

integer, parameter :: ncirs = 101
real, dimension(ncirs) :: T_cirs, p_cirs, z_cirs, rho_cirs

!== Executable statements:

  luni = 10

  select case( lat )
   case( 85 ) 
    open(unit=luni,file=carma_data_path//'85NT4_ptprof.out',status='old')
    read(luni,*) !Skip file header
    do k=1,ncirs
      read(luni,*) z_cirs(k), T_cirs(k), p_cirs(k), rho_cirs(k), qskip
    enddo
   case( 10 )
    open(unit=luni,file=carma_data_path//'10NT26_ptprof.out',status='old')
    do k=1,ncirs
      read(luni,*) z_cirs(k), T_cirs(k), p_cirs(k), rho_cirs(k)
    enddo
   case( 15 ) 
    open(unit=luni,file=carma_data_path//'15ST17_ptprof.out',status='old')
    read(luni,*) !Skip file header
    do k=1,ncirs
      read(luni,*) z_cirs(k), T_cirs(k), p_cirs(k), rho_cirs(k), qskip
    enddo
   case( 55 )
    open(unit=luni,file=carma_data_path//'55ST6_ptprof.out',status='old')
    do k=1,ncirs
      read(luni,*) z_cirs(k), T_cirs(k), p_cirs(k), rho_cirs(k)
    enddo
  end select
  close(luni)

  z_cirs(:) = z_cirs*1000. !convert to m

  !Convert pressure [mbar]->[Pa]
   p_cirs(:) = p_cirs(:) * 100. ! -> mb_to_Pa

  !Convert density [#/cm3] -> [kg/m3]
  !airdens(:) = 1.e6 * wtmol_air/avogadnum * airdens(:) ! ->cm3_to_m3

  do k=1,ncirs
    p(k)=p_cirs(k)
    t(k)=t_cirs(k)
  enddo
  n_profile=ncirs

end subroutine get_atmos_profiles_cirs

!```````````````````````````````````````````````````````````````````
subroutine get_atmos_profiles_vinatier(lat,T,p,n_profile)

!+ Temperature/density profiles from CIRS data (supplied by S. Vinatier) 

!~ ELB 11/22/2015

implicit none

!== Arguments:
integer, intent(IN) :: lat
real, dimension(:), intent(OUT) :: T, p
integer, intent(out) :: n_profile

!== Local declarations:
integer :: luni, k

integer, parameter :: ncirs = 96
real, dimension(ncirs) :: T_cirs, p_cirs, z_cirs
real :: rskip

!== Executable statements:

  luni = 10

 !Use old CIRS pressure-temperature file to get altitude values
 !(Note: Newer files have one additional data point, at the surface)
  open(unit=luni,file=carma_data_path//'profT_87S_May_2015.txt',status='old')
  do k=1,4
    read(luni,*) !Skip file header
  enddo
  do k=1,ncirs-1
    !                z - km; p - mbar; T - K
    read(luni,*) z_cirs(k), p_cirs(k), T_cirs(k)
  enddo
  close(luni)

  z_cirs(:) = z_cirs*1000. !convert to m

 !Latitude files below are written "surface" to "TOA"
 !Reverse direction at read-in
  if( lat > 0 ) then
    select case (lat)
      case ( 87 ) 
        open(unit=luni,file=carma_data_path//'1305_T_nad_87S_v3p2_21.dat',status='old')
      case ( 78 ) 
        open(unit=luni,file=carma_data_path//'1305_T_nad_78S_v4p3_6.dat',status='old')
      case ( 73 ) 
        open(unit=luni,file=carma_data_path//'1305_T_nad_73S_v4p3_5.dat',status='old')
    end select

    z_cirs(ncirs) = 0.
    do k=ncirs,1,-1
      !                p - mbar; T - K; x  ;  x
      read(luni,*) p_cirs(k), T_cirs(k), rskip, rskip
    enddo
    close(luni)

  endif

 !Convert pressure units
  p_cirs(:) = p_cirs(:)*100. !-> mb_to_Pa

  do k=1,ncirs
    p(k)=p_cirs(k)
    t(k)=t_cirs(k)
  enddo
  n_profile=ncirs
 
end subroutine get_atmos_profiles_vinatier

#ifdef USE_CARMA
!```````````````````````````````````````````````````````````````````````````````
! check routine needed by netcdf code in read_carma_tp_profile

  subroutine check(status)

    use netcdf

    implicit none

    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      WRITE( wrf_err_message , * ) "Terminal NETCDF ERROR: ",trim(nf90_strerror(status))
      CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    end if
  end subroutine check  

!```````````````````````````````````````````````````````````````````````````````

subroutine read_carma_tp_profile(T,p,n_profile,ingest_carma_outfile)

use netcdf

implicit none

!== Arguments:
character(len=*), intent(in) :: ingest_carma_outfile
real, dimension(:), intent(OUT) :: T, p
integer, intent(out) :: n_profile

!== Local declarations:
integer :: k
integer :: ncid, varid
integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in
integer :: dimid_t, dimid_z, n_t, n_z
integer, parameter :: nalt_expected = 60
character(len=90) :: a  ! dummy array for storing stuff we don't need
real, dimension(nalt_expected) :: z_in, p_in, t_in

!== Executable statements:

! Open the file read-only:
  call check( nf90_open(trim(ingest_carma_outfile), NF90_NOWRITE, ncid) )
  call check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )

! Examine the dimensions:
  call check( nf90_inq_dimid(ncid, "ntime", dimid_t))
  call check( nf90_inq_dimid(ncid, "nz_mic", dimid_z))
  call check( nf90_inquire_dimension(ncid, dimid_t, a, n_t))
  call check( nf90_inquire_dimension(ncid, dimid_z, a, n_z))

  if( (n_z /= nalt_expected) ) then
    WRITE( wrf_err_message , * ) &
'CARMA input altitude range not expected. Compare defaults in read_carma_tp_profile in module_read_soundings.F'
    CALL wrf_error_fatal( TRIM( wrf_err_message ) )
  endif

! Read in height, pressure, and temperature:
  call check( nf90_inq_varid(ncid, 'zc', varid) )
  call check( nf90_get_var(ncid, varid, z_in, start=(/ 1/)) )
  call check( nf90_inq_varid(ncid, 'p', varid) )
  call check( nf90_get_var(ncid, varid, p_in, start=(/ 1/)) )
  call check( nf90_inq_varid(ncid, 't', varid) )
  call check( nf90_get_var(ncid, varid, t_in, start=(/ 1/)) )

  ! Close file
  call check( nf90_close(ncid) )

  n_profile=nalt_expected

  do k=1,n_profile
    t(k)=t_in(k)
    p(k)=p_in(k)
  enddo

end subroutine read_carma_tp_profile
#endif

end module module_read_soundings
