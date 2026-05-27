      subroutine proc_path(water,land,day,ice,snow,snglnt,coast,
     +                     desert,smoke,shadow,testbits)

      implicit none
      save

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c ... routine for determining processing path through
c ... the algorithm and setting the appropriate bit
c ... flags. This is where the appropriate day/night
c ... and land/water bits are set. Also sets bit for
c ... smoke and shadow.
c
c!Input parameters:
c water         Logical variable - true if water background
c land          Logical variable - true if land background
c day           Logical variable - true if sza < 80
c ice           Logical variable - true if ice background
c snow          Logical variable - true if snow background
c snglnt        Logical variable - true if pixel is sun glint 
c               contaminated
c coast         Logical variable - true if coast background
c desert        Logical variable - true if desert background
c smoke         Logical variable - true if smoke was detected
c shadow        Logical variable - true if shadow was detected
c
c!Output Parameters:
c testbits      6 cell byte array containing cloud mask bit results
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

c ... scalar arguments ..
      logical water,land,day,ice,snow,snglnt,coast,desert,smoke,
     +         shadow

c ... array arguments
      byte testbits(6)

c ... local scalars ..
      integer debug,h_output

c     external subroutines
      external set_bit,clear_bit

c ... Common statement for debug purposes
!      common / bug / debug, h_output

c ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within proc_path testing routine '',/)')
!      endif
c ...............................................................

c ...
c     Set snow/ice bit.
      if((.not. snow) .and. (.not. ice))  then
         call set_bit(testbits,5)
      end if

c ... Set day/night flag
      if (day) then
         call set_bit(testbits,3)
      endif

c ... Set sunglint flag
      if( (.not. snglnt) .or. (.not. water) ) then
c        If pixel is not in geometric sun glint region, set bit
         call set_bit(testbits,4)
      end if

c     Set coast, desert, or land processing path flags.  Default is
c     water (00), which is set during initialization of the bit fla
c ... Now set land/sea bits
      if (coast) then
c       Set coast bit.
        call set_bit(testbits,6)
      else if (desert) then
c       Set desert bit.
        call set_bit(testbits,7)
      else if (land) then
c        Set "land" bits.
         call set_bit(testbits,6)
         call set_bit(testbits,7)
      end if

c     Now set the NCO and the shadow bits, if appropriate
      if (shadow) then
         call clear_bit(testbits,10)
      endif

c     Set the Non-cloud obstruction flag (currently related to smoke)
      if (smoke) then
         call clear_bit(testbits,8)
      endif

c ......  debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(15x,'' PROC_PATH variables '',/)')
!        write(h_output,'(2x,'' Day  Ice  Desert  Snow  Snglnt Land Water
!     + Coast Smoke'',/,3L5,2x,8L6)') day,ice,desert,snow,snglnt,land,
!     + water,coast,smoke
!      endif
c ................................................................


      return
      end
