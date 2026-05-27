      subroutine thin_ci_chk_ir(pxldat,vza,cirrus_ir,qa_bits,testbits)


      implicit none
      save

!---------------------------------------------------------------------
!!F77
!
!!Description:
! ... Routine to test for thin cirrus using IR channels.  This 
! ... will indicate whether the we believe the cirrus is thin
! ... enough for most tests to be applied without affecting
! ... results.  It will allow PI's with algorithms which are
! ... very sensitive to thin cirrus contamination to see if
! ... it might be there without affecting the final cloud mask
! ... confidence.
!
!!Input parameters:
! pxldat        Array containing reflectance or brightness temperatures
!               for all bands for a single pixel
! vza           Current pixel viewing angle
!
!!Output Parameters:
! cirrus_ir     Logical variable flagging thin cirrus contaminated
!               scenes in the infrared
! qa_bits       10 byte array contining qa bit results
! testbits      6 byte array containing bit results
!
!!Revision History:
! 10/04  Collection 5b   R. Frey
! Corrected thresholds.
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!END
!---------------------------------------------------------------------

      include 'global.inc'

! ... scalar arguments ..
      real vza
      logical cirrus_ir

! ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

! ... local scalars
      real cosvza,pi,dtr,masdf1,masir11,masir12,schi,dfthrsh,ci1,ci2,   &
           diftemp
      logical code
      integer debug,h_output
      real, parameter :: Rel_equality_EPS = 0.000001

! ... external subroutines ..
      external tview,clear_bit,set_qa_bit

! ... intrinsic functions
      intrinsic cos

! ... Common statement for debug purposes
!      common / bug / debug, h_output

!     Routine which checks for the presence of thin cirrus. This check 
!     is made independently of other spectral tests which may check
!     for similar conditions.
!
!     Check to see if IR thin cirrus bit should be set.
!     Right now this is based upon the APOLLO thin cirrus 11-12 BTDIF.
!     This test has been fairly robust over all but snow covered 
!     regions.

! ... assignment statements
      masir11 = pxldat(24)
      masir12 = pxldat(25)

! ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0
      masdf1 = 0.0
      cosvza = 0.0
      diftemp = 0.0
      schi = 0.0 
      dfthrsh = 0.0
      code = .true.
      ci1 = 0.0
      ci2 = 0.0

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Processing thin_ci_chk_ir routine'',/)')
!      endif
! ................................................................


! ... 11-12um brightness temperature difference test
! ... for low clouds).
      if (nint(masir11) .ne. nint(bad_data) .and.   &
          nint(masir12) .ne. nint(bad_data) .and.   &
          vza .gt. 0.0) then

        masdf1 = masir11 - masir12

! ...   11-12um brightness temperature difference test
! ...   for thin cirrus).
! ...   added apollo viewing angle/masir11 regressed threshold.
! ...   calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza).gt.Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if

! ...   interpolate look-up table values of 11 - 12 micron bt
! ...   difference thresholds (function of viewing zenith
! ...   and 11 micron brightness temperature).
        call tview(1,schi,masir11,diftemp)

! ...   if a threshold was determined by apollo, then use this
! ...   as the thin cirrus test, otherwise use a standard threshold
! ...   else don't use this threshold
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          code = .false.
        else
          dfthrsh = diftemp
        endif
      
! ...   Want to use a threshold range of very thin cirrus.
        if (code) then
          call set_qa_bit(qa_bits,11)
          ci1 = dfthrsh
          ci2 = dfthrsh + (0.3 * dfthrsh) 
          if (masdf1 .gt. ci1 .and. masdf1 .le. ci2) then
             call clear_bit(testbits,11)
             cirrus_ir = .true.
          endif
        endif

      endif

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x,'' thin_cirrus_ir vars:'',/,6f10.2,/)')
!     *        masir11,masir12,masdf1,vza,schi,diftemp
!        write(h_output,'(10x,'' more variables:'',/,l4,3f10.2,l4,/)')
!     *        code,dfthrsh,ci1,ci2,cirrus_ir
!      endif
! .................................................................

      return
      end

