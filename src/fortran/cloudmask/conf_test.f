      subroutine conf_test(val,locut,hicut,power,midpt,nmval,
     *                      conflev)

      implicit none
      save

c---------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine for determining the level of confidence of a
c     particular clear sky test.  Input single threshold value
c     ('midpt') with associated confidence limits ('locut', 'hicut').
c     Low" and "high" cutoffs
c     refer to low or high confidence ends of an interval and 
c     not necessarily to absolute value.  Routine calculates the
c     confidence based on an "S" function.  One may change the shape
c     of the function by changing 'power' and/or 'midpt'.
c     This version currently accepts only 1 set of limits to
c     process
c
c!Input Parameters:
c val      current individual test value
c locut    low confidence cutoff value (less then this is 0 conf.)
c hicut    high confidence cutoff value (greather than this is
c               100% conf.)
c power    S function curve power
c midpt    midpoint of S curve (currently 50% - test straight threshold)
c nmval    1 or 2 threshold values for this test?
c
c!Output Parameters:
c conflev  calculated confidence that fov is unobstructed for this test
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!End
c----------------------------------------------------------------------

c     scalar arguments
      integer nmval
      real power,val,conflev,locut,hicut,midpt

c     local scalars
      real alpha,gamma,range,coeff,s1,c,beta
      logical flipped

      coeff = 2.0 ** (power - 1.0)

c     Check if testing a single threshold or a range of values.
      if(nmval .eq. 1) then

c        Single threshold.
         if(hicut .gt. locut) then
           gamma = hicut
           alpha = locut
           flipped = .false.
         else
           gamma = locut
           alpha = hicut
           flipped = .true.
         end if
          beta = midpt

c        Check for value beyond function range.   
         if(.not. flipped .and. val .gt. gamma) then
            c = 1.0

         else if(.not. flipped .and. val .lt. alpha) then
            c = 0.0

         else if(flipped .and. val .gt. gamma) then
            c = 0.0

         else if(flipped .and. val .lt. alpha) then
            c = 1.0

          else

c           Value is within the range of the function.
            if(val .le. beta) then

               range = 2.0 * (beta - alpha)
               s1 = (val - alpha) / range
               if(.not.flipped) c = coeff * s1**power
               if(flipped) c = 1.0 - (coeff * s1**power)

            else

               range = 2.0 * (beta - gamma)
               s1 = (val - gamma) / range
               if(.not.flipped) c = 1.0 - (coeff * s1**power)
               if(flipped) c = coeff * s1**power

            end if

         end if

      end if

c     Force confidence values to be between 0 and 1.              
      if(c .gt. 1.0) c = 1.0
      if(c .lt. 0.0) c = 0.0

      conflev = c

      return
      end
