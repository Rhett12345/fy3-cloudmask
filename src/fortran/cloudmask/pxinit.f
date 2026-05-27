      subroutine pxinit(testbits,qa_bits,precip_water,vza,
     +                  sfctmp,pmsl,u_wind,v_wind,plat,plon,lsf,polar,
     +                  day,night,land,water,coast,snglnt,visusd,
     +                  vrused,snow,ice,desert,bad_value,bad_geo,
     +                  uniform,shadow,smoke,cirrus_ir,cirrus_vis,
     +                  nmtests,nbands,nbad_1km,nbad_250,hi_elev,
     +                  antarctic,sh_ocean,sg_bad_data,map_ice,
     +                  map_snow,sh_lake)

      implicit none
      save

c ... Common statement for debug purposes
      !common / bug / debug, h_output

c----------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for initializing variables used in processing
c     individual MAS pixels.
c
c!Input parameters:
c testbits      Cloud Mask bit flag holder
c qa_bits       Byte array containing qa bit results
c geo_flag      Integer array containing geolocation good/bad flags
c               1-lat,2-lon,3-szen,4-vzen,5-rel_angle
c precip_water  value of precipitable water (g/cm2) at pixel location
c vza           Viewing zenith angle at pixel location
c sfctmp        Model surface temperate for current pixel
c pmsl          Model mean sea level pressure for current pixel
c u_wind        Model u wind component for current pixel
c v_wind        Model v wind component for current pixel
c plat          Pixel latitude value
c plon          Pixel longitude value
c lsf           Numeric land/sea flag (integer)
c polar         Logical variable flagging polar scenes
c day           Logical variable flagging day scenes
c night         Logical variable flagging night scenes
c land          Logical variable flagging land scenes
c water         Logical variable flagging water scenes
c coast         Logical variable flagging coastal scenes
c snglnt        Logical variable flagging sunglint contaminated scenes  
c visusd        Logical variable flagging scenes where visible
c               data were used
c vrused        Logical variable flagging scenes where vis. ratio test
c               can be applied
c snow          Logical variable flagging snow background scenes
c ice           Logical variable flagging ice background scenes
c desert        Logical variable flagging desert background scenes
c bad_value     Logical variable flagging bad input radiances or 
c               reflectances.
c bad_geo       Logical variable flagging bad lat/long data
c uniform       Logical variable flagging uniform scenes
c               (Places where all pixels in context are similar)
c shadow	Logical variable flagging shadow contaminated scenes
c smoke		Logical variable flagging smoke contaminated scenes
c cirrus_ir     Logical variable flagging thin cirrus contaminated 
c 		scenes in the infrared
c cirrus_vis    Logical variable flagging thin cirrus contaminated 
c 		scenes in the visible
c nmtests       Number of tests applied to this pixel
c nbands        Number of bands successfully read for this pixel
c nbad_1km      Number of bands with bad data for this pixel (1km)
c nbad_250      Number of bands with bad data for this pixel (250 
c hi_elev       Logical flag indicating elevation > 2000 meters
c antarctic     Logical flag indicating regions south of -60 latitude
c sh_ocean      Logical flag indicating shallow ocean
c sg_bad_data   Logical flag indicating bad data in one or more channels
c               needed for sun-glint processing
c map_ice       Logical flag indicating ice background (from ancillary
c               data)
c map_snow      Logical flag indicating snow background (from ancillary
c               data)
c sh_lake       Logical flag indicating shallow inland lakes
c
c!Output Parameters:
c None.
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

      !include 'global.inc'

c     scalar arguments
      real vza,precip_water,plat,plon,sfctmp,pmsl,u_wind,v_wind
      integer lsf,nmtests,nbands,nbad_1km,nbad_250,geo_flag(5)
      logical polar,day,night,land,water,snglnt,visusd,snow,ice,uniform,
     +        bad_value,shadow,coast,desert,vrused,smoke,cirrus_vis,
     +        cirrus_ir,bad_geo,hi_elev,antarctic,sh_ocean,sg_bad_data,
     +        map_ice,map_snow,sh_lake

c     array arguments
      byte testbits(6),qa_bits(10)

c     local scalars
      integer i !,debug,h_output

c     external subroutines
      external set_bit

c ... debug statement ............................................
      !if (debug .gt. 0) then
      !  write(h_output,'(10x/,''Subroutine pxinit.f '',/)')
      !endif
c ...............................................................
       
c     First initialize the bit holder
c     Initialize output "bit" array.  All bits are initialized to
c     0 with the exception of the thin cirrus, shadow and non-cloud
c     obstruction flags, which are set to 1.
c     Only initialize first 7 bits of qa_bits each time.  The ancillary
c     data bits are filled once for each granule.
      do i = 1 , 7
        if (i .le. 6) testbits(i) = 0 
        qa_bits(i) = 0
      end do

c ... Initialize the geolocation bad data flag holder
c     do i = 1 , 5
c        geo_flag(i) = 0
c     end do

c     NCO bit
      call set_bit(testbits,8)
c     Thin cirrus solar bit
      call set_bit(testbits,9)
c     Shadow bit
      call set_bit(testbits,10)
c     Thin cirrus infrared bit
      call set_bit(testbits,11)
c     Suspended dust bit
      call set_bit(testbits,28)

c     Initialize precipitable water holder
      precip_water = 0.0

c     Initialize pixel latitude and longitude holders
      plat = -999.0
      plon = -999.0

c     Initialize land/sea flag 
      lsf = 0
c     Initialize spectral test counter
      nmtests = 0
c     Initialize bands read counter
      nbands = 0
c     Initialize 1km bad channel data counter
      nbad_1km = 0
c     Initialize 250m bad channel data counter
      nbad_250 = 0

c     Initialize viewing zenith angle
      vza = 0.0

c     Initialize surface temperature holder
      sfctmp = 0.0
c     Initialize mean sea level pressure holder
      pmsl = 0.0
c     Initialize u_wind component holder
      u_wind = 0.0
c     Initialize v_wind component holder
      v_wind = 0.0

c     Initialize logical variables which will decide what processing
c     path to take
      polar = .false.
      day = .false.
      night = .false.
      land = .false.
      water = .false.
      coast = .false.
      snglnt = .false.
      visusd = .true.
      vrused = .true.
      snow = .false.
      ice = .false.
      desert = .false.
      bad_value = .false.
      uniform = .true.
      shadow = .false.
      smoke = .false.
      cirrus_ir = .false.
      cirrus_vis = .false.
      bad_geo = .false.
      hi_elev = .false.
      antarctic = .false.
      sh_ocean = .false.
      sh_lake = .false.
      sg_bad_data = .false.
      map_ice = .false.
      map_snow = .false.

      return
      end
