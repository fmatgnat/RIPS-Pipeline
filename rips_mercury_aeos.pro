pro RIPS_Mercury_AEOS, part

imaging_statsec = '[183 : 809, 39 : 444]'
spectra_statsec = '[ 99 : 990,502 : 978]'

if part eq 0 then begin ;break into managable data cubes
    cube = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Mercury - 500 kinetic 1 sec scan - Na 3A + ND imaging - Na order (2).fits', 0, header, /unsigned, /silent )
    imaging_cube = reform(cube[183 : 809, 39 : 444, *])
    spectra_cube = reform(cube[ 99 : 990,502 : 978, *])
    MWRFITS, imaging_cube, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\imaging_cube.fits', header, /create, /silent
    MWRFITS, spectra_cube, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\spectra_cube.fits', header, /create, /silent
    Sky_cube = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Sky flats - Na 3A + ND imaging - Na order.fits', 0, header, /unsigned, /silent )
    imaging_sky_cube = reform(sky_cube[183 : 809, 39 : 444, *])
    spectra_sky_cube = reform(sky_cube[ 99 : 990,502 : 978, *])
    MWRFITS, imaging_sky_cube, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\imaging_sky_cube.fits', header, /create, /silent
    MWRFITS, spectra_sky_cube, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\spectra_sky_cube.fits', header, /create, /silent
endif

if part eq 1 then begin ;Find the centroids by cross-correlation with the last image
    imaging_cube = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\imaging_cube.fits', 0, header, /fscale, /silent )
    frame = reform(imaging_cube[*,*,0])
    s = size(imaging_cube, /dimensions) 
    aligned_imaging_cube = fltarr(s) 
    window, 0, xs = s[0], ys = s[1]

    Dark = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 21 HST\Dark - 200 sec - Na 3A + ND imaging - Na order.fits', 0, header, /fscale, /silent )
    dark = dark[183 : 809, 39 : 444]
    dark = sigma_filter( dark, 5, N_sigma=5)
    flat = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 21 HST\Tungsten lamp flat - 200 sec - Na 3A + ND imaging - Na order_Autosave_recovered.fits', 0, header, /fscale ) 
    flat = (flat[183 : 809, 39 : 444] - dark) / mean(flat[183 : 809, 39 : 444]) 
    flat[WHERE(flat lt .01, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    flat[WHERE(flat gt 2.0, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    flat = smart_shift(flat, 1.25, .5,  /interp) ;hack
    gooddata = where(Finite(flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then flat[baddata] = interpol(flat[gooddata], gooddata, baddata)

    reference = (reform(imaging_cube[*,*,499]) - dark)
    shift_array = intarr(500,2)
    for i = 0, s[2]-1 do begin
        frame = (reform(imaging_cube[*,*,i]) - dark)
        CORREL_OPTIMIZE, reference, frame, xoffset_optimum, yoffset_optimum, NUMPIX = 150
        shift_array[i,*] = [xoffset_optimum, yoffset_optimum]
    endfor
    shift_array[21:25,0] = -221                            ;hack
    shift_array = smooth(shift_array, [5,1], /EDGE_MIRROR) ;hack
    for i = 0, s[2]-1 do begin
        frame = (reform(imaging_cube[*,*,i]) - dark) / flat ;make this into a movie!
        aligned_imaging_cube[*,*,i] = shift(frame, [shift_array[i,*]])
        tv, bytscl(shift(frame, [shift_array[i,*]]), 0, 5000) 
        wait, 0.05
    endfor    
    save, shift_array, aligned_imaging_cube, filename = 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\shift_array.sav'
endif

if part eq 2 then begin ;Isolate the sodium emission in every frame
    spectra_cube = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\spectra_cube.fits', 0, header, /fscale, /silent )
    frame = reform(spectra_cube[*,*,0])
    s = size(spectra_cube, /dimensions) 
    window, 0, xs = s[0], ys = s[1]
    tv, bytscl(Frame)
    
    ;For this data there is no dark or bias calibration frames
    ;Therefore, use the conventional mode date from the following night
    
    Dark = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 21 HST\Dark - 200 sec - Na 3A + ND imaging - Na order.fits', 0, header, /fscale, /silent )
    dark = dark[ 99 : 990,502 : 978]
    dark = sigma_filter( dark, 5, N_sigma=5)

    ;For this data there is no reference spectrum without sodium
    ;The sky spectrum won't work, the lunar spectrum the following night has a wider slit
    ;Therefore, let's hack our way through using Mercury itself.
    
    reference = mean(spectra_cube, dimension = 3) - (dark)
    window, 1, xs = s[0], ys = s[1]
    tv, bytscl(reference)
    ref_spectrum = total(reference[*,120:220], 2)    
    amplifer_correction = mean(reference[*,426:476], dimension = 2)
    amplifier_contribution = rebin(amplifer_correction, s[0], s[1])
    dark = dark + amplifier_contribution ; correct for the fact that the dark correction from conventional mode is clearly insufficient for EMCCD mode
    
    ;get a spectral flat
    flat = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 21 HST\Tungsten lamp flat - 200 sec - Na 3A + ND imaging - Na order_Autosave_recovered.fits', 0, header, /fscale ) 
    flat = (flat[ 99 : 990,502 : 978] - dark) 
    flat = flat / median(flat) 
    
    reference = (mean(spectra_cube, dimension = 3) - (dark)) / flat
    window, 1, xs = s[0], ys = s[1]
    tv, bytscl(reference, 0, 100) 
    ref_spectrum = total(reference[*,120:220], 2)    
    window, 2, xs = s[0], ys = s[1]
    cgplot, ref_spectrum, color = 'black', xr = [250,600]
    dummy = ref_spectrum
    dummy[286:297] = !values.f_Nan ;remove Mercury D2 emission
    dummy[528:538] = !values.f_Nan ;reject unusual counts for centroid
    gooddata = where(Finite(dummy), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then dummy[baddata] = interpol(dummy[gooddata], gooddata, baddata)
    cgplot, dummy, color = 'red', /overplot
    ref_spectrum = rebin(dummy, s[0], s[1]) ; this is what a flat-fielded lunar drift scan reference spectrum might have looked like of they had been given the time to take one with this night  
    
    ;now scale and subtract 
    window, 2, xs = s[0], ys = s[1]
    window, 3, xs = s[0], ys = s[1]
    window, 4, xs = s[0], ys = s[1]
    exosphere_spectra_cube = fltarr( size(spectra_cube, /dimensions) )
    for i = 0, s[2]-1 do begin
      frame = reform(spectra_cube[*,*,i]) - dark ;dark subtract the spectrum
      frame = frame / flat                       ;flatfield the spectrum
      
      illum = frame / ref_spectrum    ;scale the illumination against Mercury's disk
      illum[528:538, *] = !values.f_Nan  ;carefully avoid sodium emissions when scaling to the disk
      illum[286:297, *] = !values.f_Nan  ;carefully avoid sodium emissions when scaling to the disk
      
      illum_along_slit = mean(illum, dimension=1, /nan)
      illum_along_slit = MEDSMOOTH( illum_along_slit, 5 ) ;helps with hot pixels from bad flat-fielding where the slit has some dust
      This_Frames_Fake_Mercury_Spectrum = ref_spectrum * rebin(transpose(illum_along_slit), s[0], s[1])
      Just_exosphere = Frame - This_Frames_Fake_Mercury_Spectrum
      wset, 2
      tv, bytscl(frame, 0, 200)
      wset, 3
      tv, bytscl(This_Frames_Fake_Mercury_Spectrum, 0, 200)
      wset, 4
      tv, bytscl(Just_exosphere, 0, 200)
      wait, 0.02
      exosphere_spectra_cube[ *, *, i] = Just_exosphere
    endfor
    save, exosphere_spectra_cube, filename = 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\exosphere_spectra_cube.sav'
endif
if part eq 3 then begin; extract and place the 1D spectra. 
  restore, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\exosphere_spectra_cube.sav'
  s = size(exosphere_spectra_cube, /dimensions) 
  
  ;Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
  dispersion = (5895.92424 - 5889.95095) / (632. - 390.) ;Angstroms per pixel
  dispersion_Velocity = 299792.*dispersion / 5889.95095  ;Km/s per pixel at Na D2

   ;-------------------------------------------Extraction----------------------------------------------------------- 

    ;run MPFIT once to find the D2 line center
      img = total(exosphere_spectra_cube, 3)
      search = 16           ;pixels to search over
      expected_pixel = 291  ;Rough D2 pixel location
      D2_height = fltarr(n_elements(img[0,*]))
      D2_center = fltarr(n_elements(img[0,*]))
      D2_width  = fltarr(n_elements(img[0,*]))
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[expected_pixel-search:expected_pixel+search, i], a, STATUS = STATUS, /positive) 
        if STATUS ge 1 then D2_height[i] = A[0] else D2_height[i] = !values.F_nan
        if STATUS ge 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
        if STATUS ge 1 then D2_width[i] = A[2] else D2_width[i] = !values.F_nan
      endfor
      y = findgen(n_elements(img[0,*]))
      height = smooth(D2_height, 5, /edge_mirror) / s[2]
      height[0:25] = 1. & height[280:*] = 1.
      real = where(finite(D2_center), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
      location = poly(findgen(n_elements(img[0,*])), coeff) + expected_pixel-search 
      real = where(finite(D2_width), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_width[real], 2)
      width = poly(findgen(n_elements(img[0,*])), coeff) 
      dummy = img
      dummy[location, y] = max(dummy) 
      tv, bytscl(dummy, -100, 200)

    ;setup MPFIT Gaussian parameters
      guess = [80.,location[0],1.5,0.]
      A = guess
      parinfo = replicate( {value:0.D, fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 3)
      parinfo[0].limited = 1b                                 ;limit amplitude 
      parinfo[0].limits = [0.1, 1000. ]                       ;positive only
      parinfo[1].limited = 1b                                 ;limit sigma   
      parinfo[1].limits = [14., 18.]                          ;limit sigma width in pixels      
      parinfo[2].limited = 1b                                 ;limit sigma   
      parinfo[2].limits = [1., 5.]                            ;limit sigma width in pixels
  
  brightness = fltarr(s[2], s[1]) & err_brightness = fltarr(s[2], s[1])
  linewidth = fltarr(s[2], s[1]) & err_linewidth = fltarr(s[2], s[1])
  for n = 0, s[2]-1 do begin
      img = exosphere_spectra_cube[*,*,n]
      err_img = sqrt( abs(img) )     
      for i = 0, s[1]-1 do begin
               a[0] = height[i]
               a[1] = location[i] - expected_pixel + search
               a[2] = Width[i]
               
               ;cgplot, findgen(search*2. + 1), img[location[i]-search:location[i]+search, i]
               result = mpfitpeak(findgen(search*2. + 1), img[location[i]-search:location[i]+search, i], A, PERROR = Err_A, $
                                  error = err_img[location[i]-search:location[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
               ;cgoplot, findgen(search*2. + 1), result, color = 'red'
               if STATUS ge 1 then begin
                  brightness[n,i] = A[0]*A[2]*SQRT(2*!DPI) 
                  err_brightness[n,i] = brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
                  linewidth[n,i] = dispersion_Velocity*sqrt((2.355*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
                  err_linewidth[n,i] = dispersion_Velocity*2.355*Err_A[2]
               endif else begin
                  brightness[n,i] = !values.F_nan
                  linewidth[n,i] = !values.F_nan
               endelse 
      endfor
  endfor  
  tv, bytscl(brightness, 0, 1000)
  save, brightness, linewidth, filename = 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\brightness.sav'
endif  
if part eq 4 then begin ;put everything together and build an image of the exosphere. 
   restore, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\brightness.sav'
   restore, 'H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\shift_array.sav'
   imaging_cube = MRDFITS('H:\DATA\RIPS\2018 - AEOS with RIPS\June 20 HST\Processed\imaging_cube.fits', 0, header, /fscale, /silent )
   s = size(imaging_cube, /dimensions)
   
   brightness = congrid(brightness, s[2], s[1]) ;HACK NEED PROPER SCALING OF THE SPECTRUM TO IMAGING PLATESCALES
   home = [341,238] ;Mercury centroid of the last frame number 500, to which all the other frames are aligned
   
   big_cube = fltarr(s[0], 2*s[1], s[2])
   window, 0, xs = s[0], ys = 2*s[1]
   for i = 0, s[2]-1 do begin
        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        big_frame[*, s[1]:*] = rotate(empty_frame, 7) * 5. ;the factor of 5 is a HACK FOR NOW
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i])
        tv, bytscl(big_frame, 0, 1000) 
        big_cube[*,*,i] = big_frame
        wait, 0.02   
        ;stop
   endfor     
   window, 2, xs = s[0], ys = 2*s[1]
   test = median(big_cube, dimension = 3, /even)
   ;tv, bytscl(smooth(test,[3,3], /NaN), 0, 10000)
   tv, bytscl(test, 0, 10000)  
   stop
endif  
end

