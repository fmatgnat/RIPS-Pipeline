pro RIPS_Mercury_Perkins_v03, PART=part, NIGHT=night

; *********************************************************************************************************************
; *********************************************************************************************************************
; Routine to extract Na emission from spectral scans over Mercury's disk with Perkins/RIPS in order to create a 2-D
; image of Na above and near Mercury's disk.
; 
; The program is split into 4 "parts":
;   0 = break kinetic series into manageable datacube
;       input:   kinetic series FITS file(s)
;       output:  RIPS imaging (x,y,t) and spectral (wl,y,t) datacubes 
;   1 = find Mercury centroids by cross-corellation with previous image
;       input:   imaging datacube from Part 0
;       output:  co-aligned imaging datacube (x,y,t)
;   2 = isolate sodium emission in every frame
;       input:   spectral datacube from Part 0
;       output:  the calibrated spectral datacube (wl,y,t) for Na exosphere
;   3 = extract Na signal and place into 1D spectra
;       input:   spectral datacube from Part 2 (wl,y,t)
;       output:  Na brightness (wl,y) and linewidth (wl,y)
;   4 = put everything together and build exosphere image
;       input:   Na brightness and linewidth from Part 3 and imaging cube from Part 0
;       output:  none yet, a combined image called "test" for now...
;  99 = do all of 0-4
;
; Input variables
;   part           - the "part" of the program to execute
;   night          - the night of RIPS data to use; only '15' is viable from the March 2018 run, but we've kept
;                    the option to specify night in case future runs are more successful!  There are numerous 
;                    sub-variables associated with each night: e.g., mercury_dir, dark_dir, flat_dir, sky_dir, 
;                    Mercury_file, sky_file, dark_file, flat_file, Na_D_rngs, ims, sps
;   thresh         - hotpixel removal threshold (sigma above background; Default=12; acre routine)
;   width          - median smoothing width to get local background reference (for pixels above thresh; Default=5)
; 
; Created:  18 Oct 2018 - copied from RIPS_Mercury_AEOS_v03.pro to start (CS and LM)
; Edits:    25 Oct 2018 - renamed to "v02", cleaned up with the plan to leave v02 as a minimal, yet still hard-coded
;                         version, and then to improve the automation in v03 (plus add additional improvements) (LM)    
;           30 Oct 2018 - changed approach to reference images, we now use either a pre-generated reference from
;                         RIPS_Mercury_Perkins_reference.pro or from a previous run (LM)
;            1 Nov 2018 - renamed to "v03" as part of integration with GitHub (LM) 
; *********************************************************************************************************************
; =====================================================================================================================
; Define variables
; =====================================================================================================================
SetDefaultValue, part, 4
SetDefaultValue, night, '15'                                              ; the night of the run; default '15' (only kept for possible future runs with multiple nights)
SetDefaultValue, framerate, 30                                            ; default framerate for video output
SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh) 
SetDefaultValue, minfac, 0.9                                              ; maximum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, maxfac, 1.4                                              ; minimum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, dfac, 0.01                                               ; factor increment for minfac->maxfac
SetDefaultValue, ct, 20                                                   ; default color table
SetDefaultValue, max_spec, 1550.  ;~1900 or 1750 when smoothed             ; partial hack --> the max value of the spectral Na cube after averaging (used for "prettier" movies)
SetDefaultValue, max_img, 1.6e6                                           ; partial hack --> the max value of the image cube after coadding (used for "prettier" movies)
SetDefaultValue, contour_outline, 0.3                                     ; contour level (between 0-1) of the outer Mercury "disk" (this is plotted on top of Na images to give a sense of where Mercury is for pics and movies)
SetDefaultValue, big_frame_max, 1500;                                     ; partial hack --> the max value of big_frame after scaling (used for "prettier" movies)
SetDefaultValue, img_extraction, [145,295,85,235]                         ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
SetDefaultValue, spec_extraction, [145,295,77,227]                        ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
SetDefaultValue, movie_scale, 4                                           ; scale the x-y dimensions of the extracted "Na image" and "Mercury disk" arrays for movies by this value (used for bigger and therefore "prettier" movies)
SetDefaultValue, do_rotate, 0                                             ; 0=no rotation of image, 1=rotate Na image (inverted, =5 -- see below)
SetDefaultValue, do_median, 0                                             ; 0=use mean when constructing Na Mercury image, 1=use median instead
SetDefaultValue, do_smooth, 1                                             ; set to 1 to smooth reference and frame images before coaligning
SetDefaultValue, do_plot, 1                                               ; set to 1 to plot progress (slightly slows things down)
SetDefaultValue, do_bad, 1                                                ; set to 1 to also make a couple "bad frame" images
SetDefaultValue, smooth_width, 6                                          ; width of smoothing kernel if do_smooth = 1
SetDefaultValue, stddev_cutoff, 2.                                        ; only use Mercury frames with stddev > (median(std_devs) - [stddev_cutoff]_sigma), these are defined as "good" 

mercury_dir   = 'Y:\obs_18\Perkins_RIPS_March\15\renamed\'                ; directory with Mercury FITS file(s)
dark_dir      = 'Y:\obs_18\Perkins_RIPS_March\14\Carl_Keep\'              ; directory with dark file
flat_dir      = 'Y:\obs_18\Perkins_RIPS_March\15\'                        ; directory with flat file
sky_dir       = 'Y:\obs_18\Perkins_RIPS_March\15\'                        ; directory with sky file
Mercury_file  = 'Mercury_Na_' + strfix(['791','792','793','794','796','798','800','803']) + '.fits'  ; Mercury kinetic series to use (note that these have been renamed from the "RIPS1_..." defaults)
Mercury_file  = 'Mercury_Na_796.fits' ; TEMP for faster testing
ref_file      = 'Mercury_ref_800frames.fits'                              ; a reference image of Mercury (combined from "good" Mercury data) -- to be used for finding alignments below
sky_file      = 'RIPS1_Thu Mar 15 2018_01.12.12_785.fits'                 ; sky frame
dark_file     = 'Dark.fits'                                               ; dark frame
flat_file     = 'Flat_NaSpectra_624slitwidth_NaND1imaging1.fits'          ; flat frame
Na_D_rng      = 20                                                        ; number of pixels around Na D Mercury emission to avoid (in illum scaling for part 2)
do_realign    = 0                                                         ; 0=run part 1 as normal, 1=use previous cross-correlation as master template
offset_maxes  = [80,40]                                                   ; if the absolute value of the [x,y] offset is larger than these values then that frame is marked "bad", and we try an alternative approach
ims           = [283,759,39,374]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec")
sps           = [99,899,502,941]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec")
iref_factor   = 1.                                                        ; factor to scale reference spectrum by
platescale_fac = (275./350.) / (250./300.)                                ; purported imaging to spectral platescale ratio
count_factor  = 1.e4                                                      ; this is the scaling factor to apply to the white light Mercury image (i.e., the slit covers some regions of the disk, so we scale based on total slit coverage at each location)

outdir        =  Mercury_dir + 'Processed\'                               ; output directory
nfiles        =  n_elements(Mercury_file)                                 ; # of different Mercury kinetic frames

; =====================================================================================================================
; Find display size, set window defaults, reset plotting variables just in case
; =====================================================================================================================
DEVICE, GET_SCREEN_SIZE = ss                                              ; find screen size
winpos_x      = ss(0)/2                                                   ; default initial x window location
winpos_y      = 0                                                         ; default initial y window location
!p.multi      = 0.
!p.position   = 0.
cgloadct, ct



; =====================================================================================================================
; Part 0 : break kinetic series into manageable datacube
; =====================================================================================================================
if part eq 0 or part eq 99 then begin 
    icube = MRDFITS(Mercury_dir + Mercury_file(0), 0, header, /unsigned, /silent )
    n_frames_per_file = sxpar(header, 'NUMKIN')                           ; # of frames for each "Mercury_file"
    s         = size(icube, /dimensions) 
    cube = fltarr(s(0),s(1),s(2)*nfiles)                                  ; define cube for multiple kinetic series
    for ifile = 0, nfiles-1 do begin                                      ; loop over kinetic series and combine
      icube   = MRDFITS(Mercury_dir + Mercury_file(ifile), 0, header, /unsigned, /silent )
      cube(*,*,s(2)*ifile:s(2)*ifile+n_frames_per_file-1) = icube
    endfor
    imaging_cube = reform(cube[ims(0):ims(1),ims(2):ims(3), *])           ; extract the IMAGING portion of the kinetic series and combine
    spectra_cube = reform(cube[sps(0):sps(1),sps(2):sps(3), *])           ; extract the SPECTRAL portion of the kinetic series and combine
    MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent
    MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent
    Sky_cube  = MRDFITS(sky_dir + sky_file, 0, sky_header, /unsigned, /silent )
    imaging_sky_cube = reform(sky_cube[ims(0):ims(1),ims(2):ims(3), *])
    spectra_sky_cube = reform(sky_cube[sps(0):sps(1),sps(2):sps(3), *])
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', sky_header, /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', sky_header, /create, /silent
    beep
endif



; =====================================================================================================================
; Part 1 : find Mercury centroids by cross-corellation with previous image
; =====================================================================================================================
if part eq 1 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
    if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )  ; don't need to read in if we're continuing
    if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
    frame     = reform(imaging_cube[*,*,0])                               ; individual frame from kinetic series
    s         = size(imaging_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
    ss        = size(spectra_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame
    
    ; Define variables based on frame sizes
    shift_array = intarr(s[2],2)                                          ; array of x,y shift values for aligning images
    std_devs  = fltarr(s[2])                                              ; standard deviations of imaging frames (TEMPORARILY used to define image quality)
    
    ; Set up windows based on image sizes
    if do_plot then begin
        window, 0, xs = s[0], ys = s[1], xpos=winpos_x,         ypos=winpos_y
        window, 1, xs = s[0], ys = s[1], xpos=winpos_x,         ypos=winpos_y+s[1]+40
        window, 2, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20, ypos=0
        window, 3, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20, ypos=s[1]+40
    endif

    iDark = MRDFITS(dark_dir + dark_file, 0, dark_header, /fscale, /silent )   ; read in dark
    iflat = MRDFITS(flat_dir + flat_file, 0, flat_header, /fscale )            ; read flat, extract imaging portion, scale
    dark      = idark[ims(0):ims(1),ims(2):ims(3)]                        ; extract and normalize to counts/sec
    dark      = sigma_filter( dark, width, N_sigma=thresh)                ; get rid of any hot pixels
    iflat     = iflat[ims(0):ims(1),ims(2):ims(3)]                        ; extract 
    flat      = iflat / max(iflat)
    flat      = flat + (1. - median(flat))                                ; this normalizes such that the median of the flat is 1
    flat[WHERE(flat lt .01, /NULL)] = !values.f_Nan                       ; reject unusual counts for centroid
    flat[WHERE(flat gt 2.0, /NULL)] = !values.f_Nan                       ; reject unusual counts for centroid
    gooddata = where(Finite(flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then flat[baddata] = interpol(flat[gooddata], gooddata, baddata)
    if do_realign then begin
        restore, outdir+'shift_array.sav'                                 ; contains shift_array, std_devs, and aligned_imaging_cube
        shift_array_old = shift_array                                     ; keep old shift_array for comparison
        reference = total(aligned_imaging_cube,3)                         ; create "master Mercury" template for cross-correlation
    endif else begin
        ireference = readfits(Mercury_dir + ref_file, ref_header, /silent)
        acre, reform(ireference(*,*,0)), reference, thresh, width
    endelse

    ; Find indices of the pixels under the slit in the imaging channel
    junk = min(total(flat, 2), slit_center)                         ; roughly find the slit center
    h = histogram(flat, binsize = .05, REVERSE_INDICES=ri)          ; bin the flat into 0.05 bins
    slit_indices = ri[ri[0]:ri[9]-1]                               ; find inidices of pixels in the lowest few bins of the flat's histogram 
    slit_indices = array_indices(flat, slit_indices)              ; convert into x & y indices
    keep = where(abs(slit_indices[0,*] - slit_center) lt 10, /Null); keep only inidicees within 10 pixels of slit center       
    slit_indices = slit_indices[*,keep]
    dummy = flat
    dummy[ slit_indices[0,*], slit_indices[1,*] ] = max(flat)     ; show the location of pixels fully under the slit.
    if do_plot then cgimage, hist_equal(dummy), title='flat'

    ; Find slit width -- NOTE that we need to replace these next few lines with Carl's updated slit method (30 Oct 2018, LM)
;    xpoly = poly_fit(indgen(n_elements(total(flat,2))),total(flat,2),3,yfit=yfit)   ; fit the 3D trend across the flat
;    stddev_yfit = stddev(total(flat,2) - yfit)                            ; find standard deviation of detrended x behavior
;    xslit_pixels = where( (total(flat,2) - yfit) le -stddev_yfit )        ; find x pixels more than 1 stddev below mean
;    slit_width = n_elements(xslit_pixels)/2 - 1                           ; define slit width 
;    xslit     = (where( total(flat,2) eq min(total(flat,2)) ))[0]             ; the "middle" of the slit
;    if do_plot then begin
;        cgimage, hist_equal(flat), /axes, title='flat'
;        for i = -1, 1, 2 do cgoplot, [1,1]*(xslit+i*slit_width), !Y.Crange
;        wait, 2.
;    endif

    ; Begin loop through each frame, calculate alignment
    wset, 1
    cgimage, reference, /axes, title='Reference image for alignment'
    for i = 0, s[2]-1 do begin
        iframe = (reform(imaging_cube[*,*,i]) - dark) / flat              ; "raw" imaging frame
        iframe[ slit_indices[0,*], slit_indices[1,*] ] = !values.F_Nan    ; define slit pixels as NaN
        acre, iframe, frame, thresh, width                                ; clean any hot pixels
;        frame[xslit-slit_width:xslit+slit_width,*] = !values.F_Nan        ; define slit pixels as NaN
        fill_missing_rand, frame, !values.F_Nan, 1                        ; replace slit pixels based on nearby pixels
        if do_smooth then begin                                           ; set up bandpass filter to aid image alignment
          frame_bandpass = bandpass_filter(smooth(frame,smooth_width,/edge_truncate), 0., 0.15, /butterworth)
          ref_bandpass = bandpass_filter(smooth(reference,smooth_width,/edge_truncate), 0., 0.15, /butterworth)
        endif else begin
          frame_bandpass = bandpass_filter(frame, 0., 0.15, /butterworth)
          ref_bandpass = bandpass_filter(reference, 0., 0.15, /butterworth)
        endelse
        CORREL_OPTIMIZE, ref_bandpass, frame_bandpass, xoffset_optimum, yoffset_optimum, NUMPIX=150
        shift_array[i,*] = [xoffset_optimum, yoffset_optimum]             ; track required alignment offsets

        ; If CORREL fails then we need to either define these frames as "bad" (i.e., set shift_array[i,*] = -666) or else try a new method.
        xcenters = [ (where(total(ref_bandpass^2.,2) eq max(total(ref_bandpass^2.,2))))[0], $
                     (where(total(frame_bandpass^2.,2) eq max(total(frame_bandpass^2.,2))))[0] ]
        ycenters = [ (where(total(ref_bandpass^2.,1) eq max(total(ref_bandpass^2.,1))))[0], $
                     (where(total(frame_bandpass^2.,1) eq max(total(frame_bandpass^2.,1))))[0] ]
        if abs(xoffset_optimum) gt offset_maxes[0] or abs(yoffset_optimum) gt offset_maxes[1] then begin  ; for now we try a new method
          xoffset_optimum = xcenters[0] - xcenters[1]
          yoffset_optimum = ycenters[0] - ycenters[1]
          shift_array[i,*] = [xoffset_optimum, yoffset_optimum]
          print, 'CORREL failed for frame ' + strfix(i) + ' -- new offsets (x,y): ' + strfix(xoffset_optimum) + ', ' + strfix(yoffset_optimum)
        endif else print, i, shift_array(i,0), shift_array(i,1)           ; update progress
        
        ; Plot progress to track quality of alignments
        if do_plot then begin
            ref_bkg  = mean(ref_bandpass(s[0]-100:*,*))
            xrange = fix(where(normalize(total(ref_bandpass-ref_bkg,2)) gt 0.1)) ; x range centered on Mercury
            yrange = fix(where(normalize(total(ref_bandpass-ref_bkg,1)) gt 0.1)) ; y range centered on Mercury
            wset, 0
            cgcontour, ref_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(ref_bandpass), title='initial bandpass comparison', levels=[0.2,0.4,0.6,0.8], /fill
            cgcontour, ref_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(ref_bandpass), title='initial bandpass comparison', levels=[0.2,0.4,0.6,0.8], /onimage
            cgcontour, frame_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(frame_bandpass), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='red'
            wset, 2
            cgcontour, ref_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(ref_bandpass), title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8], /fill
            cgcontour, ref_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(ref_bandpass), title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8], /onimage
            f = shift(frame_bandpass, [shift_array[i,*]])
            cgcontour, f(min(xrange):max(xrange),min(yrange):max(yrange))/max(f), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='red'
            wset, 1
            cgimage, reference(min(xrange):max(xrange),min(yrange):max(yrange)), /axes, title='raw reference'
            wset, 3
            f = shift(frame, [shift_array[i,*]])
            cgimage, f(min(xrange):max(xrange),min(yrange):max(yrange)), /axes, title='aligned raw frame'
            cgcontour, ref_bandpass(min(xrange):max(xrange),min(yrange):max(yrange))/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], /onimage
        endif ;do_plot==1
    endfor  ; end loop to determine shift_array values to align images

    good_frames = where(shift_array(*,0) ne -666)                         ; keep track of good frames (if applicable -- we won't use the bad ones later on)
    aligned_imaging_cube = fltarr(s[0],s[1],n_elements(good_frames))      ; imaging cube after everything is aligned
    aligned_spectra_cube = fltarr(ss[0],ss[1],n_elements(good_frames))    ; spectral cube after realignment (note that +y in the spatial dimension is -y in the spectral!)

    ; Now loop over good frames to align everything
    for i = 0, n_elements(good_frames)-1 do begin
        iframe = (reform(imaging_cube[*,*,i]) - dark) / flat              ; "raw" imaging frame
        acre, iframe, frame, thresh, width                                ; clean any hot pixels
        specframe = (reform(spectra_cube[*,*,i]))                         ; "raw" spectral frame -- missing specdark and specflat ??
        aligned_imaging_cube[*,*,i] = shift(frame, [shift_array[i,*]])    ; apply shift_array values to imaging frame
        aligned_spectra_cube[*,*,i] = shift(specframe, [0,-shift_array[i,1]])  ; note that we are taking the negative of the y shift from the image!
        img_frame = reform(frame(min(xrange):max(xrange),min(yrange):max(yrange))) ; extract just the image part of the frame
        std_devs(i) = stddev(img_frame)                                   ; higher std dev --> sharper image (TEMP, we can do better!)
    endfor    
;    save, shift_array, aligned_imaging_cube, aligned_spectra_cube, std_devs, xslit, slit_width, filename = outdir + 'shift_array.sav'
    save, shift_array, aligned_imaging_cube, aligned_spectra_cube, std_devs, slit_indices, filename = outdir + 'shift_array.sav'
    beep
endif



; =====================================================================================================================
; Part 2 : Isolate the sodium emission in every frame
; =====================================================================================================================
if part eq 2 or part eq 99 then begin 
    img_cube  = MRDFITS(outdir + 'imaging_cube.fits', 0, header2, /fscale, /silent )
    if part ne 99 then restore, outdir + 'shift_array.sav'                ; contains shift_array, std_devs, aligned_imaging_cube, aligned_spectra_cube, xslit, slit_width
    frame     = reform(aligned_spectra_cube[*,*,0])
    s         = size(aligned_spectra_cube, /dimensions) 
    si        = size(img_cube, /dimensions)
    if do_plot then begin                                                 ; set up window for plotting if needed
        window, 0, xpos=winpos_x,         ypos=winpos_y,           xs=s[0],  ys=s[1],  title='IDL 0 - FRAME'
        window, 1, xpos=winpos_x,         ypos=winpos_y+s[1]+40,   xs=s[0],  ys=s[1],  title='IDL 1 - SCALED REFERENCE'
        window, 2, xpos=winpos_x,         ypos=winpos_y+2*s[1]+80, xs=s[0],  ys=s[1],  title='IDL 2 - EXOSPHERE'
        window, 3, xpos=winpos_x+s[0]+20, ypos=winpos_y,           xs=si[0], ys=si[1], title='IDL 3 - FRAME'
    endif

    iDark     = MRDFITS(dark_dir + dark_file, 0, dark_header, /fscale, /silent )   ; read in dark
    iflat     = MRDFITS(flat_dir + flat_file, 0, flat_header, /fscale )            ; read flat, extract imaging portion, scale
    dark      = idark[sps(0):sps(1),sps(2):sps(3)]
    dark      = sigma_filter( dark, width, N_sigma=thresh)
    acre, dark, dark, thresh, width                                       ; try acre, still hot pixels in dark

    if night eq '15' then begin
        spectra_sky_cube = readfits(outdir + 'spectra_sky_cube.fits', header, /silent )
        ireference = spectra_sky_cube - dark
        acre, ireference, reference, thresh, width

        ; get a spectral flat
        iflat = (iflat[sps(0):sps(1),sps(2):sps(3)] - dark)
        iflat = iflat / median(iflat)
        acre, iflat, flat, thresh, width
        flat  = flat + (1. - median(flat))                                 ; this normalizes again such that median(flat) = 1

        reference = reference / flat                                      ; this is our dark and flat corrected spectral reference

        ; Now find the xoffset in the reference spectrum (likely due to slightly different grating angles)
        y_Mercury       = where( total(frame,1)-min(total(frame,1)) ge 0.5*(max(total(frame,1))-min(total(frame,1))) ) ; FWHM range of mercury
        frame_spectrum  = total(frame(*,y_Mercury),2)                     ; the spectrum for this frame, integrated over spatial extent of Mercury disk
        ref_spectrum    = total(reference(*,y_Mercury), 2)                ; reference spectrum, integrated over spatial extent of Mercury disk
        frame_D2_x      = ( where( frame_spectrum(0:s[0]/2) eq min(frame_spectrum(0:s[0]/2))) ) [0]    ; hack (we only want the D2 line, so we just check the 1st half of spectrum)
        ref_D2_x        = ( where( ref_spectrum(0:s[0]/2)   eq min(ref_spectrum(0:s[0]/2))   ) ) [0]   ; hack (we only want the D2 line, so we just check the 1st half of spectrum)

        reference       = shift(reference, frame_D2_x - ref_D2_x)         ; correct the reference based on the Na D alignment
        ref_spectrum    = shift(ref_spectrum, frame_D2_x - ref_D2_x)      ; and the same for the spectrum

        ; Find ranges of pixels surrounding Mercury Na D lines
        D1_x            = (where(frame_spectrum eq min(frame_spectrum)))[0] ; x pixel for first Na D line
        frame_temp      = frame_spectrum
        frame_temp(D1_x-Na_D_rng:D1_x+Na_D_rng) = mean(frame_spectrum)    ; get rid of Na emission peak to find next Na D line
        D2_x            = (where(frame_temp eq min(frame_temp)))[0]       ; x pixel for second Na D line
        if D1_x lt D2_x then Na_D_rngs = [D1_x-Na_D_rng, D1_x+Na_D_rng, D2_x-Na_D_rng, D2_x+Na_D_rng] $  ; define ranges based on whether first D found is D1 or D2
                        else Na_D_rngs = [D2_x-Na_D_rng, D2_x+Na_D_rng, D1_x-Na_D_rng, D1_x+Na_D_rng]

        exosphere_spectra_cube = fltarr( size(aligned_spectra_cube, /dimensions) )  ; variable to hold the Mercury spectral info

        ; Now scale reference spectrum and subtract from frame spectrum 
        for i = 0, s[2]-1 do begin                                        ; LOOP over each frame in the series
            frame       = reform(aligned_spectra_cube[*,*,i]) - dark      ; dark subtract the spectrum
            frame       = frame / flat                                    ; flatfield the spectrum
      
            ; Now find best scaling for ref_spectrum and subtract (based on solar reflected emission between D lines        
            illum = frame / reference                                     ; scale the illumination against Mercury's disk
            for iD = 0, 2, 2 do illum[Na_D_rngs(iD+0):Na_D_rngs(iD+1),*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
      
            illum_along_slit = MEDIAN(illum, dimension=1)                 
            illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )               ;helps with hot pixels from bad flat-fielding where the slit has some dust
            scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
            Just_exosphere = Frame - scaled_reference * iref_factor       ; at this point the spectral info along the spatial extent of Mercury should hopefully be indistinguishable from the background, leaving only exosphere emission

            ; Check that we've done a good job of minimizing spectrum between Na D lines
            between_D = smooth( just_exosphere(Na_D_rngs(1):Na_D_rngs(2),*), smooth_width, /edge_truncate )
            mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                              mean(between_D(*,10:min(Y_Mercury))), $             ; mean value below mercury's disk (avoid edge effects)
                              mean(between_D(*,max(Y_Mercury):s[1]-10)) ])        ; mean value above Mercury's disk (avoid edge effects)

            if mean_vals(0) gt 2.*mean(mean_vals(1:2)) then begin         ; find the best "ref_factor" if we haven't
                nfacs = (maxfac - minfac)/dfac + 1                        ; # of factors to consider
                factors = cgscalevector(indgen(nfacs),minfac,maxfac)
                fac_diffs = fltarr(nfacs)                                 ; keep track of the differences in mean_vals for each factor
                for ifac = 0, nfacs-1 do begin                            ; iterate through the factors to confirm the best match
                    Just_exosphere = frame - scaled_reference*factors(ifac)
                    between_D = smooth( just_exosphere(Na_D_rngs(1)+Na_D_rng:Na_D_rngs(2)-Na_D_rng,*), smooth_width )
                    mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                                      mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                                      mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
                    fac_diffs(ifac) = abs( mean_vals(0) - mean(mean_vals(1:2)) )            
                endfor ;ifac
                best_fac = factors(( where( fac_diffs eq min(fac_diffs) ) )[0])
                if best_fac eq maxfac or best_fac eq minfac then $
                  print, 'Frame ' + strfix(i) + ', best_fac = ' + strfix(best_fac) + ' -- END OF FAC RANGE!' $  ; print out so we know if min/maxfac are reasonable
                  else print, 'Frame ' + strfix(i) + ', best_fac = ' + strfix(best_fac)   ; print out so we know if min/maxfac are reasonable
                just_exosphere = frame - scaled_reference * best_fac
            endif
        
            if do_plot then begin
                wset, 0
                cgimage, bytscl(frame, 0, 200), /axes                          ; plot the original spectral frame
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange     ; show the Na exclusion zones
                wset, 1
                if do_smooth then cgimage, bytscl(scaled_reference, 0, 200), /axes, title='smoothed and scaled by ' + strfix(best_fac) $
                             else cgimage, bytscl(scaled_reference, 0, 200), /axes, title='scaled by ' + strfix(best_fac)
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange
                wset, 2
                cgimage, bytscl(Just_exosphere, 0, 200), /axes                 ; plot the Mercury Na emission after subtracting disk
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange
                wset, 3
                img_frame = shift(reform(img_cube(*,*,i)), [shift_array[i,*]]) ; show frame just for visual reference
                cgimage, bytscl( img_frame, 0, 3500 )
                wait, 0.02
                exosphere_spectra_cube[ *, *, i] = Just_exosphere              ; combine into cube for saving
                cgtext, 1, 1, 'Frame ' + strfix(1+i)
            endif ;do_plot==1
        endfor ;i (loop over frames)
        save, exosphere_spectra_cube, filename = outdir + 'exosphere_spectra_cube.sav'
    endif ;night=='15'

    beep
endif



; =====================================================================================================================
; Part 3 : extract and place the 1D spectra.
; =====================================================================================================================
if part eq 3 or part eq 99 then begin 
    if part ne 99 then restore, outdir + 'exosphere_spectra_cube.sav'
    s = size(exosphere_spectra_cube, /dimensions) 
  
    ;Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
    dispersion = (5895.92424 - 5889.95095) / (632. - 390.) ;Angstroms per pixel
    dispersion_Velocity = 299792.*dispersion / 5889.95095  ;Km/s per pixel at Na D2

    ; Find ranges of pixels surrounding Mercury Na D lines
    frame_spectrum = total(total(exosphere_spectra_cube,3),2)
    D1_x           = (where(frame_spectrum eq max(frame_spectrum)))[0]    ; x pixel for D1 line
    frame_temp     = frame_spectrum
    frame_temp(D1_x-Na_D_rng:D1_x+Na_D_rng) = mean(frame_spectrum)        ; get rid of Na emission peak to find next Na line
    D2_x           = (where(frame_temp eq max(frame_temp)))[0]            ; x pixel for D2 line
    Na_D_rngs      = [D1_x-Na_D_rng, D1_x+Na_D_rng, D2_x-Na_D_rng, D2_x+Na_D_rng]

    ;-------------------------------------------Extraction----------------------------------------------------------- 

    ; run MPFIT once to find the D2 line center
    img            = total(exosphere_spectra_cube, 3)
    search         = 16                                                   ; pixels to search over
    expected_pixel = mean(Na_D_rngs(0:1))                                 ; Rough D2 pixel location
    D2_height = fltarr(n_elements(img[0,*]))
    D2_center = fltarr(n_elements(img[0,*]))
    D2_width  = fltarr(n_elements(img[0,*]))
    for i = 0, n_elements(img[0,*]) - 1 do begin
      result = mpfitpeak(findgen(search*2. + 1), img[expected_pixel-search:expected_pixel+search, i], a, STATUS = STATUS, /positive) 
      if STATUS ge 1 then D2_height[i] = A[0] else D2_height[i] = !values.F_nan
      if STATUS ge 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
      if STATUS ge 1 then D2_width[i] = A[2] else D2_width[i] = !values.F_nan
    endfor
    y              = findgen(n_elements(img[0,*]))
    height         = smooth(D2_height, smooth_width, /edge_mirror) / s[2]
    height[0:25]   = 1. & height[280:*] = 1.
    real           = where(finite(D2_center), /NULL)
    COEFF          = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
    location       = poly(findgen(n_elements(img[0,*])), coeff) + expected_pixel-search 
    real           = where(finite(D2_width), /NULL)
    COEFF          = ROBUST_POLY_FIT(y[real], D2_width[real], 2)
    width          = poly(findgen(n_elements(img[0,*])), coeff) 
    dummy          = img
    dummy[location, y] = max(dummy) 
    tv, bytscl(dummy, -100, 200)

    ; setup MPFIT Gaussian parameters
    guess          = [80.,location[0],1.5,0.]
    A              = guess
    parinfo        = replicate( {value:0.D, fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 3)
    parinfo[0].limited = 1b                                               ;limit amplitude 
    parinfo[0].limits = [0.1, 1000. ]                                     ;positive only
    parinfo[1].limited = 1b                                               ;limit sigma   
    parinfo[1].limits = [14., 18.]                                        ;limit sigma width in pixels      
    parinfo[2].limited = 1b                                               ;limit sigma   
    parinfo[2].limits = [1., 5.]                                          ;limit sigma width in pixels
  
    brightness     = fltarr(s[2], s[1]) & err_brightness = fltarr(s[2], s[1])
    linewidth      = fltarr(s[2], s[1]) & err_linewidth = fltarr(s[2], s[1])

    for n = 0, s[2]-1 do begin
        if do_smooth then img = smooth( exosphere_spectra_cube[*,*,n], [2,smooth_width] ) $
                     else img = exosphere_spectra_cube[*,*,n]
        err_img    = sqrt( abs(img) )     
        for i = 0, s[1]-1 do begin
            a[0]   = height[i]
            a[1]   = location[i] - expected_pixel + search
            a[2]   = Width[i]
            result = mpfitpeak(findgen(search*2. + 1), img[location[i]-search:location[i]+search, i], A, PERROR = Err_A, $
                               error = err_img[location[i]-search:location[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
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
    save, brightness, linewidth, filename = outdir + 'brightness.sav'
    beep
endif  ; part==3



; =====================================================================================================================
; Part 4 : put everything together and build images and movies of the exosphere.
; =====================================================================================================================
if part eq 4 or part eq 99 then begin  
    if part ne 99 then restore, outdir + 'brightness.sav'                 ; contains brightness, linewidth
    if part ne 99 then restore, outdir + 'shift_array.sav'                ; contains shift_array, std_devs, aligned_imaging_cube, aligned_spectra_cube, xslit, slit_width
    aligned_raw    = aligned_imaging_cube                                 ; keep the "raw" aligned imaging cube, as we'll fill the slit pixels in aligned_imaging_cube itself
    imaging_cube   = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
    ireference     = readfits(Mercury_dir + ref_file, ref_header, /silent)
    integration    = sxpar(header, 'EXPOSURE')                            ; integration time for individual frames
    s              = size(imaging_cube, /dimensions)                      ; size of frame
    ssc            = size(aligned_spectra_cube, /dimensions)              ; initial size of aligned_spectra_cube
    sb             = size(brightness, /dimensions)                        ; initial size of brightness

    ; reform arrays based on difference in spectral and imaging platescales
    s1_new         = round(s[1] * platescale_fac)
    brightness     = congrid(brightness, s[2], s1_new)
    rotated_brightness = rotate(brightness,7)                             ; rotate brightness for 1-to-1 mapping into movie frames (note 7 = transposed and rotated ccw 270 deg)
    aligned_spectra_cube = congrid(aligned_spectra_cube, ssc[0], s1_new, s[2])
    
    ; Find "home", or the approximate location of the slit (to map extracted Na to, just used for big_frame)
    iflat          = MRDFITS(flat_dir + flat_file, 0, header, /fscale )   ; read flat
    flat           = iflat[ims(0):ims(1),ims(2):ims(3)]                   ; extract imaging portion
    home           = (where(total(flat,2) eq min(total(flat,2))))[0]      ; "home" position in x
    
    ; Fill the slit pixels (for nicer plotting outputs)
    for iframe = 0, ssc[2]-1 do begin
        frame      = reform(aligned_imaging_cube(*,*,iframe))             ; extract individual frame
        frame[ slit_indices[0,*], slit_indices[1,*] ] = !values.F_Nan     ; define slit pixels as NaN
;        frame[xslit-slit_width:xslit+slit_width,*] = !values.F_Nan        ; define slit pixels as NaN
        fill_missing_rand, frame, !values.F_Nan, 1                        ; replace slit pixels based on nearby pixels
        aligned_imaging_cube[*,*,iframe] = frame                          ; update with "fixed" frame 
    endfor
    
    big_cube_spec  = fltarr(s[0], s1_new, s[2]) + !Values.F_NAN             ; initialize cube of spectral image data
    img_cube       = fltarr(s[0], s[1])
    img_count      = fltarr(s[0], s[1])                                   ; array to track number of "real" pixels (ie not obscured by slit)
    spec_cube      = fltarr(s[0], s[1])
    spec_dwell     = fltarr(s[0], s[1])                                   ; keep track of the # of frames in each pixel
    spec_temp      = fltarr(s[0], s[1])                                   ; spectral array for temporary use 
    big_frame_total = fltarr(s[0], s[1]+sb[1])                            ; combined imaging and spectral frame

    ; Find image and spec x/y regions to extract and define shorthand variables for x and y sizes
    ymax_spec      = fix(mean( where(normalize(total(rotated_brightness,1,/nan)) ge 0.05) ))
    ymax_img       = (where(total(total(aligned_imaging_cube,3),1) eq max(total(total(aligned_imaging_cube,3),1))))[0]
    xmax_img       = (where(total(total(aligned_imaging_cube,3),2) eq max(total(total(aligned_imaging_cube,3),2))))[0]
    ytot_img       = normalize(total(total(aligned_imaging_cube,3),1) - mean(total(total(aligned_imaging_cube,3),1)))
    img_width      = n_elements(where(ytot_img ge 0.1))
    img_extraction = [xmax_img-img_width,xmax_img+img_width,ymax_img-img_width,ymax_img+img_width]
    spec_extraction = [xmax_img-img_width,xmax_img+img_width,ymax_spec-img_width,ymax_spec+img_width]
    xs_spec        = spec_extraction(1) - spec_extraction(0) + 1
    ys_spec        = spec_extraction(3) - spec_extraction(2) + 1
    xs_img         = img_extraction(1) - img_extraction(0) + 1
    ys_img         = img_extraction(3) - img_extraction(2) + 1

    ; Extract Mercury reference image for overplotting
    reference      = reform(ireference(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3),0))
    if do_smooth then reference = smooth(reference,smooth_width,/edge_truncate)
   
    ; Set up windows for plotting (we always plot in part 4 since the movies use those plot windows)
    window, 0, xs = s[0], ys = 2*s[1], xpos=winpos_x, ypos=winpos_y
    window, 1, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x+s[0]+20, ypos=winpos_y   
    window, 2, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x+s[0]+20, ypos=winpos_y+ys_img*movie_scale+40
    window, 3, xs = s[0], ys = 2*s[1], xpos=winpos_x, ypos=winpos_y+2*s[1]+40

    ; Some video related initializations
    mpgFilename    = outdir + 'Mercury Na image.mp4'
    mpgFilename2   = outdir + 'Mercury disk.mp4'
    mpgFilename3   = outdir + 'Mercury scan.mp4'
    video          = IDLffVideoWrite(mpgFilename, format='mp4')
    video2         = IDLffVideoWrite(mpgFilename2, format='mp4')
    video3         = IDLffVideoWrite(mpgFilename3, format='mp4')
    framedims      = [xs_img*movie_scale, ys_img*movie_scale]
    framedims2     = [xs_img*movie_scale, ys_img*movie_scale]
    framedims3     = [s[0], 2*s[1]]
    stream         = video.AddVideoStream(framedims[0], framedims[1], framerate)
    stream2        = video2.AddVideoStream(framedims2[0], framedims2[1], framerate)
    stream3        = video3.AddVideoStream(framedims3[0], framedims3[1], framerate)

    good_frames    = where(std_devs ge (median(std_devs)-stddev_cutoff*stddev(std_devs)), complement=bad_frames)  ; these are the good frames we'll use below (in terms of good seeing)
    ngood_shifts   = n_elements(where(shift_array(*,0) ne -666))          ; these are the good frames (in terms of a reasonable shift_array value)

    for ig = 0, n_elements(good_frames)-1 do begin                        ; loop over good frames
    
        i          = good_frames(ig)
        if shift_array[i,0] eq -666 then goto, skipbad                    ; skip over bad frames if they exist (ie don't add to movie)
        
        big_frame  = fltarr(s[0], s[1]+s1_new) + !Values.F_NaN
        empty_frame = make_array(s[0], s1_new, /float, value = !values.F_NaN)
        empty_frame[home+shift_array[i,0]-1:home+shift_array[i,0]+1,*] = [rotated_brightness[i, *], rotated_brightness[i, *], rotated_brightness[i, *]] ; need to use a more realistic representation of the slit width...
        big_frame[*, s[1]:*] = empty_frame
        big_frame[*, 0:s[1]-1] = reform(aligned_raw[*,*,i]) * max(empty_frame,/nan) / max(reform(aligned_raw[*,*,i]),/nan)  ; the max ratios scale the bottom of big_frame based on the relative values in imaging and spectral channels
        big_frame_total = big_frame_total + big_frame

        ; Create frames for "Mercury scan" movie -- this shows the image and extracted Na spectrum frame-by-frame
        wset, 0
        cgloadct, ct, /silent
        cgimage, bytscl(big_frame, 0, big_frame_max)
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_scale
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, .75, .97, 'Frame # ' + strfix(i+1), /data, col='black', /normal, charthick=1.5
        movie_png = cgsnapshot(filename = outdir + 'Mercury scan.png', /nodialog)
        image = read_png(outdir + 'Mercury scan.png')
        void = video3 -> Put(stream3, image)

        ; Create frames for "Mercury Na image" movie -- this shows the Na values along the slit frame-by-frame
        wset, 1
        spec_cut_for_big_cube_spec = reform(big_frame(*,s[1]:*))
        if do_smooth then big_cube_spec[*,*,i] = smooth( spec_cut_for_big_cube_spec, [0,smooth_width], /nan ) $
                     else big_cube_spec[*,*,i] = spec_cut_for_big_cube_spec
        img_cube = img_cube + reform(big_frame(*,0:s[1]))                 ; running total of image cube extraction
        spec_frame = reform(big_frame(*,s[1]:*))                          ; the spectral frame portion
        spec_real = where(finite(spec_frame),complement=spec_bkg)         ; find elements with real values, otherwise assign to background
        spec_summed = total(finite(spec_frame),2)                         ; by summing rows we can find the columns corresponding to the slit location
        slit_columns = where(spec_summed eq max(spec_summed))             ; columns corresponding to slit location
        if n_elements(slit_columns) le 3 then begin                       ; sometimes we have a bad frame, in which case slit_columns thinks everything is beneath the slit
          img_count = img_count + 1.                                      ; increment "real" pixels
          img_count(slit_columns,*) = img_count(slit_columns,*) - 1.      ; but also account for the pixels covered by the slit
          spec_dwell(min(slit_columns):max(slit_columns),*) = spec_dwell(min(slit_columns):max(slit_columns),*) + integration
        endif
        spec_frame(spec_bkg) = 0d
        if do_smooth then spec_frame = smooth(spec_frame, [0,smooth_width])
        spec_cube(spec_real) = spec_cube(spec_real) + spec_frame(spec_real)
        spec_dwell_num = where(spec_dwell gt 0.)                          ; find pixels with spectral information included
        if ig lt n_elements(good_frames)-1 then $
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num) + spec_frame(spec_dwell_num) else $   
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; ie just don't plot slit contribution for last frame
        spec_frame_slit = spec_frame(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
        spec_frame_slit_value = where(spec_frame_slit gt 0., complement=spec_frame_slit_bkg)
        spec_frame_slit_2 = spec_frame_slit
        spec_frame_slit_2(spec_frame_slit_bkg) = 0.;!Values.F_NAN
        blah = big_cube_spec(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3),*)
        if do_median then spec_movie_temp = median(blah, dimension=3, /even) $
                     else spec_movie_temp = mean(blah, dimension=3, /nan)
        print, 'Frame ' + strfix(i) + ', max Na image value = ', max(spec_movie_temp,/nan)
        spec_movie = bytscl(spec_movie_temp, 0., max_spec)
        axis_format = {XTicklen:0.0001, YTickLen:0.0001, Xthick:4, Ythick:4}
        if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                     else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
        if do_rotate then cgimage, rotate(spec_frame_slit_2,5), transparent=20, missing_Value=0.0, ctindex=0 $
                     else cgimage, spec_frame_slit_2,  transparent=20, missing_Value=0.0, ctindex=0
        img_cube_cut = img_cube(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
        cgcontour, normalize(reference), /onimage, levels=[contour_outline], label=0
        cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
        cgtext, .68, .05, '# frames = ' + strfix(ig+1), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.4*movie_Scale
        cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.5*movie_Scale, font=-1
        movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
        image = read_png(outdir + 'Mercury Na.png')
        void = video -> Put(stream, image)

        ; Create frames for "Mercury disk" movie -- this shows buildup of the white light disk image of Mercury
        wset, 2
        if do_smooth then img_count_cut = gauss_smooth(img_count(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3)),3,width=smooth_width,/edge_truncate) $
                     else img_count_cut = img_count(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
        img_movie = congrid(img_cube_cut * (count_factor*(max(img_count_cut)/img_count_cut)), xs_img*movie_scale, ys_img*movie_scale)
        cgloadct, 0, /silent
        cgimage, img_movie
        cgcolorbar, /vertical, charsize=1.e-3
        cgtext, .68, .05, '# frames = ' + strfix(ig+1), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.4*movie_Scale
        cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.5*movie_Scale
        movie_png = cgsnapshot(filename = outdir + 'Mercury disk.png', /nodialog)
        image = read_png(outdir + 'Mercury disk.png')
        void = video2 -> Put(stream2, image)

        wset, 3
        big_frame_total_bottom = big_frame_total * 0. + !Values.F_NaN
        big_frame_total_bottom(*,0:s[1]-1) = big_frame_total(*,0:s[1]-1)
        big_frame_total_top = big_frame_total * 0. + !Values.F_NaN
        big_frame_total_top(*,s[1]:*) = big_frame_total(*,s[1]:*)
        if do_smooth then begin
          cgimage, smooth( bytscl(big_frame_total_top, 0, big_frame_max*float(ig+1)), [2,smooth_width] )
          cgimage, smooth(big_frame_total_bottom, smooth_width, /nan), missing_value=!Values.F_NaN
        endif else begin
          cgimage, bytscl(big_frame_total_top, 0, big_frame_max*float(ig+1))
          cgimage, big_frame_total_bottom, missing_value=!Values.F_NaN
        endelse
        skipbad: wait, 0.02   
    endfor ;ig
   
    ; LOOP over bad frames too for comparison figure ???
    if do_bad then begin
       img_cube_bad = fltarr(s[0], s[1])
       for ib = 0, n_elements(bad_frames)-1 do begin                         ; loop over bad frames
           i = bad_frames(ib)
           if shift_array[i,0] eq -666 then goto, skipbad2
           big_frame  = fltarr(s[0], s[1]+s1_new)
           empty_frame = make_array(s[0], s1_new, /float, value = !values.F_NaN)
           empty_frame[home+shift_array[i,0]-2:home+shift_array[i,0]+2,*] = [brightness[i,*], brightness[i,*], brightness[i,*], brightness[i,*], brightness[i,*]] ;Make it bigger?
           big_frame[*, s[1]:*] = rotate(empty_frame, 7) ;* 5. ;the factor of 5 is a HACK FOR NOW
           big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i < ngood_shifts-1])
           img_cube_bad = img_cube_bad + reform(big_frame(*,0:s[1]))
           skipbad2: wait, 0.02
        endfor ;ib     
        img_cube_cut_bad = img_cube_bad(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
    endif ;do_bad==1
   
    spec_dwell_num = where(spec_dwell gt 0.)                               ; find pixels with spectral information included
    spec_cube(spec_dwell_num) = spec_cube(spec_dwell_num) / ( spec_dwell(spec_dwell_num) > 1.) ; normalize by # of frames summed, don't divide by zero
   
    ; Now add some blank frames to the Na image movie before plotting a smoothed frame
    wset, 1
    spec_raw = spec_cube(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
    spec_raw = congrid(spec_raw, xs_spec*movie_scale, ys_spec*movie_scale) ; embiggen
    loadct, ct, /silent
    if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                 else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
    cgcontour, normalize(reference), /onimage, levels=[contour_outline], label=0
    cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
    cgtext, .68, .05, '# frames = ' + strfix(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.4*movie_Scale
    cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.5*movie_Scale, font=-1
    movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
    image = read_png(outdir + 'Mercury Na.png')

    for i = 0, 30 do void = video -> Put(stream, image)
    if do_rotate then cgimage, rotate(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate), 7) $
                 else cgimage, smooth(spec_raw,smooth_width*movie_scale,/edge_truncate)
    cgcontour, normalize(reference), /onimage, levels=[contour_outline], label=0
    cgtext, .68, .05, '# frames = ' + strfix(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.4*movie_Scale
    cgtext, 0.02, 0.96, 'Mercury Na (smoothed=' + strfix(smooth_width) + ')', /normal, charthick=0.75*movie_Scale, charsize=0.5*movie_Scale, font=-1
    movie_png = cgsnapshot(filename = outdir + 'Mercury Na - smoothed.png', /nodialog)
    image = read_png(outdir + 'Mercury Na - smoothed.png')
    for i = 0, 30 do void = video -> Put(stream, image)                    ; add a bunch of identical frames so we can get a sense of how it looks before smoothing
      
    video -> Cleanup
    video2 -> Cleanup
    video3 -> Cleanup

    ; Overplot Na contours on top of Mercury disk
    loadct, 0, /silent
    cgimage, bytscl(img_movie), /axes
    spec_smoothed = bytscl(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate))
    cgcontour, normalize(double(spec_smoothed)), /onimage, color='gold', levels=[0.2,0.4,0.6,0.8]
    al_legend, ['Disk','Na'], box=0, /bottom, /right, charsize=0.5*movie_scale, line=[0,0], linsize=0.2, thick=[6,2], color=['gray','gold'], textcolor='white'
    combined = cgsnapshot(filename = outdir + 'Mercury Na + disk.png', /nodialog)

   ; Generate the "bad images" Mercury disk for comparison
   if do_bad then begin
        bad_img = bytscl(img_cube_cut_bad, 0.2*max(img_cube_cut_bad), max(img_cube_cut_bad))
        cgimage, bad_img
        cgcolorbar, /vertical, charsize=1.e-3
        cgtext, .65, .05, '# frames = ' + strfix(ib), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_Scale
        cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
        bad_png = cgsnapshot(filename = outdir + 'Mercury disk - bad frames.png', /nodialog)

        ; And might as well compare the good vs bad contours while we're here...
        cgimage, img_movie, /axes
        cgcontour, normalize(img_movie), /onimage, levels=[0.2,0.4,0.6,0.8], color='gold'
        cgcontour, normalize(img_cube_cut_bad), /onimage, levels=[0.2,0.4,0.6,0.8], color='red'
        al_legend, ['good','bad'], box=0, /bottom, /right, charsize=2, line=[0,0], linsize=0.2, thick=[4,4], color=['orange','red'], textcolor=['orange','red']
        compare = cgsnapshot(filename = outdir + 'Mercury disk - good vs bad contours.png', /nodialog)
    endif ;do_bad==1
   
    ; Create image showing the spatial coverage of the slit position
    wset, 2
    dwell_img = spec_dwell(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
    loadct, 22, /silent
    cgimage, bytscl(dwell_img, 0., max(dwell_img)), /axes, title='Slit Coverage'
    cgcontour, normalize(reference), /onimage, levels=[contour_outline], label=0
    cgloadct, 22, ncolors=max(spec_dwell)/integration, bottom=0, /silent  
    cgcolorbar, /vertical, maxrange=max(spec_dwell), title='Total Na data coverage (sec)', /discrete, ncolors=max(spec_dwell)/integration, position=[0.9,0.2,0.95,0.8], tickinterval=integration*4.
    combined = cgsnapshot(filename = outdir + 'Mercury slit coverage.png', /nodialog)

    ; Write aligned images as a new reference
    time = strsplit(systime(/UTC),' ',/extract)
    MWRFITS, total(aligned_imaging_cube,3), Mercury_dir + 'Mercury_ref_image_' + time(2) + time(1) + time(4) + '.fits', ref_header, /create
      
    beep
    stop
endif ;part==4

  
end



