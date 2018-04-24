
;+
; NAME:
;   PHI
;
; PURPOSE:
;   A Bayesian MCMC method for 2D photometric decompositions of
;   galaxies and other astronomical objects.    
;
;
;
;
;
;
;
;-
function psf_gaussian, parameters, NPIXEL=npixel, NDIMENSION=ndim, FWHM=fwhm,  $
                        DOUBLE = double, CENTROID=cntrd, ST_DEV=st_dev,  $
                        XY_CORREL=xy_corr, NORMALIZE=normalize
;+
; NAME:
;       PSF_GAUSSIAN
;
; PURPOSE:
;       Create a 1-d, 2-d, or 3-d Gaussian with specified FWHM, center 
; EXPLANATION:
;       Return a point spread function having Gaussian profiles,
;       as either a 1D vector, a 2D image, or 3D volumetric-data.
;
; CALLING SEQUENCE:
;       psf = psf_Gaussian( NPIXEL=, FWHM= , CENTROID = 
;                     [ /DOUBLE, /NORMALIZE, ST_DEV=,  NDIMEN= ] ) 
; or:
;       psf = psf_Gaussian( parameters, NPIXEL = ,NDIMEN = )
;
; REQUIRED INPUT KEYWORD:
;       NPIXEL = number pixels for each dimension, specify as an array,
;               or just one number to make all sizes equal.
;

;
; INPUTS (optional):
;
;       parameters = an NDIMEN by 3 array giving for each dimension:
;                       [ maxval, center, st_dev ],  overrides other keywords.
;
; EXAMPLE:
;       (1) Create a 31 x 31 array containing a normalized centered Gaussian 
;       with an X FWHM = 4.3 and a Y FWHM = 3.6
;
;       IDL> array = PSF_GAUSSIAN( Npixel=31, FWHM=[4.3,3.6], /NORMAL )
;
;       (2) Create a 50 pixel 1-d Gaussian vector with a maximum of 12, 
;          centered at  pixel 23 with a sigma of 19.2
;
;       IDL> psf = psf_gaussian([12,23,19.2],npixel=50)

        On_error,2
	compile_opt idl2

        if (N_params() LT 1 ) and $
            ~(keyword_set( FWHM) || keyword_set(ST_DEV)) then begin
                print,'Syntax - psf = PSF_GAUSSIAN( parameters, NPIXEL = )'
                print, $
       'or       psf = PSF_GAUSSIAN( FWHM = ,ST_DEV = ,NPIXEL = ,[CENTROID = ])'
                return, -1
        endif

        sp = size( parameters )
        if sp[0] EQ 1 then begin               ;Vector supplied?
                ndim = 1
                factor = parameters[0]
                cntrd = parameters[1]
                st_dev = parameters[2] 
         endif  else  if (sp[0] GE 1) then begin    ;Ndimen x 3 array supplied?
                 ndim = sp[1]
                 factor = total( parameters[*,0] )/float( ndim )
                cntrd = parameters[*,1]
                st_dev = parameters[*,2]
           endif

        double = keyword_set(double)
        if double then idltype = 5 else idltype = 4
        if N_elements( ndim ) NE 1 then ndim=2
        ndim = ndim>1

        if N_elements( npixel ) LE 0 then begin
                message,"must specify size of result with NPIX=",/INFO
                return,(-1)
          endif else begin 
	      npix = npixel
	      if N_elements( npix ) LT ndim then npix = replicate( npix[0], ndim )
         endelse

        if (N_elements( cntrd ) LT ndim) && (N_elements( cntrd ) GT 0) then $
                        cntrd = replicate( cntrd[0], ndim )

        if N_elements( cntrd ) LE 0 then cntrd=(npix-1)/2. 
        if N_elements( fwhm ) GT 0 then begin 
               st_dev = fwhm/( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
               if ~double then st_dev  = float(st_dev)
        endif 

        if N_elements( st_dev ) LE 0 then begin
                message,"must specify ST_DEV= or FWHM=",/INFO
                return,(-1)
          endif

        if N_elements( st_dev ) LT ndim then $
                        st_dev = replicate( st_dev[0], ndim )

        CASE ndim OF

        1: BEGIN
                x = findgen( npix[0] ) - cntrd[0]
                psf = gaussian( x, [1,0,st_dev] )
             END

        2: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]

                if N_elements( xy_corr ) EQ 1 then begin
                        sigfac = 1 / (2. * st_dev^2 )
                        y2 = sigfac[1] * y^2
                        x1 = sigfac[0] * x
                        yc = y * ( xy_corr/(st_dev[0]*st_dev[1]) )
                        for j=0,npix[1]-1 do begin
                                zz = x * (yc[j] + x1) + y2[j]
                                w = where( zz LT 86, nw )
                                if (nw GT 0) then psf[w,j] = exp( -zz[w] )
                          endfor
                  endif else begin
                        psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE=double )
                        psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE=double )
                        error = check_math(/print, MASK=32)
                        save_except = !EXCEPT & !EXCEPT = 0
                        for j=0,npix[1]-1 do psf[0,j] = psfx * psfy[j]
                        error = check_math(MASK=32)    ;Clear floating underflow
                        !EXCEPT = save_except  
                   endelse
             END

        3: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]
                z = make_array( npix[2], /INDEX, TYPE=idltype ) - cntrd[2]
                psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE = double )
                psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE = double)
                psfz = gaussian( z, [ 1, 0, st_dev[2] ], DOUBLE = double )
                error = check_math(MASK=32,/PRINT)
                save_except = !EXCEPT & !EXCEPT = 0
                for k=0,npix[2]-1 do begin
                    for j=0,npix[1]-1 do psf[0,j,k] = psfx * psfy[j] * psfz[k]
                 endfor
                 error = check_math(MASK=32)
                 !EXCEPT = save_except  
             END

        ENDCASE

        if keyword_set( normalize ) then return, psf/total( psf )

        if N_elements( factor ) EQ 1 then begin
                if (factor NE 1) then return,factor*psf else return,psf
           endif else return, psf
     end

function psf_construct, ndimension=imag_dim, func=func, fwhm=fwhm , betaa=betaa
  
  xc=imag_dim/2. 
  yc=imag_dim/2.
  nx=imag_dim
  ny=imag_dim
  imag_tot=nx*ny
  i=lindgen(imag_tot)
  rx=(double(i mod nx)-xc)
  ry=(double(i/nx)-yc)
  r=sqrt(rx^2D + ry^2D)
  
  ;Gaussian PSF
  if (func eq 0) then begin
     fwhm_gauss=fwhm
     sigma=fwhm/(2D *sqrt(2D *alog(2D)))
     gauss=(1D /(2D *sigma^2D *!DPI))*exp(-0.5D *(r/sigma)^2D)
     gauss=reform(gauss,nx,ny)
     return,gauss
  endif 
  ;Moffat PSF
  if (func eq 1) then begin
     alpha=fwhm/(2D *sqrt(2D ^(1D /betaa) -1D))
     moffat=((betaa-1D)/(!DPI*alpha^2D))*(1D +(r/alpha)^2D)^(-1D*betaa)
     moffat=reform(moffat,nx,ny)
     return,moffat
  endif
  ;Double-gaussian PSF
  if (func eq 2) then begin
     gauss1 =  psf_gaussian(NPIXEL = imag_dim,FWHM=[FWHM[0],FWHM[0]],/normal)
     gauss2 =  psf_gaussian(NPIXEL = imag_dim,FWHM=[FWHM[1],FWHM[1]],/normal)
     gauss = gauss1 + gauss2 
     gauss = gauss/total(gauss)
     
     return,gauss
  endif 
end 

;PERCENTILES
function percentiles,data,value=value
 
  result = -1
  n = n_elements(data)
  if (n le 0) then return,result ; error : data not defined
 
; check if speficic percentiles requested - if not: set standard
  if(not keyword_set(value)) then value = [ 0., 0.25, 0.5, 0.75, 1.0 ]
 
  ix = sort(data)
 
; loop through percentile values, get indices and add to result
; This is all we need since computing percentiles is nothing more
; than counting in a sorted array.
  for i=0,n_elements(value)-1 do begin
 
     if(value(i) lt 0. OR value(i) gt 1.) then return,-1
     if(value(i) le 0.5) then ind = long(value(i)*n)    $
     else ind = long(value(i)*(n+1))
     if (ind ge n) then ind = n-1 ; small fix for small n
                                   ; (or value eq 1.)

     if(i eq 0) then result = data(ix(ind))  $
     else result = [result, data(ix(ind)) ]
  endfor
 
  return,result
end
FUNCTION postinfo,data

  num = n_elements(data)
;FIND medium 
  sdata = sort(data)
  data1 = data[sdata]
  meddata = median(data1)
  
  
;Find 1 sigma regime 
  indlow = round(0.16*num)
  indhi = num - indlow

  datalow = data1[indlow]
  datahi = data1[indhi]
  
;Return the medium, and lower and upper 1sigma limits 
  struc = {med:meddata,siglow:datalow,sighigh:datahi} 
  return, struc
END

;AUTO-CORRELATION
Function autocorr, data, ncols = ncols, nlags = nlags
  
  ;Removes unwanted cols 
  data = REFORM(data)

  ndims = N_elements(data)
  
  
  if (n_elements(nlags) eq 0) then nlags = 50D
  if ((ndims - 2l) ge nlags) then lag = range(0,nlags,nlags+1) else lag = range(0,ndims,ndims+1)
  
  ACF = A_CORRELATE(data, lag)
 
  
  return, ACF
END

;PRIORS AND STEP FUNCTIONS
Function Uprior, N, A, B, seed
  x = fltarr(n_elements(A),N)
  for i=0, n_elements(A)-1 do begin
     x[i,*] = A[i] + (B[i]-A[i])*RANDOMU(seed, N,/DOUBLE)
  endfor
  return, x
END

FUNCTION Nrandom, AVG=AVG, SIG=SIG
                                ;Accepts a grid of sigmas or average
                                ;values and returns a grid with the
                                ;same dimensions containing normal
                                ;random numbers
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(SIG) then SIG = 1D
  
  avgnum = n_elements(AVG)
  signum = n_elements(SIG)
  
  
  sigsize = size(SIG,/DIMENSIONS)
  if (sigsize[0] eq 0) then sigsize = 1
  avgsize = size(AVG,/DIMENSIONS)
  if (sigsize[0] eq 0) then avgsize = 1

  R = make_array(signum,/DOUBLE)
  for i=0, signum-1 do begin
     R[i] = RANDOMN(seed)
  endfor
  
  if (n_elements(sigsize) gt 1) then R = reform(R,sigsize[0],sigsize[1])
  Result = (R * SIG) + AVG

  return, Result

END
Function Nprior, N, AVG, SIG, seed
  R = fltarr(N)
  rx = fltarr(N)
  ry = fltarr(N)
  
  for i=0, N-1 do begin 
     repeat begin 
        x = 2*RANDOMU(seed, 1)-1
        y = 2*RANDOMU(seed, 1)-1
        R_s = x^2 + y^2
     endrep until R_s LT 1.0 
     R(i) = sqrt(R_s)          
     rx(i) = x
     ry(i) = y
  endfor
  
  G1 = sqrt(-ALOG(R))*2*rx/R
  ;G2 = sqrt(-ALOG(R))*2*ry/R
  G1 = AVG+(G1*SIG)
  ;G2 = AVG+(G2*SIG)
  
  return, G1
END

Function Gprior, x, AVG=AVG, SIG=SIG, seed=seed
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(SIG) then SIG = 1D
  
  Px = (1D/(sqrt(2D*!DPI)*SIG)) * EXP(-0.5*(((x - AVG)/sig)^2D))
  
  return,Px
END

Function Uniprior, x, A=A, B=B, seed=seed
   if NOT keyword_set(A) then A = 0D
   if NOT keyword_set(B) then B = 1D
  
   n = n_elements(x)
   Px = fltarr(n)
   for i=0, n-1 do begin
      if (x[i] le B) and (x[i] ge A) then begin
         Px[i] = 1D / (B - A)
      endif else begin
         Px[i] = 0D
      endelse
   endfor
   
   return, Px
END

Function vprior,index,x,a=a
  case index of
     ;0: x = Uprior(1, (x - (x/10D)),  (x + (x/10D)))
     0: begin
        x1 = x - a
        x2 = x + a
        xnew = Uprior(1, x1, x2)
     end
     1:begin 
        xnew = Nprior(1, x, a)
     end
  ENDCASE
  return, xnew
END


FUNCTION priors, pars
  common setup1, tmcmc, bagal_like, nchain, nstep, tstep, pridis, f1dfit
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  common priorss,  pri_pars
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  sorder2 = where(nlist eq 'ser2X0', sernn2)
  sorder3 = where(nlist eq 'ser3X0', sernn3)
  sorder4 = where(nlist eq 'ser4X0', sernn4)
  beorder = where(nlist eq 'bexpX0', bexpnn)
  forder = where(nlist eq 'ferX0', fernn)
  forder2 = where(nlist eq 'fer2X0', fernn2)
  
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)


  Priserx0_temp = !Values.F_NAN
  Prisery0_temp = !Values.F_NAN
  Priser2x0_temp = !Values.F_NAN
  Priser2y0_temp = !Values.F_NAN
  Priser3x0_temp = !Values.F_NAN
  Priser3y0_temp = !Values.F_NAN
  Priser4x0_temp = !Values.F_NAN
  Priser4y0_temp = !Values.F_NAN
     
  Priexpx0_temp = !Values.F_NAN
  Priexpy0_temp = !Values.F_NAN

  Pribexpx0_temp = !Values.F_NAN
  Pribexpy0_temp = !Values.F_NAN

  Priferx0_temp = !Values.F_NAN
  Prifer2x0_temp = !Values.F_NAN
  
  if (f1dfit eq 0) then begin 
     skip = 1+i
     readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
              PriserX0, PriserY0, PriIe, PriRe, Prin, PriserE, PriserPA,$
              PriexpX0, PriexpY0, PriI0, Prih, PriexpE, PriexpPA,$
              Priser2X0, Priser2Y0, PriIe2, PriRe2, Prin2, Priser2E, Priser2PA,$
              Priser3X0, Priser3Y0, PriIe3, PriRe3, Prin3, Priser3E, Priser3PA,$
              Priser4X0, Priser4Y0, PriIe4, PriRe4, Prin4, Priser4E, Priser4PA,$
              FORMAT='A,F,F,F,'+$
              'F,F,F,F,F,F,F,'+$
              'F,F,F,F,F,F,'+$
              'F,F,F,F,F,F,F,'+$
              'F,F,F,F,F,F,F,'+$
              'F,F,F,F,F,F,F',$
              SKIPLINE=skip,/SILENT,NUMLINE=1

     if (sernn eq 1) then PriserPA = cos(PriserPA*!DPI/180D)
     if (sernn2 eq 1) then Priser2PA = cos(Priser2PA*!DPI/180D)
     if (sernn3 eq 1) then Priser3PA = cos(Priser3PA*!DPI/180D)
     if (sernn4 eq 1) then Priser4PA = cos(Priser4PA*!DPI/180D)
     if (expnn eq 1) then PriexpPA = cos(PriexpPA*!DPI/180D)
     
     
  endif else begin
     
    
     if (sernn eq 1) then begin
        PriserX0 = pri_pars[sorder]
        PriserY0 = pri_pars[sorder+1]
        
        PriIe= 10D^(pri_pars[sorder+2])
        PriRe = 10D^(pri_pars[sorder+3])
        Prin = pri_pars[sorder+4]
        PriserE = pri_pars[sorder+5]
        PriserPA = (pri_pars[sorder+6])
       
     endif
      if (sernn2 eq 1) then begin
        Priser2X0 = pri_pars[sorder2]
        Priser2Y0 = pri_pars[sorder2+1]
        
        PriIe2= 10D^(pri_pars[sorder2+2])
        PriRe2 = 10D^(pri_pars[sorder2+3])
        Prin2 = pri_pars[sorder2+4]
        Priser2E = pri_pars[sorder2+5]
        Priser2PA = (pri_pars[sorder2+6])
       
     endif
       if (sernn3 eq 1) then begin
        Priser3X0 = pri_pars[sorder3]
        Priser3Y0 = pri_pars[sorder3+1]
        
        PriIe3= 10D^(pri_pars[sorder3+2])
        PriRe3 = 10D^(pri_pars[sorder3+3])
        Prin3 = pri_pars[sorder3+4]
        Priser3E = pri_pars[sorder3+5]
        Priser3PA = (pri_pars[sorder3+6])
       
     endif
        if (sernn4 eq 1) then begin
        Priser4X0 = pri_pars[sorder4]
        Priser4Y0 = pri_pars[sorder4+1]
        
        PriIe4= 10D^(pri_pars[sorder4+2])
        PriRe4 = 10D^(pri_pars[sorder4+3])
        Prin4 = pri_pars[sorder4+4]
        Priser4E = pri_pars[sorder4+5]
        Priser4PA = (pri_pars[sorder4+6])
       
     endif

    
     if (expnn eq 1) then begin
        PriexpX0 = pri_pars[eorder]
        PriexpY0 = pri_pars[eorder+1]
        
        PriI0 = 10D^(pri_pars[eorder+2])
        Prih = 10D^(pri_pars[eorder+3])
        PriexpE = pri_pars[eorder+4]
        PriexpPA = (pri_pars[eorder+5])
     endif

  endelse
  
  
  

 
  
  prior = MAKE_ARRAY(1,/STRING)

  
  ;SERSIC PROFILE
  if sernn gt 0 then begin
     
     
     serX0 = pars[sorder]
     serY0 = pars[sorder+1]
     fserx0 = flist[0,sorder]
     fsery0 = flist[0,sorder+1]
     pserx0 = flist[1,sorder]
     psery0 = flist[1,sorder+1]
     Priserx0_temp = Priserx0
     Prisery0_temp = Prisery0

     Ie = 10D^(pars[sorder+2])
     fIe = flist[0,sorder+2]
     pIe = flist[1,sorder+2]

     Re = 10D^(pars[sorder+3])
     fRe = flist[0,sorder+3]
     pRe = flist[1,sorder+3]
     
     n = (pars[sorder+4])
     fn = flist[0,sorder+4]
     pn = flist[1,sorder+4]
     
     serE = pars[sorder+5]
     fserE = flist[0,sorder+5]
     pserE = flist[1,sorder+5]
     
     serPA = (pars[sorder+6])
     fserPA = flist[0,sorder+6]
     pserPA = flist[1,sorder+6]
     
     if fserX0 eq 1 or fserY0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fserX0 eq 1 then begin 
           CASE pserx0 of 
              0: begin
                 AserX0 = (PriserX0)-3
                 BserX0 = (PriserX0/2D)+3
                 if serX0 ge Aserx0 and serX0 le BserX0 then begin
                    Prserx0 = 1D ;Uniprior(serX0, A=Aserx0, B=BserX0)
                 endif else PrserX0 = 0.
              end
              1: Prserx0 = Gprior(serX0, AVG=PriserX0, SIG=(PriserX0/3D))
           ENDCASE
           prior = [prior, 'Prserx0']
        endif 
        if fserY0 eq 1 then begin 
           CASE pserx0 of 
              0: begin
                 AserY0 = (PriserY0/2D)-3
                 BserY0 = (PriserY0/2D)+3
                 if serY0 ge AserY0 and serY0 le BserY0 then begin
                    PrserY0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else PrserY0 = 0.
              end
              1: PrserY0 = Gprior(serY0, AVG=PriserY0, SIG=(s[1]/3D))
           ENDCASE
        endif 
        prior = [prior, 'Prsery0']
     endif
    

     if fIe eq 1 then begin
        PriIe = PriIe
        CASE pIe of 
           0: begin
              AIe = 15D ;Primue + 1D ;10D^(-1D * 0.4 * ((Primue +0.2D) - zcalib - 2.5*alog10((scale^2D))))
              BIe = 25D ; Primue - 1D ;10D^(-1D * 0.4 * ((Primue -0.2D) - zcalib - 2.5*alog10((scale^2D))))
              if Ie ge AIe and Ie le BIe then begin
                 PrIe = -1.*ALOG(BIe - AIe);Uniprior(Ie, A=AIe, B=BIe)
              endif else PrIe = 0.
           end
           1: Begin
              Ieval = Ie
              Ievar = 10D
              PrIe = -0.5*( ( ((Ieval-PriIe)^2D)/Ievar) + ALOG(2D*!DPI*Ievar) )
             
           end
           2: Begin
              AIe = 15D
              BIe = 25D
             
              if (Ie gt AIe) and (Ie lt BIe) then Ieval = 15D else Ieval = Ie
              
              Ievar = 25D
              PrIe = -0.5*( ( ((Ieval-20D)^2D)/Ievar) + ALOG(2D*!DPI*Ievar) )
              
           end
        ENDCASE
        prior = [prior, 'PrIe']
     endif 
     
     if fRe eq 1 then begin 
        PriRe = PriRe
        CASE pRe of 
           0: BEGIN
              ARe = 0D ;PriRe/2D
              ;if ARe lt 1. then ARe = 1.
              BRe = radi ;(PriRe + ARe)
              if Re ge ARe and Re le BRe then begin
                 PrRe = -1.*ALOG(BRe - ARe);;1D ;Uniprior(Re, A=ARe, B=BRe)
              endif else PrRe = 0.
           end
           1: begin

              reval = Re
              Revar = 10D;radi/2D
              PrRe = -0.5*( ( ((Reval-PriRe)^2D)/Revar) + ALOG(2D*!DPI*Revar) )
           end
           2: Begin
              ARe = 0
              BRe = radi
              
              priRe = radi/2D
              if (Re lt BRe) then Reval = radi else Reval = Re
              
              Revar = radi/2D
              PrRe = -0.5*( ( ((Reval-PriRe)^2D)/Revar) + ALOG(2D*!DPI*Revar) )
              
           end

        ENDCASE
        prior = [prior,'PrRe']
     endif 
     
     if fn eq 1 then begin 
        CASE pn of 
           0: begin
              An = 0.2
              Bn = 12.
              if n ge An and n le Bn then begin
                 Prn = -1.*ALOG(Bn - An);;1D ;Uniprior(n, A=An, B=Bn)
              endif else Prn = 0.
           end
           1: begin
              nval = n
              nvar = 6.
              Prn = -0.5*( ( ((nval-Prin)^2D)/nvar) + ALOG(2D*!DPI*nvar) )

           end
           2: Begin
              An = 0D
              Bn = 12D
              
              prin = 6D
              if (n gt An) and (n lt Bn) then nval = 12D else nval = n
              
              nvar = 6.
              Prn = -0.5*( ( ((nval-Prin)^2D)/nvar) + ALOG(2D*!DPI*nvar) )
              
           end

        ENDCASE
        prior = [prior, 'Prn']
     endif 
     
     if fserE eq 1 then begin 
        CASE pserE of 
           0: begin
              AserE = 0D
              BserE = 1D
              if serE ge AserE and serE le BserE then begin
                 PrserE = -1.*ALOG(BserE - AserE);;1D ;Uniprior(serE, A=AserE, B=BserE)
              endif else PrserE = 0.
           end
           1: begin
              serEvar = 2D;mean((serE - PriserE)^2D)
              PrserE = -0.5*( ( ((serE-PriserE)^2D)/serEvar) + ALOG(2D*!DPI*serEvar) )
              if serE gt 1 or serE lt 0. then PrserE = 0.
           end
           2: Begin
              PrserE = 0
           end

        ENDCASE
        prior = [prior,'PrserE']
     endif 
     
     if fserPA eq 1 then begin 
        CASE pserPA of 
           0: begin
              AserPA = PriserPA - (priserPA/2D)
              BserPA = PriserPA + (priserPA/2D)
              if serPA ge AserPa and serPa le BserPa then begin 
                 PrserPA = -1.*ALOG(BserPA - AserPA);;1D ;Uniprior(serPA, A=AserPA, B=BserPA)
              endif else PrserPa = 0.
           end
           1: begin
              serPAvar = 1D    ;mean((serPA - PriserPA)^2D)
              PrserPA = -0.5*( ( ((serPA-PriserPA)^2D)/serPAvar) + ALOG(2D*!DPI*serPAvar) )
           end
           2: Begin
              PrserPA = 0
           end
           ENDCASE
        prior = [prior,'PrserPA']
     endif 
     
     
  endif

  ;SERSIC 2 PROFILE
  if sernn2 gt 0 then begin
     
     
     ser2X0 = pars[sorder2]
     ser2Y0 = pars[sorder2+1]
     fser2x0 = flist[0,sorder2]
     fser2y0 = flist[0,sorder2+1]
     pser2x0 = flist[1,sorder2]
     pser2y0 = flist[1,sorder2+1]
     Priser2x0_temp = Priser2x0
     Priser2y0_temp = Priser2y0

     Ie2 = 10D^(pars[sorder2+2])
     fIe2 = flist[0,sorder2+2]
     pIe2 = flist[1,sorder2+2]

     Re2 = 10D^(pars[sorder2+3])
     fRe2 = flist[0,sorder2+3]
     pRe2 = flist[1,sorder2+3]
     
     n2 = (pars[sorder2+4])
     fn2 = flist[0,sorder2+4]
     pn2 = flist[1,sorder2+4]
     
     ser2E = pars[sorder2+5]
     fser2E = flist[0,sorder2+5]
     pser2E = flist[1,sorder2+5]
     
     ser2PA = (pars[sorder2+6])
     fser2PA = flist[0,sorder2+6]
     pser2PA = flist[1,sorder2+6]
     
     if fser2X0 eq 1 or fser2Y0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fser2X0 eq 1 then begin 
           CASE pser2x0 of 
              0: begin
                 Aser2X0 = (Priser2X0)-3
                 Bser2X0 = (Priser2X0/2D)+3
                 if ser2X0 ge Aser2x0 and ser2X0 le Bser2X0 then begin
                    Prser2x0 = 1D ;Uniprior(serX0, A=Aserx0, B=BserX0)
                 endif else Prser2X0 = 0.
              end
              1: Prser2x0 = Gprior(ser2X0, AVG=Priser2X0, SIG=(Priser2X0/3D))
           ENDCASE
           prior = [prior, 'Prser2x0']
        endif 
        if fser2Y0 eq 1 then begin 
           CASE pser2x0 of 
              0: begin
                 Aser2Y0 = (Priser2Y0/2D)-3
                 Bser2Y0 = (Priser2Y0/2D)+3
                 if ser2Y0 ge Aser2Y0 and ser2Y0 le Bser2Y0 then begin
                    Prser2Y0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else Prser2Y0 = 0.
              end
              1: Prser2Y0 = Gprior(ser2Y0, AVG=Priser2Y0, SIG=(s[1]/3D))
           ENDCASE
        endif 
        prior = [prior, 'Prser2y0']
     endif
    

     if fIe2 eq 1 then begin
        PriIe2 = PriIe2
        CASE pIe2 of 
           0: begin
              AIe2 = 15D ;Primue + 1D ;10D^(-1D * 0.4 * ((Primue +0.2D) - zcalib - 2.5*alog10((scale^2D))))
              BIe2 = 25D ; Primue - 1D ;10D^(-1D * 0.4 * ((Primue -0.2D) - zcalib - 2.5*alog10((scale^2D))))
              if Ie2 ge AIe2 and Ie2 le BIe2 then begin
                 PrIe2 = -1.*ALOG(BIe2 - AIe2);Uniprior(Ie, A=AIe, B=BIe)
              endif else PrIe2 = 0.
           end
           1: Begin
              Ieval2 = Ie2
              Ievar2 = 10D
              PrIe2 = -0.5*( ( ((Ieval2-PriIe2)^2D)/Ievar2) + ALOG(2D*!DPI*Ievar2) )
             
           end
           2: Begin
              AIe2 = 15D
              BIe2 = 25D
             
              if (Ie2 gt AIe2) and (Ie2 lt BIe2) then Ieval2 = 15D else Ieval2 = Ie
              
              Ievar2 = 25D
              PrIe2 = -0.5*( ( ((Ieval2-20D)^2D)/Ievar2) + ALOG(2D*!DPI*Ievar2) )
              
           end
        ENDCASE
        prior = [prior, 'PrIe2']
     endif 
     
     if fRe2 eq 1 then begin 
        PriRe2 = PriRe2
        CASE pRe2 of 
           0: BEGIN
              ARe2 = 0D ;PriRe/2D
              ;if ARe lt 1. then ARe = 1.
              BRe2 = radi ;(PriRe + ARe)
              if Re2 ge ARe2 and Re2 le BRe2 then begin
                 PrRe2 = -1.*ALOG(BRe2 - ARe2);;1D ;Uniprior(Re, A=ARe, B=BRe)
              endif else PrRe2 = 0.
           end
           1: begin

              reval2 = Re2
              Revar2 = 10D;radi/2D
              PrRe2 = -0.5*( ( ((Reval2-PriRe2)^2D)/Revar2) + ALOG(2D*!DPI*Revar2) )
           end
           2: Begin
              ARe2 = 0
              BRe2 = radi
              
              priRe2 = radi/2D
              if (Re2 lt BRe2) then Reval2 = radi else Reval2 = Re2
              
              Revar2 = radi/2D
              PrRe2 = -0.5*( ( ((Reval2-PriRe2)^2D)/Revar2) + ALOG(2D*!DPI*Revar2) )
              
           end

        ENDCASE
        prior = [prior,'PrRe2']
     endif 
     
     if fn2 eq 1 then begin 
        CASE pn2 of 
           0: begin
              An2 = 0.2
              Bn2 = 12.
              if n2 ge An2 and n2 le Bn2 then begin
                 Prn2 = -1.*ALOG(Bn2 - An2);;1D ;Uniprior(n, A=An, B=Bn)
              endif else Prn2 = 0.
           end
           1: begin
              nval2 = n2
              nvar2 = 6.
              Prn2 = -0.5*( ( ((nval2-Prin2)^2D)/nvar2) + ALOG(2D*!DPI*nvar2) )

           end
           2: Begin
              An2 = 0D
              Bn2 = 12D
              
              prin2 = 6D
              if (n2 gt An2) and (n2 lt Bn2) then nval2 = 12D else nval2 = n2
              
              nvar2 = 6.
              Prn2 = -0.5*( ( ((nval2-Prin2)^2D)/nvar2) + ALOG(2D*!DPI*nvar2) )
              
           end

        ENDCASE
        prior = [prior, 'Prn2']
     endif 
     
     if fser2E eq 1 then begin 
        CASE pser2E of 
           0: begin
              Aser2E = 0D
              Bser2E = 1D
              if ser2E ge Aser2E and ser2E le Bser2E then begin
                 Prser2E = -1.*ALOG(Bser2E - Aser2E);;1D ;Uniprior(serE, A=AserE, B=BserE)
              endif else Prser2E = 0.
           end
           1: begin
              ser2Evar = 2D;mean((serE - PriserE)^2D)
              Prser2E = -0.5*( ( ((ser2E-Priser2E)^2D)/ser2Evar) + ALOG(2D*!DPI*ser2Evar) )
              if ser2E gt 1 or ser2E lt 0. then Prser2E = 0.
           end
           2: Begin
              Prser2E = 0
           end

        ENDCASE
        prior = [prior,'Prser2E']
     endif 
     
     if fser2PA eq 1 then begin 
        CASE pser2PA of 
           0: begin
              Aser2PA = Priser2PA - (priser2PA/2D)
              Bser2PA = Priser2PA + (priser2PA/2D)
              if ser2PA ge Aser2Pa and ser2Pa le Bser2Pa then begin 
                 Prser2PA = -1.*ALOG(Bser2PA - Aser2PA);;1D ;Uniprior(serPA, A=AserPA, B=BserPA)
              endif else Prser2Pa = 0.
           end
           1: begin
              ser2PAvar = 1D    ;mean((serPA - PriserPA)^2D)
              Prser2PA = -0.5*( ( ((ser2PA-Priser2PA)^2D)/ser2PAvar) + ALOG(2D*!DPI*ser2PAvar) )
           end
           2: Begin
              Prser2PA = 0
           end
           ENDCASE
        prior = [prior,'Prser2PA']
     endif 
     
     
  endif

  
  ;SERSIC 3 PROFILE
  if sernn3 gt 0 then begin
     
     
     ser3X0 = pars[sorder3]
     ser3Y0 = pars[sorder3+1]
     fser3x0 = flist[0,sorder3]
     fser3y0 = flist[0,sorder3+1]
     pser3x0 = flist[1,sorder3]
     pser3y0 = flist[1,sorder3+1]
     Priser3x0_temp = Priser3x0
     Priser3y0_temp = Priser3y0

     Ie3 = 10D^(pars[sorder3+2])
     fIe3 = flist[0,sorder3+2]
     pIe3 = flist[1,sorder3+2]

     Re3 = 10D^(pars[sorder3+3])
     fRe3 = flist[0,sorder3+3]
     pRe3 = flist[1,sorder3+3]
     
     n3 = (pars[sorder3+4])
     fn3 = flist[0,sorder3+4]
     pn3 = flist[1,sorder3+4]
     
     ser3E = pars[sorder3+5]
     fser3E = flist[0,sorder3+5]
     pser3E = flist[1,sorder3+5]
     
     ser3PA = (pars[sorder3+6])
     fser3PA = flist[0,sorder3+6]
     pser3PA = flist[1,sorder3+6]
     
     if fser3X0 eq 1 or fser3Y0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fser3X0 eq 1 then begin 
           CASE pser3x0 of 
              0: begin
                 Aser3X0 = (Priser3X0)-3
                 Bser3X0 = (Priser3X0/2D)+3
                 if ser3X0 ge Aser3x0 and ser3X0 le Bser3X0 then begin
                    Prser3x0 = 1D ;Uniprior(serX0, A=Aserx0, B=BserX0)
                 endif else Prser3X0 = 0.
              end
              1: Prser3x0 = Gprior(ser3X0, AVG=Priser3X0, SIG=(Priser3X0/3D))
           ENDCASE
           prior = [prior, 'Prser3x0']
        endif 
        if fser3Y0 eq 1 then begin 
           CASE pser3x0 of 
              0: begin
                 Aser3Y0 = (Priser3Y0/2D)-3
                 Bser3Y0 = (Priser3Y0/2D)+3
                 if ser3Y0 ge Aser3Y0 and ser3Y0 le Bser3Y0 then begin
                    Prser3Y0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else Prser3Y0 = 0.
              end
              1: Prser3Y0 = Gprior(ser3Y0, AVG=Priser3Y0, SIG=(s[1]/3D))
           ENDCASE
        endif 
        prior = [prior, 'Prser3y0']
     endif
    

     if fIe3 eq 1 then begin
        PriIe3 = PriIe3
        CASE pIe3 of 
           0: begin
              AIe3 = 15D ;Primue + 1D ;10D^(-1D * 0.4 * ((Primue +0.2D) - zcalib - 2.5*alog10((scale^2D))))
              BIe3 = 25D ; Primue - 1D ;10D^(-1D * 0.4 * ((Primue -0.2D) - zcalib - 2.5*alog10((scale^2D))))
              if Ie3 ge AIe3 and Ie3 le BIe3 then begin
                 PrIe3 = -1.*ALOG(BIe3 - AIe3);Uniprior(Ie, A=AIe, B=BIe)
              endif else PrIe3 = 0.
           end
           1: Begin
              Ieval3 = Ie3
              Ievar3 = 10D
              PrIe3 = -0.5*( ( ((Ieval3-PriIe3)^2D)/Ievar3) + ALOG(2D*!DPI*Ievar3) )
             
           end
           2: Begin
              AIe3 = 15D
              BIe3 = 25D
             
              if (Ie3 gt AIe3) and (Ie3 lt BIe3) then Ieval3 = 15D else Ieval3 = Ie
              
              Ievar3 = 25D
              PrIe3 = -0.5*( ( ((Ieval3-20D)^2D)/Ievar3) + ALOG(2D*!DPI*Ievar3) )
              
           end
        ENDCASE
        prior = [prior, 'PrIe3']
     endif 
     
     if fRe3 eq 1 then begin 
        PriRe3 = PriRe3
        CASE pRe3 of 
           0: BEGIN
              ARe3 = 0D ;PriRe/2D
              ;if ARe lt 1. then ARe = 1.
              BRe3 = radi ;(PriRe + ARe)
              if Re3 ge ARe3 and Re3 le BRe3 then begin
                 PrRe3 = -1.*ALOG(BRe3 - ARe3);;1D ;Uniprior(Re, A=ARe, B=BRe)
              endif else PrRe3 = 0.
           end
           1: begin

              reval3 = Re3
              Revar3 = 10D;radi/2D
              PrRe3 = -0.5*( ( ((Reval3-PriRe3)^2D)/Revar3) + ALOG(2D*!DPI*Revar3) )
           end
           2: Begin
              ARe3 = 0
              BRe3 = radi
              
              priRe3 = radi/2D
              if (Re3 lt BRe3) then Reval3 = radi else Reval3 = Re3
              
              Revar3 = radi/2D
              PrRe3 = -0.5*( ( ((Reval3-PriRe3)^2D)/Revar3) + ALOG(2D*!DPI*Revar3) )
              
           end

        ENDCASE
        prior = [prior,'PrRe3']
     endif 
     
     if fn3 eq 1 then begin 
        CASE pn3 of 
           0: begin
              An3 = 0.2
              Bn3 = 12.
              if n3 ge An3 and n3 le Bn3 then begin
                 Prn3 = -1.*ALOG(Bn3 - An3);;1D ;Uniprior(n, A=An, B=Bn)
              endif else Prn3 = 0.
           end
           1: begin
              nval3 = n
              nvar3 = 6.
              Prn3 = -0.5*( ( ((nval3-Prin3)^2D)/nvar3) + ALOG(2D*!DPI*nvar3) )

           end
           2: Begin
              An3 = 0D
              Bn3 = 12D
              
              prin3 = 6D
              if (n3 gt An3) and (n3 lt Bn3) then nval3 = 12D else nval3 = n3
              
              nvar3 = 6.
              Prn3 = -0.5*( ( ((nval3-Prin3)^2D)/nvar3) + ALOG(2D*!DPI*nvar3) )
              
           end

        ENDCASE
        prior = [prior, 'Prn3']
     endif 
     
     if fser3E eq 1 then begin 
        CASE pser3E of 
           0: begin
              Aser3E = 0D
              Bser3E = 1D
              if ser3E ge Aser3E and ser3E le Bser3E then begin
                 Prser3E = -1.*ALOG(Bser3E - Aser3E);;1D ;Uniprior(serE, A=AserE, B=BserE)
              endif else Prser3E = 0.
           end
           1: begin
              ser3Evar = 2D;mean((serE - PriserE)^2D)
              Prser3E = -0.5*( ( ((ser3E-Priser3E)^2D)/ser3Evar) + ALOG(2D*!DPI*ser3Evar) )
              if ser3E gt 1 or ser3E lt 0. then Prser3E = 0.
           end
           2: Begin
              Prser3E = 0
           end

        ENDCASE
        prior = [prior,'Prser3E']
     endif 
     
     if fser3PA eq 1 then begin 
        CASE pser3PA of 
           0: begin
              Aser3PA = Priser3PA - (priser3PA/2D)
              Bser3PA = Priser3PA + (priser3PA/2D)
              if ser3PA ge Aser3Pa and ser3Pa le Bser3Pa then begin 
                 Prser3PA = -1.*ALOG(Bser3PA - Aser3PA);;1D ;Uniprior(serPA, A=AserPA, B=BserPA)
              endif else Prser3Pa = 0.
           end
           1: begin
              ser3PAvar = 1D    ;mean((serPA - PriserPA)^2D)
              Prser3PA = -0.5*( ( ((ser3PA-Priser3PA)^2D)/ser3PAvar) + ALOG(2D*!DPI*ser3PAvar) )
           end
           2: Begin
              Prser3PA = 0
           end
           ENDCASE
        prior = [prior,'Prser3PA']
     endif 
     
     
  endif

  
  ;SERSIC 4 PROFILE
  if sernn4 gt 0 then begin
     
     
     ser4X0 = pars[sorder4]
     ser4Y0 = pars[sorder4+1]
     fser4x0 = flist[0,sorder4]
     fser4y0 = flist[0,sorder4+1]
     pser4x0 = flist[1,sorder4]
     pser4y0 = flist[1,sorder4+1]
     Priser4x0_temp = Priser4x0
     Priser4y0_temp = Priser4y0

     Ie4 = 10D^(pars[sorder4+2])
     fIe4 = flist[0,sorder4+2]
     pIe4 = flist[1,sorder4+2]

     Re4 = 10D^(pars[sorder4+3])
     fRe4 = flist[0,sorder4+3]
     pRe4 = flist[1,sorder4+3]
     
     n4 = (pars[sorder4+4])
     fn4 = flist[0,sorder4+4]
     pn4 = flist[1,sorder4+4]
     
     ser4E = pars[sorder4+5]
     fser4E = flist[0,sorder4+5]
     pser4E = flist[1,sorder4+5]
     
     ser4PA = (pars[sorder4+6])
     fser4PA = flist[0,sorder4+6]
     pser4PA = flist[1,sorder4+6]
     
     if fser4X0 eq 1 or fser4Y0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fser4X0 eq 1 then begin 
           CASE pser4x0 of 
              0: begin
                 Aser4X0 = (Priser4X0)-3
                 Bser4X0 = (Priser4X0/2D)+3
                 if ser4X0 ge Aser4x0 and ser4X0 le Bser4X0 then begin
                    Prser4x0 = 1D ;Uniprior(serX0, A=Aserx0, B=BserX0)
                 endif else Prser4X0 = 0.
              end
              1: Prser4x0 = Gprior(ser4X0, AVG=Priser4X0, SIG=(Priser4X0/3D))
           ENDCASE
           prior = [prior, 'Prser4x0']
        endif 
        if fser4Y0 eq 1 then begin 
           CASE pser4x0 of 
              0: begin
                 Aser4Y0 = (Priser4Y0/2D)-3
                 Bser4Y0 = (Priser4Y0/2D)+3
                 if ser4Y0 ge Aser4Y0 and ser4Y0 le Bser4Y0 then begin
                    Prser4Y0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else Prser4Y0 = 0.
              end
              1: Prser4Y0 = Gprior(ser4Y0, AVG=Priser4Y0, SIG=(s[1]/3D))
           ENDCASE
        endif 
        prior = [prior, 'Prser4y0']
     endif
    

     if fIe4 eq 1 then begin
        PriIe4 = PriIe4
        CASE pIe4 of 
           0: begin
              AIe4 = 15D ;Primue + 1D ;10D^(-1D * 0.4 * ((Primue +0.2D) - zcalib - 2.5*alog10((scale^2D))))
              BIe4 = 25D ; Primue - 1D ;10D^(-1D * 0.4 * ((Primue -0.2D) - zcalib - 2.5*alog10((scale^2D))))
              if Ie4 ge AIe4 and Ie4 le BIe4 then begin
                 PrIe4 = -1.*ALOG(BIe4 - AIe4);Uniprior(Ie, A=AIe, B=BIe)
              endif else PrIe4 = 0.
           end
           1: Begin
              Ieval4 = Ie4
              Ievar4 = 10D
              PrIe4 = -0.5*( ( ((Ieval4-PriIe4)^2D)/Ievar4) + ALOG(2D*!DPI*Ievar4) )
             
           end
           2: Begin
              AIe4 = 15D
              BIe4 = 25D
             
              if (Ie4 gt AIe4) and (Ie4 lt BIe4) then Ieval4 = 15D else Ieval4 = Ie4
              
              Ievar4 = 25D
              PrIe4 = -0.5*( ( ((Ieval4-20D)^2D)/Ievar4) + ALOG(2D*!DPI*Ievar4) )
              
           end
        ENDCASE
        prior = [prior, 'PrIe4']
     endif 
     
     if fRe4 eq 1 then begin 
        PriRe4 = PriRe4
        CASE pRe4 of 
           0: BEGIN
              ARe4 = 0D ;PriRe/2D
              ;if ARe lt 1. then ARe = 1.
              BRe4 = radi ;(PriRe + ARe)
              if Re4 ge ARe4 and Re4 le BRe4 then begin
                 PrRe4 = -1.*ALOG(BRe4 - ARe4);;1D ;Uniprior(Re, A=ARe, B=BRe)
              endif else PrRe4 = 0.
           end
           1: begin

              reval4 = Re4
              Revar4 = 10D;radi/2D
              PrRe4 = -0.5*( ( ((Reval4-PriRe4)^2D)/Revar4) + ALOG(2D*!DPI*Revar4) )
           end
           2: Begin
              ARe4 = 0
              BRe4 = radi
              
              priRe4 = radi/2D
              if (Re4 lt BRe4) then Reval4 = radi else Reval4 = Re4
              
              Revar4 = radi/2D
              PrRe4 = -0.5*( ( ((Reval4-PriRe4)^2D)/Revar4) + ALOG(2D*!DPI*Revar4) )
              
           end

        ENDCASE
        prior = [prior,'PrRe4']
     endif 
     
     if fn4 eq 1 then begin 
        CASE pn4 of 
           0: begin
              An4 = 0.2
              Bn4 = 12.
              if n4 ge An4 and n4 le Bn4 then begin
                 Prn4 = -1.*ALOG(Bn4 - An4);;1D ;Uniprior(n, A=An, B=Bn)
              endif else Prn4 = 0.
           end
           1: begin
              nval4 = n4
              nvar4 = 6.
              Prn4 = -0.5*( ( ((nval4-Prin4)^2D)/nvar4) + ALOG(2D*!DPI*nvar4) )

           end
           2: Begin
              An4 = 0D
              Bn4 = 12D
              
              prin4 = 6D
              if (n4 gt An4) and (n4 lt Bn4) then nval4 = 12D else nval4 = n
              
              nvar4 = 6.
              Prn4 = -0.5*( ( ((nval4-Prin4)^2D)/nvar4) + ALOG(2D*!DPI*nvar4) )
              
           end

        ENDCASE
        prior = [prior, 'Prn4']
     endif 
     
     if fser4E eq 1 then begin 
        CASE pser4E of 
           0: begin
              Aser4E = 0D
              Bser4E = 1D
              if ser4E ge Aser4E and ser4E le Bser4E then begin
                 Prser4E = -1.*ALOG(Bser4E - Aser4E);;1D ;Uniprior(serE, A=AserE, B=BserE)
              endif else Prser4E = 0.
           end
           1: begin
              ser4Evar = 2D;mean((serE - PriserE)^2D)
              Prser4E = -0.5*( ( ((ser4E-Priser4E)^2D)/ser4Evar) + ALOG(2D*!DPI*ser4Evar) )
              if ser4E gt 1 or serE lt 0. then Prser4E = 0.
           end
           2: Begin
              Prser4E = 0
           end

        ENDCASE
        prior = [prior,'Prser4E']
     endif 
     
     if fser4PA eq 1 then begin 
        CASE pser4PA of 
           0: begin
              Aser4PA = Priser4PA - (priser4PA/2D)
              Bser4PA = Priser4PA + (priser4PA/2D)
              if ser4PA ge Aser4Pa and ser4Pa le Bser4Pa then begin 
                 Prser4PA = -1.*ALOG(Bser4PA - Aser4PA);;1D ;Uniprior(serPA, A=AserPA, B=BserPA)
              endif else Prser4Pa = 0.
           end
           1: begin
              ser4PAvar = 1D    ;mean((serPA - PriserPA)^2D)
              Prser4PA = -0.5*( ( ((ser4PA-Priser4PA)^2D)/ser4PAvar) + ALOG(2D*!DPI*ser4PAvar) )
           end
           2: Begin
              Prser4PA = 0
           end
           ENDCASE
        prior = [prior,'Prser4PA']
     endif 
     
     
  endif
  

  
  if expnn gt 0 then begin
     
     
     expX0 = pars[eorder]
     expY0 = pars[eorder+1]
     fexpx0 = flist[0,eorder]
     fexpy0 = flist[0,eorder+1]
     pexpx0 = flist[1,eorder]
     pexpy0 = flist[1,eorder+1]
     Priexpx0_temp = Priexpx0
     Priexpy0_temp = Priexpy0
  
     I0 = 10D^(pars[eorder+2])
     fI0 = flist[0,eorder+2]
     pI0 = flist[1,eorder+2]

     h = 10D^(pars[eorder+3])
     fh = flist[0,eorder+3]
     ph = flist[1,eorder+3]
     
     expE = pars[eorder+4]
     fexpE = flist[0,eorder+4]
     pexpE = flist[1,eorder+4]
     
     expPA = (pars[eorder+5]);*180D/!DPI
     fexpPA = flist[0,eorder+5]
     pexpPA = flist[1,eorder+5]
     
     if fexpX0 eq 1 or fexpY0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fexpX0 eq 1 then begin 
           CASE pexpx0 of 
              0: begin
                 AexpX0 = (PriexpX0)-3
                 BexpX0 = (PriexpX0/2D)+3
                 if expX0 ge Aexpx0 and expX0 le BexpX0 then begin
                    Prexpx0 = 1D ;Uniprior(expX0, A=Aexpx0, B=BexpX0)
                 endif else PrexpX0 = 0.
              end
              1: BEGIN
                 Prexpx0 = -0.5*ALOG(2D*!DPI*serPAvar) + 0.5*(n-PriserPA)*(n-PriserPA)/serPAvar
              end
           ENDCASE
           prior = [prior, 'Prexpx0']
        endif 
        if fexpY0 eq 1 then begin 
           CASE pexpx0 of 
              0: begin
                 AexpY0 = (PriexpY0/2D)-3
                 BexpY0 = (PriexpY0/2D)+3
                 if expY0 ge AexpY0 and expY0 le BexpY0 then begin
                    PrexpY0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else PrexpY0 = 0.
              end
              1: PrexpY0 = Gprior(expY0, AVG=PriexpY0, SIG=(s[1]/3D))
              
           ENDCASE
        endif 
        prior = [prior, 'Prexpy0']
     endif

     if fI0 eq 1 then begin 
        CASE pI0 of 
           0: begin
              AI0 = PriI0 +0.5D ;10D^(-1D * 0.4 * ((Primu0 +0.2D) - zcalib - 2.5*alog10((scale^2D))))
              BI0 = PriI0 -0.5D ;10D^(-1D * 0.4 * ((Primu0 -0.2D) - zcalib - 2.5*alog10((scale^2D))))
              if I0 ge AI0 and I0 le BI0 then begin
                 PrI0 = -1.*ALOG(BI0 - AI0);;1D ;Uniprior(Ie, A=AIe, B=BIe)
              endif else PrI0 = 0.
           end
           1: Begin            
              I0val = I0
              I0var = 10D;25.
              PrI0 = -0.5*( ( ((I0val-20D)^2D)/I0var) + ALOG(2D*!DPI*I0var) )
           end
           2: Begin
              AI0 = 15D
              BI0 = 25
             
              if (I0 ge AI0) and (I0 le BI0) then I0val = 15D else I0val = I0
              
              I0var = 25D
              PrI0 = -0.5D*( ( ((I0val-20D)^2D)/I0var) + ALOG(2D*!DPI*I0var) )
           end
              
        ENDCASE
        prior = [prior, 'PrI0']
     endif 

     if fh eq 1 then begin 
        CASE ph of 
           0: BEGIN
              Ah = Prih/2D
              if Ah lt 1. then Ah = 1.
              Bh = (Prih - Ah) + Prih
              if h ge Ah and h le Bh then begin
                 Prh = -1.*ALOG(Bh - Ah);;1D ;Uniprior(Re, A=ARe, B=BRe)
              endif else Prh = 0.
           end
           1: begin
              
              hval = h
              hvar = 10D;radi/2D
              Prh = -0.5D*( ( ((hval-Prih)^2D)/hvar) + ALOG(2D*!DPI*hvar) )
           end
           2: Begin
              Ah = 0
              Bh = radi
              
              prih = radi/2D
              if (h lt Bh) then hval = radi else hval = h

              hvar = radi/2D
              Prh = -0.5*( ( ((hval-Prih)^2D)/hvar) + ALOG(2D*!DPI*hvar) )
           end
        ENDCASE
        prior = [prior,'Prh']
     endif 
          
     if fexpE eq 1 then begin 
        CASE pexpE of 
           0: begin
              AexpE = 0.
              BexpE = 1.
              if expE ge AexpE and expE le BexpE then begin
                 PrexpE = -1.*ALOG(BexpE - AexpE);;1D ;Uniprior(expE, A=AexpE, B=BexpE)
              endif else PrexpE = 0.
           end
           1: begin
              expEvar = 2D
              PrexpE = -0.5*( ( ((expE-PriserE)^2D)/expEvar) + ALOG(2D*!DPI*expEvar) )
              if expE gt 1 or expE lt 0. then PrexpE = 0.
           end
           2: Begin
              PrexpE = 0
           end
        ENDCASE
        prior = [prior,'PrexpE']
     endif 
     
     if fexpPA eq 1 then begin 
        CASE pexpPA of 
           0: begin
              AexpPA = PriexpPA - (priexpPA/2D)
              BexpPA = PriexpPA + (priexpPA/2D)
              if expPA ge AexpPa and expPa le BexpPa then begin 
                 PrexpPA = -1.*ALOG(BexpPA - AexpPA);;1D ;Uniprior(expPA, A=AexpPA, B=BexpPA)
              endif else PrexpPa = 0.
           end
           1: begin
              expPAvar = 45D
              PrexpPA = -0.5*( ( ((expPA-PriexpPA)^2D)/expPAvar) + ALOG(2D*!DPI*expPAvar) )
           end
           2: Begin
              PrexpPA = 0
           end
        ENDCASE
        prior = [prior,'PrexpPA']
     endif 
     
  endif 
  
  if (cennn eq 1) then begin
     
     Prix0 = mean([Priserx0_temp,Priexpx0_temp],/NAN)
     Priy0 = mean([Prisery0_temp,Priexpy0_temp],/NAN)
     
  endif
  
  if cennn gt 0 then begin
     
     X0 = pars[cenorder]
     Y0 = pars[cenorder+1]

     fx0 = flist[0,cenorder]
     fy0 = flist[0,cenorder+1]
     px0 = flist[1,cenorder]
     py0 = flist[1,cenorder+1]

      if fX0 eq 1 or fY0 eq 1 then begin
        s = DOUBLE(size(ima,/DIMENSIONS))
        if fX0 eq 1 then begin 
           CASE px0 of 
              0: begin
                 AX0 = (PriX0)-3
                 BX0 = (PriX0/2D)+3
                 if X0 ge Ax0 and X0 le BX0 then begin
                    Prx0 = 1D ;Uniprior(expX0, A=Aexpx0, B=BexpX0)
                 endif else PrX0 = 0.
              end
              1: BEGIN
              end
           ENDCASE
           prior = [prior, 'Prx0']
        endif 
        if fY0 eq 1 then begin 
           CASE py0 of 
              0: begin
                 AY0 = (PriY0/2D)-3
                 BY0 = (PriY0/2D)+3
                 if Y0 ge AY0 and Y0 le BY0 then begin
                    PrY0 = 1D ;Uniprior(serY0, A=AserY0, B=BserY0)
              endif else PrY0 = 0.
              end
              1: PrY0 = Gprior(Y0, AVG=PriY0, SIG=(s[1]/3D))
              
           ENDCASE
        endif 
        prior = [prior, 'Pry0']
     endif

   endif

  if skynn gt 0 then begin
     
     sky = pars[skyorder]
     fsky = flist[0,skyorder]
     psky = flist[1,skyorder]

  endif
  
  REMOVE,0,prior
  pri = strjoin(prior,'+')
  
  ;excprior = execute('pri = ' + pri)
return,pri 
END

;Functions for the Newton-raphson method
function func1,pars
  x = pars[0]
  ie = pars[1]
  re = pars[2]
  n = pars[3]
  i0 = pars[4]
  h = pars[5]

  fexp = (I0*exp(-x/h))
  
  invn = 1D / n
  bn = (1.9992 * n) - 0.3271
  
  fser = Ie * EXP(-bn * ( ((x/Re)^(invn)) - 1D) )
  
  f = fser - fexp
  return, f
end

function func2,pars
  x = pars[0]
  ie = pars[1]
  re = pars[2]
  n = pars[3]
  i0 = pars[4]
  h = pars[5]

  fexp = (I0*exp(-x/h)) / h
  
  invn = 1D / n
  bn = (1.9992 * n) - 0.3271
  
  fser = Ie * bn * EXP(-bn * ( ((x/Re)^(invn)) - 1D) ) * ((x/Re)^(invn)) / (n * x)

  f = fexp - fser
  return, f
end

;Newton-Raphson method
function cross, x,res,func=func,ffunc=ffunc,iter=iter,pars=pars
  
  oldres = res
  ;res = ABS(res)
  if (res[0] lt 0) then index = where(res ge 0 ,nn) else  index = where(res le 0 ,nn) 
  crX = [0D,0D]

  if (nn le 0) then begin
     ;return, crX                ;No crossing points 
     goto, theEnd
  endif 
  
  ;Perform Newton-Raphson 
  x0 = x[min(index)]
  
  if (NOT keyword_set(iter)) then iter = 10
  count = 0
  f_pars = [x0,pars]
  f0 = call_function(func,f_pars)
  ff0 = call_function(ffunc,f_pars)
  xn = ABS(x0 - (f0/ff0))
  fn = f0
  ffn = ff0
  delx = ABS(xn - x0)
  xold = x0
  
  if (ABS(fn) GT 0.001D) then begin
     REPEAT BEGIN
        f_pars = double([xold,pars])  
        fn = call_function(func,f_pars)
        ffn = call_function(ffunc,f_pars)
        xn = ABS(xold - (fn/ffn))
        delx = ABS(xn - xold)
        ;print, xold,xn,fn,ffn,delx
        xold = xn
        count = count + 1
     ENDREP UNTIL (ABS(fn) LT 0.001D) or (count gt 100)
  endif 
  crX[0] = xn                   ;First crossing point
  if (count ge 100) then print, 'N-R thingy is being shit'
  
  newx = x[min(index):*]
  newres = res[min(index):*]
  if (res[0] lt 0) then newindex = where(newres le 0,nn) else newindex = where(newres ge 0,nn) 
  if (nn le 0) then goto, theEnd ;One crossing points 
  
  ;Perform Newton-Raphson 
  x0 = newx[min(newindex)]
  

  count = 0
  f_pars = [x0,pars]
  f0 = call_function(func,f_pars)
  ff0 = call_function(ffunc,f_pars)
  xn = ABS(x0 - (f0/ff0))
  fn = f0
  ffn = ff0
  
  delx = ABS(xn - x0)
  xold = x0
  if (ABS(fn) GT 0.001D) then begin
     REPEAT BEGIN
        f_pars = [xold,pars]
        fn = call_function(func,f_pars)
        ffn = call_function(ffunc,f_pars)
        xn = ABS(xold - (fn/ffn))
        delx = ABS(xn - xold)
        xold = xn
        ;print, xold,xn,fn,ffn,delx
        count = count + 1
     ENDREP UNTIL ABS(fn) LT 0.001D or (count gt 100)
     if (count ge 100) then print, 'N-R thingy is being shit'
  endif 
  crX[1] = xn                ;Second crossing point


  theEnd:
  res = oldres
  return, crX   ;Two crossing points 
end

function bagal_flag, pars
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  sorder2 = where(nlist eq 'ser2X0', sernn2)
  sorder3 = where(nlist eq 'ser3X0', sernn3)
  sorder4 = where(nlist eq 'ser4X0', sernn4)
  beorder = where(nlist eq 'bexpX0', bexpnn)
  forder = where(nlist eq 'ferX0', fernn)
  forder2 = where(nlist eq 'fer2X0', fernn2)
  
  cenorder = where(nlist eq 'X0', cennn)

  if sernn le 0 then sflag = 1
  if sernn2 le 0 then sflag2 = 1
  if sernn3 le 0 then sflag3 = 1
  if sernn4 le 0 then sflag4 = 1
  if expnn le 0 then eflag = 1
  if bexpnn le 0 then beflag = 1
  if fernn le 0 then fflag = 1
  if fernn2 le 0 then fflag2 = 1
  if cennn le 0 then cenflag = 1

  flag = 1
  ;SERSIC PROFILE
  if sernn gt 0 then begin
     sflag = 1
     ;Centre coordinates x0 and y0
     serX0 = pars[sorder]
     serY0 = pars[sorder+1]
     
     if serx0 lt ((nx/2D) - 20) then sflag = 0
     if serx0 gt ((nx/2D) + 20) then sflag = 0
     
     if sery0 lt ((ny/2D) - 20) then sflag = 0
     if sery0 gt ((ny/2D) + 20) then sflag = 0
     
     
     ;Effective intensity 
     Ie = (10D^(pars[sorder+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if Ie lt 0D then sflag = 0
     
     ;Effective Radius 
     Re = (10D^(pars[sorder+3]))[0]
     if Re lt 0.01D then sflag = 0
     ;if Re ge 8 then sflag = 0
     ;if (Re ge nx-1) or (Re ge ny-1) then sflag = 0 
     
     ;Sersic index
     n = (10D^pars[sorder+4])[0]
     if n lt 0.4 then sflag = 0
     if n gt 8. then sflag = 0
     
     ;Sersic b/a 
     serE = pars[sorder+5]   
     if serE lt 0.2D then sflag = 0
     if serE gt 1D then sflag = 0
     
     ;Sersic PA
     ;serPA = ACOS(pars[sorder+6])*180D/!DPI 
     serPA = (pars[sorder+6]) 
     if serPa lt -360D then sflag = 0
     if serPA gt 360D then sflag = 0
     
     
     ;if (sflag eq 0) then stop
  endif
  
  ;SERSIC 2 PROFILE
  if sernn2 gt 0 then begin
     sflag2 = 1
     ;Centre coordinates x0 and y0
     ser2X0 = pars[sorder2]
     ser2Y0 = pars[sorder2+1]
     
     if ser2x0 lt ((nx/2D) - 20) then sflag2 = 0
     if ser2x0 gt ((nx/2D) + 20) then sflag2 = 0
     
     if ser2y0 lt ((ny/2D) - 20) then sflag2 = 0
     if ser2y0 gt ((ny/2D) + 20) then sflag2 = 0
     
     
     ;Effective intensity 
     Ie2 = (10D^(pars[sorder2+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if Ie2 lt 0D then sflag2 = 0
     
     ;Effective Radius 
     Re2 = (10D^(pars[sorder2+3]))[0]
     if Re2 lt 0.01D then sflag2 = 0
     ;if Re ge 8 then sflag = 0
     ;if (Re ge nx-1) or (Re ge ny-1) then sflag = 0 
     
     ;Sersic index
     n2 = (10D^pars[sorder2+4])[0]
     if n2 lt 0.4 then sflag2 = 0
     if n2 gt 8. then sflag2 = 0
     
     ;Sersic b/a 
     ser2E = pars[sorder2+5]   
     if ser2E lt 0.2D then sflag2 = 0
     if ser2E gt 1D then sflag2 = 0
     
     ;Sersic PA
     ;serPA = ACOS(pars[sorder+6])*180D/!DPI 
     ser2PA = (pars[sorder2+6]) 
     if ser2Pa lt -360D then sflag2 = 0
     if ser2PA gt 360D then sflag2 = 0
     
     
     ;if (sflag eq 0) then stop
  endif
  
  ;SERSIC 3 PROFILE
  if sernn3 gt 0 then begin
     sflag3 = 1
     ;Centre coordinates x0 and y0
     ser3X0 = pars[sorder3]
     ser3Y0 = pars[sorder3+1]
     
     if ser3x0 lt ((nx/2D) - 20) then sflag3 = 0
     if ser3x0 gt ((nx/2D) + 20) then sflag3 = 0
     
     if ser3y0 lt ((ny/2D) - 20) then sflag3 = 0
     if ser3y0 gt ((ny/2D) + 20) then sflag3 = 0
     
     
     ;Effective intensity 
     Ie3 = (10D^(pars[sorder3+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if Ie3 lt 0D then sflag3 = 0
     
     ;Effective Radius 
     Re3 = (10D^(pars[sorder3+3]))[0]
     if Re3 lt 0.01D then sflag3 = 0
     ;if Re ge 8 then sflag = 0
     ;if (Re ge nx-1) or (Re ge ny-1) then sflag = 0 
     
     ;Sersic index
     n3 = (10D^pars[sorder3+4])[0]
     if n3 lt 0.4 then sflag3 = 0
     if n3 gt 8. then sflag3 = 0
     
     ;Sersic b/a 
     ser3E = pars[sorder3+5]   
     if ser3E lt 0.2D then sflag3 = 0
     if ser3E gt 1D then sflag3 = 0
     
     ;Sersic PA
     ;serPA = ACOS(pars[sorder+6])*180D/!DPI 
     ser3PA = (pars[sorder3+6]) 
     if ser3Pa lt -360D then sflag3 = 0
     if ser3PA gt 360D then sflag3 = 0
     
     
     ;if (sflag eq 0) then stop
  endif
  
  ;SERSIC 4 PROFILE
  if sernn4 gt 0 then begin
     sflag4 = 1
     ;Centre coordinates x0 and y0
     ser4X0 = pars[sorder4]
     ser4Y0 = pars[sorder4+1]
     
     if ser4x0 lt ((nx/2D) - 20) then sflag4 = 0
     if ser4x0 gt ((nx/2D) + 20) then sflag4 = 0
     
     if ser4y0 lt ((ny/2D) - 20) then sflag4 = 0
     if ser4y0 gt ((ny/2D) + 20) then sflag4 = 0
     
     
     ;Effective intensity 
     Ie4 = (10D^(pars[sorder4+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if Ie4 lt 0D then sflag4 = 0
     
     ;Effective Radius 
     Re4 = (10D^(pars[sorder4+3]))[0]
     if Re4 lt 0.01D then sflag4 = 0
     ;if Re ge 8 then sflag = 0
     ;if (Re ge nx-1) or (Re ge ny-1) then sflag = 0 
     
     ;Sersic index
     n4 = (10D^pars[sorder4+4])[0]
     if n4 lt 0.4 then sflag4 = 0
     if n4 gt 8. then sflag4 = 0
     
     ;Sersic b/a 
     ser4E = pars[sorder4+5]   
     if ser4E lt 0.2D then sflag4 = 0
     if ser4E gt 1D then sflag4 = 0
     
     ;Sersic PA
     ;serPA = ACOS(pars[sorder+6])*180D/!DPI 
     ser4PA = (pars[sorder4+6]) 
     if ser4Pa lt -360D then sflag4 = 0
     if ser4PA gt 360D then sflag4 = 0
     
     
     ;if (sflag eq 0) then stop
  endif

 ;Exponential 
  if expnn gt 0 then begin
     expX0 = pars[eorder]
     expY0 = pars[eorder+1]
     eflag = 1
     if expx0 lt ((nx/2D) - 20) then eflag = 0
     if expx0 gt ((nx/2D) + 20) then eflag = 0
     
     if expy0 lt ((ny/2D) - 20) then eflag = 0
     if expy0 gt ((ny/2D) + 20) then eflag = 0
     
     ;Central intensity 
     I0 = (10D^(pars[eorder+2]))[0]
     ;if mu0 gt 30 then efalg = 0
     if I0 lt 0. then eflag = 0

     ;Scale height 
     h = (10D^(pars[eorder+3]))[0]
     if h lt 0.05D then eflag = 0
     ;if h ge (2D*radi) then eflag = 0
     if (h ge nx-1) or (h ge ny-1) then eflag = 0

     ;Exponential b/a 
     expE = pars[eorder+4]   
     if expE lt 0.2D then eflag = 0
     if expE gt 1D then eflag = 0

     ;Exponential PA
     ;expPA = ACOS(pars[eorder+5])*180D/!DPI
     expPA = (pars[eorder+5])
     if expPa lt -360D then eflag = 0
     if expPA gt 360D then eflag = 0
     
     ;if (eflag eq 0) then stop
  endif

  ;Broken Exponential 
  if bexpnn gt 0 then begin
     bexpX0 = pars[beorder]
     bexpY0 = pars[beorder+1]
     beflag = 1
     if bexpx0 lt ((nx/2D) - 20) then beflag = 0
     if bexpx0 gt ((nx/2D) + 20) then beflag = 0
     
     if bexpy0 lt ((ny/2D) - 20) then beflag = 0
     if bexpy0 gt ((ny/2D) + 20) then beflag = 0
     
     ;Central intensity 
     bI0 = (10D^(pars[beorder+2]))[0]
     ;if mu0 gt 30 then efalg = 0
     if bI0 lt 0. then beflag = 0

     ;Scale height 1  
     h1 = (10D^(pars[beorder+3]))[0]
     if h1 lt 0.05D then beflag = 0
     ;if h ge (2D*radi) then eflag = 0
     if (h1 ge nx-1) or (h1 ge ny-1) then beflag = 0

     ;Scale height 2  
     h2 = (10D^(pars[beorder+4]))[0]
     if h2 lt 0.05D then beflag = 0
     if h2 lt h1 then beflag = 0
     ;if h ge (2D*radi) then eflag = 0
     if (h2 ge nx-1) or (h2 ge ny-1) then beflag = 0
     
     ;Shapeness  
     bexpa = pars[beorder+5]   
     ;if bexpa lt 0.2D then beflag = 0
     ;if bexpa gt 1D then beflag = 0

     ;Brake radius 
     rb = (10D^pars[beorder+6])[0]
     if rb lt h1 then beflag = 0
     if rb gt h2 then beflag = 0
    
     ;Exponential b/a 
     bexpE = pars[beorder+7]   
     if bexpE lt 0.2D then beflag = 0
     if bexpE gt 1D then beflag = 0

     ;Exponential PA
     ;expPA = ACOS(pars[eorder+5])*180D/!DPI
     bexpPA = (pars[beorder+8])
     if bexpPa lt -360D then beflag = 0
     if bexpPA gt 360D then beflag = 0
     
     ;if (eflag eq 0) then stop
  endif
 
  if fernn gt 0 then begin
     fflag = 1
     ;Centre coordinates x0 and y0
     ferX0 = pars[forder]
     ferY0 = pars[forder+1]
     
     if ferx0 lt ((nx/2D) - 20) then fflag = 0
     if ferx0 gt ((nx/2D) + 20) then fflag = 0
     
     if fery0 lt ((ny/2D) - 20) then fflag = 0
     if fery0 gt ((ny/2D) + 20) then fflag = 0
     
     
     ;Central intensity 
     ferI0 = (10D^(pars[forder+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if ferI0 lt 0D then fflag = 0
     
     ;Length 
     l = (10D^(pars[forder+3]))[0]
     if l lt 0.01D then fflag = 0
     
     ;Shape parameter
     fern = (10D^pars[forder+4])[0]
     if fern lt 0.4 then fflag = 0
     if fern gt 8. then fflag = 0
     
     ;Ferrors b/a 
     ferE = pars[forder+5]   
     if ferE lt 0.1D then fflag = 0
     if ferE gt 1D then fflag = 0
     
     ;Ferrors PA
     ferPA = (pars[forder+6]) 
     if ferPa lt -360D then fflag = 0
     if ferPA gt 360D then fflag = 0

     ;Ferrors c
     ferc = pars[forder+7]   
     if ferc lt -1 then fflag = 0
     if ferc gt 1D then fflag = 0
     
     ;if (sflag eq 0) then stop
  endif

  if fernn2 gt 0 then begin
     fflag2 = 1
     ;Centre coordinates x0 and y0
     fer2X0 = pars[forder2]
     fer2Y0 = pars[forder2+1]
     
     if fer2x0 lt ((nx/2D) - 20) then fflag2 = 0
     if fer2x0 gt ((nx/2D) + 20) then fflag2 = 0
     
     if fer2y0 lt ((ny/2D) - 20) then fflag2 = 0
     if fer2y0 gt ((ny/2D) + 20) then fflag2 = 0
     
     
     ;Central intensity 
     fer2I0 = (10D^(pars[forder2+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if fer2I0 lt 0D then fflag2 = 0
     
     ;Length 
     l2 = (10D^(pars[forder2+3]))[0]
     if l2 lt 0.01D then fflag2 = 0
     
     ;Shape parameter
     fer2n = (10D^pars[forder2+4])[0]
     if fer2n lt 0.4 then fflag2 = 0
     if fer2n gt 8. then fflag2 = 0
     
     ;Ferrors b/a 
     fer2E = pars[forder2+5]   
     if fer2E lt 0.1D then fflag2 = 0
     if fer2E gt 1D then fflag2 = 0
     
     ;Ferrors PA
     fer2PA = (pars[forder2+6]) 
     if fer2Pa lt -360D then fflag2 = 0
     if fer2PA gt 360D then fflag2 = 0

     ;Ferrors c
     fer2c = pars[forder2+7]   
     if fer2c lt -1D then fflag2 = 0
     if fer2c gt 1D then fflag2 = 0
     
     ;if (sflag eq 0) then stop
  endif
 
  
  
  if cennn gt 0 then begin
     
     X0 = pars[cenorder]
     Y0 = pars[cenorder+1]
     cenflag = 1
     if x0 lt ((nx/2D) - 20) then cenflag = 0
     if x0 gt ((nx/2D) + 20) then cenflag = 0
     
     if y0 lt ((ny/2D) - 20) then cenflag = 0
     if y0 gt ((ny/2D) + 20) then cenflag = 0

     ;if (cenflag eq 0) then stop
  endif


  Bndflag = 1
  ;priorshere
  if (sernn gt 0) and (expnn gt 0) then begin
     
     bn = (1.9992D * n) - 0.3271D
     Lser = (2D * !dpi * Ie * (Re^2D) * n * (EXP(bn)) * gamma(2D *n)) / (bn^(2D * n))
     Lexp = (2D * !dpi * I0 * (h^2D))
     BD = Lser / Lexp
     hRe = Re / h
     BD_1Re = BD / (2D * igamma(2D,hRe))


     ;Bulge-to-disc ratio
     if (BD lt 0.01D) then Bndflag = 0
     if (BD gt 9D) then Bndflag = 0
     
     ;Re to h ratio
    
     if (Re gt (1.678 * h)) then Bndflag = 0
     if (Re lt (0.05D * h)) then Bndflag = 0
     

     ;Bulge-to-disc ratio < Re
     if (BD_1Re lt 0.1) then Bndflag = 0
     
     ;if Bndflag eq 0 then stop
     ;Double crossing points 
     ;x = [1D:200D]
     ;dexp = (I0*exp(-x/h))
   
     ;invn = 1D / n
     ;bn = (1.9992 * n) - 0.3271     
     ;dser = Ie * EXP(-bn * ( ((x/Re)^(invn)) - 1D) )
     
     ;pars2 = [Ie,Re,n,I0,h]
     ;res = dser - dexp
     ;crX = cross(x,res,func='func1',ffunc='func2',pars=pars2)
     
     ;if (crX[1] lt 7D*Re) and (crX[1] ne 0) then Bndflag = 0
     ;if (TOTAL(crX) eq 0) or (crX[0] eq 0) then Bndflag = 0
     ;if (TOTAL(crX) eq 0) then Bndflag = 0

     ;if (TOTAL(crX) eq 0) or (crX[0] gt 7D*Re) then Bndflag = 0
     ;if (crX[1] lt 10D*Re) and (crX[0] lt 10D*Re) and (crX[0] ne 0) and (crX[1] ne 0) then stop
     ;if BnDflag eq 0 then stop
     
     
  endif

  if (sflag eq 0) or (sflag2 eq 0) or (sflag3 eq 0) or (sflag4 eq 0) or (eflag eq 0) or (beflag eq 0) or (fflag eq 0) or (fflag2 eq 0) or (BnDflag eq 0) or (cenflag eq 0) then flag = 0
  ;if flag eq 0 and hre gt 1.05 then stop
  return, flag
END

function modify_initial, pars
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  sorder2 = where(nlist eq 'ser2X0', sernn2)
  sorder3 = where(nlist eq 'ser3X0', sernn3)
  sorder4 = where(nlist eq 'ser4X0', sernn4)
  cenorder = where(nlist eq 'X0', cennn)

  if sernn le 0 then sflag = 1
  if sernn2 le 0 then sflag2 = 1
  if sernn3 le 0 then sflag3 = 1
  if sernn4 le 0 then sflag4 = 1
  if expnn le 0 then eflag = 1
  if cennn le 0 then cenflag = 1

  ;flag = 1
  flag = [1,1,1,1,1,1,1,1]
  ;SERSIC PROFILE
  if sernn gt 0 then begin
     sflag = 1
     ;Centre coordinates x0 and y0
     serX0 = pars[sorder]
     serY0 = pars[sorder+1]
     
     if serx0 lt ((nx/2D) - 20) then sflag = 0
     if serx0 gt ((nx/2D) + 20) then sflag = 0
     
     if sery0 lt ((ny/2D) - 20) then sflag = 0
     if sery0 gt ((ny/2D) + 20) then sflag = 0
     
     
     ;Effective intensity 
     Ie = (10D^(pars[sorder+2]))[0]
     ;if Ie gt 6. then sflag = 0
     if Ie lt 0D then sflag = 0
     
     ;Effective Radius 
     Re = (10D^(pars[sorder+3]))[0]
     if Re lt 0.01D then sflag = 0
     ;if Re ge (radi) then sflag = 0
     if (Re ge nx-1) or (Re ge ny-1) then sflag = 0 
     
     ;Sersic index
     n = (10D^pars[sorder+4])[0]
     if n lt 0.4 then sflag = 0
     if n ge 8. then sflag = 0
     
     ;Sersic b/a 
     serE = pars[sorder+5]   
     if serE lt 0.2D then sflag = 0
     if serE gt 1D then sflag = 0
     
     ;Sersic PA
     ;serPA = ACOS(pars[sorder+6])*180D/!DPI 
     serPA = (pars[sorder+6])
     if serPa lt -360D then sflag = 0
     if serPA gt 360D then sflag = 0

     ;if (sflag eq 0) then stop
  endif

  if expnn gt 0 then begin
     expX0 = pars[eorder]
     expY0 = pars[eorder+1]
     eflag = 1
     if expx0 lt ((nx/2D) - 20) then eflag = 0
     if expx0 gt ((nx/2D) + 20) then eflag = 0
     
     if expy0 lt ((ny/2D) - 20) then eflag = 0
     if expy0 gt ((ny/2D) + 20) then eflag = 0
     
     ;Central intensity 
     I0 = (10D^(pars[eorder+2]))[0]
     ;if mu0 gt 30 then efalg = 0
     if I0 lt 0. then eflag = 0

     ;Scale height 
     h = (10D^(pars[eorder+3]))[0]
     if h lt 1D then eflag = 0
     ;if h ge (2D*radi) then eflag = 0
     if (h ge nx-1) or (h ge ny-1) then eflag = 0

     ;Exponential b/a 
     expE = pars[eorder+4]   
     if expE lt 0.3D then eflag = 0
     if expE gt 1D then eflag = 0

     ;Exponential PA
     expPA = (pars[eorder+5]) 
     if expPa lt -360D then eflag = 0
     if expPA gt 360D then eflag = 0
     
  endif

  if cennn gt 0 then begin
     
     X0 = pars[cenorder]
     Y0 = pars[cenorder+1]
     cenflag = 1
     if x0 lt ((nx/2D) - 20) then cenflag = 0
     if x0 gt ((nx/2D) + 20) then cenflag = 0
     
     if y0 lt ((ny/2D) - 20) then cenflag = 0
     if y0 gt ((ny/2D) + 20) then cenflag = 0

     ;if (cenflag eq 0) then stop
  endif


  Bndflag = 1
  if (sernn gt 0) and (expnn gt 0) then begin
     
     bn = (1.9992D * n) - 0.3271D
     Lser = (2D * !dpi * Ie * (Re^2D) * n * (EXP(bn)) * gamma(2D *n)) / (bn^(2D * n))
     Lexp = (2D * !dpi * I0 * (h^2D))
     BD = Lser / Lexp
     hRe = Re / h
     BD_1Re = BD / (2D * igamma(2D,hRe))


     ;Bulge-to-disc ratio
     if (BD lt 0.01D) then Bndflag = 0
     if (BD gt 9D) then Bndflag = 0
     
     ;Re to h ratio
     if (Re gt (1.678 * h)) then Bndflag = 0
     if (Re lt (0.05D * h)) then Bndflag = 0
     

     ;Bulge-to-disc ratio(< Re)
     if (BD_1Re lt 0.1) then Bndflag = 0
     

     ;Bulge-to-disc ratio < Re
     ;BD_1Re = BD / (2D * igamma(2D,hRe) * gamma(2D)) 
     ;if (BD_1Re le 1.5D) then flag[0] = 0

     
     ;Double crossing points 
     ;x = [1D:200D]
     ;dexp = (I0*exp(-x/h))
   
     ;invn = 1D / n
     ;bn = (1.9992 * n) - 0.3271     
     ;dser = Ie * EXP(-bn * ( ((x/Re)^(invn)) - 1D) )
     
     ;pars2 = [Ie,Re,n,I0,h]
     ;res = dser - dexp
     ;crX = cross(x,res,func='func1',ffunc='func2',pars=pars2)
     
     ;if (crX[1] lt 7D*Re) and (crX[1] ne 0) then flag[1] = 0
     ;if (TOTAL(crX) eq 0) or (crX[0] gt 7D*Re) then Bndflag = 0
     ;if (TOTAL(crX) eq 0) or (crX[0] eq 0) then flag[2] = 0
     ;if (crX[1] lt 10D*Re) and (crX[0] lt 10D*Re) and (crX[0] ne 0) and (crX[1] ne 0) then stop
     ;if BnDflag eq 0 then stop
     
     ;fexp = (I0*exp(-x/h))
     
     ;invn = 1D / n
     ;bn = (1.9992 * n) - 0.3271
     
     ;fser = Ie * EXP(-bn * ( ((x/Re)^(invn)) - 1D) )
     
     ;if fser[0] le fexp[0] then flag[4] = 0
     
  endif
  return, flag
END

Pro bagal_priors,pars,varmat=varmat
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  common mflags, fsersic,fsersic2,fsersic3,fsersic4,fexp,fbexp,ffer,ffer2
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  skip = i
;  format = ''
;  readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
;           serX0, serY0, Ie, Re, n, serE, serPA,$
;           expX0, expY0, I0, h, expE, expPA,$
;           ser2X0, ser2Y0, Ie2, Re2, n2, ser2E, ser2PA,$
;           ser3X0, ser3Y0, Ie3, Re3, n3, ser3E, ser3PA,$
;           ser4X0, ser4Y0, Ie4, Re4, n4, ser4E, ser4PA,$
;           bexpX0,$
;           FORMAT=format,$
;           SKIPLINE=skip,/SILENT,NUMLINE=1
;  print, format

  strlist = ['NAME','seein','seein2','betaaa',$
             'serX0', 'serY0', 'Ie', 'Re', 'n', 'serE', 'serPA',$
             'expX0', 'expY0', 'I0', 'h', 'expE', 'expPA',$
             'ser2X0', 'ser2Y0', 'Ie2', 'Re2', 'n2', 'ser2E', 'ser2PA',$
             'ser3X0', 'ser3Y0', 'Ie3', 'Re3', 'n3', 'ser3E', 'ser3PA',$
             'ser4X0', 'ser4Y0', 'Ie4', 'Re4', 'n4', 'ser4E', 'ser4PA',$
             'bexpX0','bexpY0','bI0','h1','h2','bexpa','rb','bexpE','bexpPA',$
             'ferX0','ferY0','ferI0','l','fern','ferE','ferPA','ferc',$
             'fer2X0','fer2Y0','fer2I0','l2','fer2n','fer2E','fer2PA','fer2c'$
            ]

  head_length = 1
  file = './inputs/bagal_inputs.txt'
  data_line=fltarr(n_elements(strlist))       ; Creates an empty array of floating point numbers 
  head_line=''                  ; Creates an empty string (for text) 
  header=strarr(head_length)    ; Creates an empty array of strings 

  j = 0
  
  openr,lun,file, /GET_LUN               ; open the file to read from 
  
;If there's a header then read to header string array 

  for k=0,head_length-1 do begin
     readf,lun,head_line
     header(k)=head_line
  endfor 

;Read data 
  repeat begin
     readf,lun,data_line
;Reading only specified line
     if j eq skip then data = data_line 
     j=j+1l
  endrep until eof(lun)                             ;read until end of the file (eof) 
  close,lun

  for jj=1, n_elements(data)-1 do begin
     exec = EXECUTE(strlist[jj] + ' = ' + string(data[jj]))
  endfor
  
  
  
  sorder = where(nlist eq 'serX0', sernn)
  sorder2 = where(nlist eq 'ser2X0', sernn2)
  sorder3 = where(nlist eq 'ser3X0', sernn3)
  sorder4 = where(nlist eq 'ser4X0', sernn4)
  eorder = where(nlist eq 'expX0', expnn)
  beorder = where(nlist eq 'bexpX0', bexpnn)
  forder = where(nlist eq 'ferX0', fernn)
  forder2 = where(nlist eq 'fer2X0', fernn2)
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)

  pars = MAKE_ARRAY(1,/FLOAT)
  varmat = MAKE_ARRAY(1,/FLOAT)
  
;SERSIC PROFILE
  serx0_temp = !Values.F_NAN
  sery0_temp = !Values.F_NAN
  if sernn gt 0 then begin

     
     ;Centre coordinates x0 and y0
     fserx0 = flist[0,sorder]
     fsery0 = flist[0,sorder+1]
     pserx0 = flist[1,sorder]
     psery0 = flist[1,sorder+1]
     if serX0 eq -1 then begin
        serX0 = round(nx/2.D)
        if fserx0 eq 1 then serx0 = round(vprior(pserx0,serX0,a=1.))
     endif
     if serY0 eq -1 then begin 
        serY0 = round(ny/2.D)
        if fsery0 eq 1 then sery0 = round(vprior(psery0,sery0,a=1.))
     endif
     serx0_temp = serx0
     sery0_temp = sery0

     ;The effective intensity 
     fIe = flist[0,sorder+2]
     priIe = flist[1,sorder+2]
     if Ie eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fIe eq 1 then Ie = vprior(priIe,Ie,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Ie = ALOG10(Ie)
     ;Ie = Ie / B ;NORMALISED 

     
     

     ;The effective radius 
     fRe = flist[0,sorder+3]
     priRe = flist[1,sorder+3]
     if Re eq -1 then begin
        Re = Uprior(1, 1., 3.)  ;Re in pixel
        if fRe eq 1 then Re = vprior(priRe,Re,a=0.1)
     endif 
     Re = ALOG10(Re)

     ;The Sersic index
     fn = flist[0,sorder+4]
     prin = flist[1,sorder+4]
     if n eq -1 then begin 
        n = Uprior(1, 1., 6.) 
        if fn eq 1 then n = vprior(prin,n,a=0.1)
     endif
     n = ALOG10(n)
     
     ;The ellipticity
     fserE = flist[0,sorder+5]
     priserE = flist[1,sorder+5]
     if serE eq -1 then begin 
        serE = Uprior(1, 0.7, 1.) 
        if fserE eq 1 then serE = vprior(priserE,serE,a=0.01)
     endif
     
     
     ;The position angle 
     fserPA = flist[0,sorder+6]
     priserPA = flist[1,sorder+6]
     if serPA eq -1 then begin
        serPA = Uprior(1, 0., 90.) 
        if fserPA eq 1 then serPA = vprior(priserPA,serPA,a=0.1)
     endif

     parsser = [serX0,serY0,Ie,Re,n,serE,serPA]
    
     pars = [pars,parsser]
     
     aser = [serx0/10.,sery0/10.,Ie/10.,Re/100.,n/100.,serE/10.,serPA/10.]
     dser = [serx0/50.,sery0/50.,Ie/50.,Re/500.,n/500.,serE/50.,serPA/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D]
     varmat = [varmat,parscale]
  endif

  ;SERSIC 2 PROFILE
  ser2x0_temp = !Values.F_NAN
  ser2y0_temp = !Values.F_NAN
  if sernn2 gt 0 then begin

     
     ;Centre coordinates x0 and y0
     fser2x0 = flist[0,sorder2]
     fser2y0 = flist[0,sorder2+1]
     pser2x0 = flist[1,sorder2]
     pser2y0 = flist[1,sorder2+1]
     if ser2X0 eq -1 then begin
        ser2X0 = round(nx/2.D)
        if fser2x0 eq 1 then ser2x0 = round(vprior(pser2x0,ser2X0,a=1.))
     endif
     if ser2Y0 eq -1 then begin 
        ser2Y0 = round(ny/2.D)
        if fser2y0 eq 1 then ser2y0 = round(vprior(pser2y0,ser2y0,a=1.))
     endif
     ser2x0_temp = ser2x0
     ser2y0_temp = ser2y0

     ;The effective intensity 
     fIe2 = flist[0,sorder2+2]
     priIe2 = flist[1,sorder2+2]
     if Ie2 eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fIe2 eq 1 then Ie2 = vprior(priIe2,Ie2,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Ie2 = ALOG10(Ie2)
     ;Ie = Ie / B ;NORMALISED 

     
     

     ;The effective radius 
     fRe2 = flist[0,sorder2+3]
     priRe2 = flist[1,sorder2+3]
     if Re2 eq -1 then begin
        Re2 = Uprior(1, 1., 3.)  ;Re in pixel
        if fRe2 eq 1 then Re2 = vprior(priRe2,Re2,a=0.1)
     endif 
     Re2 = ALOG10(Re2)

     ;The Sersic index
     fn2 = flist[0,sorder2+4]
     prin2 = flist[1,sorder2+4]
     if n2 eq -1 then begin 
        n2 = Uprior(1, 1., 6.) 
        if fn2 eq 1 then n2 = vprior(prin2,n2,a=0.1)
     endif
     n2 = ALOG10(n2)
     
     ;The ellipticity
     fser2E = flist[0,sorder2+5]
     priser2E = flist[1,sorder2+5]
     if ser2E eq -1 then begin 
        ser2E = Uprior(1, 0.7, 1.) 
        if fser2E eq 1 then ser2E = vprior(priser2E,ser2E,a=0.01)
     endif
     
     
     ;The position angle 
     fser2PA = flist[0,sorder2+6]
     priser2PA = flist[1,sorder2+6]
     if ser2PA eq -1 then begin
        ser2PA = Uprior(1, 0., 90.) 
        if fser2PA eq 1 then ser2PA = vprior(priser2PA,ser2PA,a=0.1)
     endif

     parsser2 = [ser2X0,ser2Y0,Ie2,Re2,n2,ser2E,ser2PA]
    
     pars = [pars,parsser2]
     
     aser2 = [ser2x0/10.,ser2y0/10.,Ie2/10.,Re2/100.,n2/100.,ser2E/10.,ser2PA/10.]
     dser2 = [ser2x0/50.,ser2y0/50.,Ie2/50.,Re2/500.,n2/500.,ser2E/50.,ser2PA/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D]
     varmat = [varmat,parscale]
  endif

  ;SERSIC 3 PROFILE
  ser3x0_temp = !Values.F_NAN
  ser3y0_temp = !Values.F_NAN
  if sernn3 gt 0 then begin

     
     ;Centre coordinates x0 and y0
     fser3x0 = flist[0,sorder3]
     fser3y0 = flist[0,sorder3+1]
     pser3x0 = flist[1,sorder3]
     pser3y0 = flist[1,sorder3+1]
     if ser3X0 eq -1 then begin
        ser3X0 = round(nx/2.D)
        if fser3x0 eq 1 then ser3x0 = round(vprior(pser3x0,ser3X0,a=1.))
     endif
     if ser3Y0 eq -1 then begin 
        ser3Y0 = round(ny/2.D)
        if fser3y0 eq 1 then ser3y0 = round(vprior(pser3y0,ser3y0,a=1.))
     endif
     ser3x0_temp = ser3x0
     ser3y0_temp = ser3y0

     ;The effective intensity 
     fIe3 = flist[0,sorder3+2]
     priIe3 = flist[1,sorder3+2]
     if Ie3 eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fIe3 eq 1 then Ie3 = vprior(priIe3,Ie3,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Ie3 = ALOG10(Ie3)
     ;Ie = Ie / B ;NORMALISED 

     
     

     ;The effective radius 
     fRe3 = flist[0,sorder3+3]
     priRe3 = flist[1,sorder3+3]
     if Re3 eq -1 then begin
        Re3 = Uprior(1, 1., 3.)  ;Re in pixel
        if fRe3 eq 1 then Re3 = vprior(priRe3,Re3,a=0.1)
     endif 
     Re3 = ALOG10(Re3)

     ;The Sersic index
     fn3 = flist[0,sorder3+4]
     prin3 = flist[1,sorder3+4]
     if n3 eq -1 then begin 
        n3 = Uprior(1, 1., 6.) 
        if fn3 eq 1 then n3 = vprior(prin3,n3,a=0.1)
     endif
     n3 = ALOG10(n3)
     
     ;The ellipticity
     fser3E = flist[0,sorder3+5]
     priser3E = flist[1,sorder3+5]
     if ser3E eq -1 then begin 
        ser3E = Uprior(1, 0.7, 1.) 
        if fser3E eq 1 then ser3E = vprior(priser3E,ser3E,a=0.01)
     endif
     
     
     ;The position angle 
     fser3PA = flist[0,sorder3+6]
     priser3PA = flist[1,sorder3+6]
     if ser3PA eq -1 then begin
        ser3PA = Uprior(1, 0., 90.) 
        if fser3PA eq 1 then ser3PA = vprior(priser3PA,ser3PA,a=0.1)
     endif

     parsser3 = [ser3X0,ser3Y0,Ie3,Re3,n3,ser3E,ser3PA]
    
     pars = [pars,parsser3]
     
     aser3 = [ser3x0/10.,ser3y0/10.,Ie3/10.,Re3/100.,n3/100.,ser3E/10.,ser3PA/10.]
     dser3 = [ser3x0/50.,ser3y0/50.,Ie3/50.,Re3/500.,n3/500.,ser3E/50.,ser3PA/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D]
     varmat = [varmat,parscale]
  endif

  ;SERSIC 4 PROFILE
  ser4x0_temp = !Values.F_NAN
  ser4y0_temp = !Values.F_NAN
  if sernn4 gt 0 then begin

     
     ;Centre coordinates x0 and y0
     fser4x0 = flist[0,sorder4]
     fser4y0 = flist[0,sorder4+1]
     pser4x0 = flist[1,sorder4]
     pser4y0 = flist[1,sorder4+1]
     if ser4X0 eq -1 then begin
        ser4X0 = round(nx/2.D)
        if fser4x0 eq 1 then ser4x0 = round(vprior(pser4x0,ser4X0,a=1.))
     endif
     if ser4Y0 eq -1 then begin 
        ser4Y0 = round(ny/2.D)
        if fser4y0 eq 1 then ser4y0 = round(vprior(pser4y0,ser4y0,a=1.))
     endif
     ser4x0_temp = ser4x0
     ser4y0_temp = ser4y0

     ;The effective intensity 
     fIe4 = flist[0,sorder4+2]
     priIe4 = flist[1,sorder4+2]
     if Ie4 eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fIe4 eq 1 then Ie4 = vprior(priIe4,Ie4,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Ie4 = ALOG10(Ie4)
     ;Ie = Ie / B ;NORMALISED 

     
     

     ;The effective radius 
     fRe4 = flist[0,sorder4+3]
     priRe4 = flist[1,sorder4+3]
     if Re4 eq -1 then begin
        Re4 = Uprior(1, 1., 3.)  ;Re in pixel
        if fRe4 eq 1 then Re4 = vprior(priRe4,Re4,a=0.1)
     endif 
     Re4 = ALOG10(Re4)

     ;The Sersic index
     fn4 = flist[0,sorder4+4]
     prin4 = flist[1,sorder4+4]
     if n4 eq -1 then begin 
        n4 = Uprior(1, 1., 6.) 
        if fn4 eq 1 then n4 = vprior(prin4,n4,a=0.1)
     endif
     n4 = ALOG10(n4)
     
     ;The ellipticity
     fser4E = flist[0,sorder4+5]
     priser4E = flist[1,sorder4+5]
     if ser4E eq -1 then begin 
        ser4E = Uprior(1, 0.7, 1.) 
        if fser4E eq 1 then ser4E = vprior(priser4E,ser4E,a=0.01)
     endif
     
     
     ;The position angle 
     fser4PA = flist[0,sorder4+6]
     priser4PA = flist[1,sorder4+6]
     if ser4PA eq -1 then begin
        ser4PA = Uprior(1, 0., 90.) 
        if fser4PA eq 1 then ser4PA = vprior(priser4PA,ser4PA,a=0.1)
     endif

     parsser4 = [ser4X0,ser4Y0,Ie4,Re4,n4,ser4E,ser4PA]
    
     pars = [pars,parsser4]
     
     aser4 = [ser4x0/10.,ser4y0/10.,Ie4/10.,Re4/100.,n4/100.,ser4E/10.,ser4PA/10.]
     dser4 = [ser4x0/50.,ser4y0/50.,Ie4/50.,Re4/500.,n4/500.,ser4E/50.,ser4PA/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D]
     varmat = [varmat,parscale]
  endif

;EXPONENTIAL PROFILE
  expx0_temp = !Values.F_NAN
  expy0_temp = !Values.F_NAN
  
  if expnn gt 0 then begin

     ;Centre coordinates x0 and y0
     fexpx0 = flist[0,eorder]
     fexpy0 = flist[0,eorder+1]
     pexpx0 = flist[1,eorder]
     pexpy0 = flist[1,eorder+1]
     if expX0 eq -1 then begin 
        expX0 = round(nx/2.D)
        if fexpx0 eq 1 then expx0 = round(vprior(pexpx0,expX0),a=1.)
     endif
     if expY0 eq -1 then begin 
        expY0 = round(ny/2.D)
        if fexpy0 eq 1 then expy0 = round(vprior(pexpy0,expy0),a=1.)
     endif
     expx0_temp = expx0
     expy0_temp = expy0


     ;The central intensity 
     fI0 = flist[0,eorder+2]
     priI0 = flist[1,eorder+2]
     if I0 eq -1 then begin 
        ;I0 = 10D^(-1D * 0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fI0 eq 1 then I0 = vprior(priI0,I0,a=1.)
     endif ;else I0 = 10D^(-1D * 0.4 * (mu0 - zcalib - 2.5*alog10((scale^2D))))
     I0 = ALOG10(I0)

     ;The effective radius 
     fh = flist[0,eorder+3]
     prih = flist[1,eorder+3]
     if h eq -1 then begin 
        h = Uprior(1, 2., 5.)   ;Re in pixel
        if fh eq 1 then h = vprior(prih,h,a=1.)
     endif
     h = ALOG10(h)

     ;The ellipticity
     fexpE = flist[0,eorder+5]
     priexpE = flist[1,eorder+5]
     if expE eq -1 then begin 
        expE = Uprior(1, 0.7, 1.) ;Re in pixel
        if fexpE eq 1 then expE = vprior(priexpE,expE,a=1.)
     endif

     ;The position angle 
     fexpPA = flist[0,eorder+6]
     priexpPA = flist[1,eorder+6]
     if expPA eq -1 then begin 
        expPA = Uprior(1, 0., 90.) ;Re in pixel
        if fexpPA eq 1 then expPA = vprior(priexpPA,expPA,a=1.)
     endif
     

     parsexp = [expX0,expY0,I0,h,expE,expPA]
     pars = [pars,parsexp]
     aexp = [expx0/10.,expy0/10.,I0/10.,h/10.,expE/10.,expPA/10.]
     dexp = [expx0/100.,expy0/100.,I0/100.,h/100.,expE/100.,expPA/100.]
     
     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 5D]
     varmat = [varmat,parscale]
     
  endif

  ;BROKEN EXPONENTIAL PROFILE
  bexpx0_temp = !Values.F_NAN
  bexpy0_temp = !Values.F_NAN
  
  if bexpnn gt 0 then begin

     ;Centre coordinates x0 and y0
     fbexpx0 = flist[0,beorder]
     fbexpy0 = flist[0,beorder+1]
     pbexpx0 = flist[1,beorder]
     pbexpy0 = flist[1,beorder+1]
     if bexpX0 eq -1 then begin 
        bexpX0 = round(nx/2.D)
        if fbexpx0 eq 1 then bexpx0 = round(vprior(pbexpx0,bexpX0),a=1.)
     endif
     if bexpY0 eq -1 then begin 
        bexpY0 = round(ny/2.D)
        if fbexpy0 eq 1 then bexpy0 = round(vprior(pbexpy0,bexpy0),a=1.)
     endif
     bexpx0_temp = bexpx0
     bexpy0_temp = bexpy0


     ;The central intensity 
     fbI0 = flist[0,beorder+2]
     pribI0 = flist[1,beorder+2]
     if bI0 eq -1 then begin 
        ;I0 = 10D^(-1D * 0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fbI0 eq 1 then bI0 = vprior(pribI0,bI0,a=1.)
     endif ;else I0 = 10D^(-1D * 0.4 * (mu0 - zcalib - 2.5*alog10((scale^2D))))
     bI0 = ALOG10(bI0)

     ;The scale length 1
     fh1 = flist[0,beorder+3]
     prih1 = flist[1,beorder+3]
     if h1 eq -1 then begin 
        h1 = Uprior(1, 2., 5.)   ;Re in pixel
        if fh1 eq 1 then h1 = vprior(prih1,h1,a=1.)
     endif
     h1 = ALOG10(h1)

     ;The  scale length 2
     fh2 = flist[0,beorder+4]
     prih2 = flist[1,beorder+4]
     if h2 eq -1 then begin 
        h2 = Uprior(1, 2., 5.)   ;Re in pixel
        if fh2 eq 1 then h2 = vprior(prih2,h2,a=1.)
     endif
     h2 = ALOG10(h2)
     
     ;Brake sharpness 
     fbexpa = flist[0,beorder+5]
     pribexpa = flist[1,beorder+5]
     if bexpa eq -1 then begin 
        bexpa = Uprior(1, 0, 1.) 
        if fbexpa eq 1 then bexpa = vprior(pribexpa,bexpa,a=1.)
     endif

      ;The breake radius 
     frb = flist[0,beorder+6]
     prirb = flist[1,beorder+6]
     if rb eq -1 then begin 
        rb = Uprior(1, 2., 5.)   
        if frb eq 1 then rb = vprior(prirb,rb,a=1.)
     endif
     rb = ALOG10(rb)
     
     ;The ellipticity
     fbexpE = flist[0,beorder+7]
     pribexpE = flist[1,beorder+7]
     if bexpE eq -1 then begin 
        bexpE = Uprior(1, 0.7, 1.) 
        if fbexpE eq 1 then bexpE = vprior(pribexpE,bexpE,a=1.)
     endif

     ;The position angle 
     fbexpPA = flist[0,beorder+8]
     pribexpPA = flist[1,beorder+8]
     if bexpPA eq -1 then begin 
        bexpPA = Uprior(1, 0., 90.) 
        if fbexpPA eq 1 then bexpPA = vprior(pribexpPA,bexpPA,a=1.)
     endif
     

     parsexp = [bexpX0,bexpY0,bI0,h1,h2,bexpa,rb,bexpE,bexpPA]
     pars = [pars,parsexp]
     abexp = [bexpx0/10.,bexpy0/10.,bI0/10.,h1/10.,h2/10.,bexpa/10, rb/10, bexpE/10.,bexpPA/10.]
     dbexp = [bexpx0/100.,bexpy0/100.,bI0/100.,h1/100.,h2/100.,bexpa/100,rb/100,bexpE/100.,bexpPA/100.]
     
     parscale = [0.01D,0.01D, 0.01D,0.001D,0.001D, 0.001D, 0.0001D, 0.01D, 5D]
     varmat = [varmat,parscale]
     
  endif
  
  ;FERRORS PROFILE
  ferx0_temp = !Values.F_NAN
  fery0_temp = !Values.F_NAN
  if fernn gt 0 then begin

     
     ;Centre coordinates x0 and y0
     fferx0 = flist[0,forder]
     ffery0 = flist[0,forder+1]
     pferx0 = flist[1,forder]
     pfery0 = flist[1,forder+1]
     if ferX0 eq -1 then begin
        ferX0 = round(nx/2.D)
        if fferx0 eq 1 then ferx0 = round(vprior(pferx0,ferX0,a=1.))
     endif
     if ferY0 eq -1 then begin 
        ferY0 = round(ny/2.D)
        if ffery0 eq 1 then fery0 = round(vprior(pfery0,fery0,a=1.))
     endif
     ferx0_temp = ferx0
     fery0_temp = fery0

     ;The central intensity 
     fferI0 = flist[0,forder+2]
     priferI0 = flist[1,forder+2]
     if ferI0 eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if fferI0 eq 1 then ferI0 = vprior(priferI0,ferI0,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     ferI0 = ALOG10(ferI0)
     ;Ie = Ie / B ;NORMALISED 

     ;Length 
     fl = flist[0,forder+3]
     pril = flist[1,forder+3]
     if l eq -1 then begin
        l = Uprior(1, 1., 3.)
        if fl eq 1 then l = vprior(pril,l,a=0.1)
     endif 
     l = ALOG10(l)

     ;The shape parameter
     ffern = flist[0,forder+4]
     prifern = flist[1,forder+4]
     if fern eq -1 then begin 
        fern = Uprior(1, 1., 6.) 
        if ffern eq 1 then fern = vprior(prifern,fern,a=0.1)
     endif
     fern = ALOG10(fern)
     
     ;The ellipticity
     fferE = flist[0,forder+5]
     priferE = flist[1,forder+5]
     if ferE eq -1 then begin 
        ferE = Uprior(1, 0.7, 1.) 
        if fferE eq 1 then ferE = vprior(priferE,ferE,a=0.01)
     endif
     
     
     ;The position angle 
     fferPA = flist[0,forder+6]
     priferPA = flist[1,forder+6]
     if ferPA eq -1 then begin
        ferPA = Uprior(1, 0., 90.) 
        if fferPA eq 1 then ferPA = vprior(priferPA,ferPA,a=0.1)
     endif

     ;c parameter
     fferc = flist[0,forder+7]
     priferc = flist[1,forder+7]
     if ferc eq -1 then begin 
        ferc = Uprior(1, 0, 1.) 
        if fferc eq 1 then ferc = vprior(priferc,ferc,a=0.01)
     endif
     
     parsser = [ferX0,ferY0,ferI0,l,fern,ferE,ferPA,ferc]
    
     pars = [pars,parsser]
     
     afer2 = [ferx0/10.,fery0/10.,ferI0/10.,l/100.,fern/100.,ferE/10.,ferPA/10.,ferc/10.]
     dfer2 = [ferx0/50.,fery0/50.,ferI0/50.,l/500.,fern/500.,ferE/50.,ferPA/50.,ferc/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D, 0.01D]
     varmat = [varmat,parscale]
  endif

  ;FERRORS 2 PROFILE
  fer2x0_temp = !Values.F_NAN
  fer2y0_temp = !Values.F_NAN
  if fernn2 gt 0 then begin

     
     ;Centre coordinates x0 and y0
     ffer2x0 = flist[0,forder2]
     ffer2y0 = flist[0,forder2+1]
     pfer2x0 = flist[1,forder2]
     pfer2y0 = flist[1,forder2+1]
     if ferX0 eq -1 then begin
        fer2X0 = round(nx/2.D)
        if ffer2x0 eq 1 then fer2x0 = round(vprior(pfer2x0,fer2X0,a=1.))
     endif
     if fer2Y0 eq -1 then begin 
        fer2Y0 = round(ny/2.D)
        if ffer2y0 eq 1 then fer2y0 = round(vprior(pfer2y0,fer2y0,a=1.))
     endif
     fer2x0_temp = fer2x0
     fer2y0_temp = fer2y0

     ;The central intensity 
     ffer2I0 = flist[0,forder2+2]
     prifer2I0 = flist[1,forder2+2]
     if fer2I0 eq -1 then begin
        ;Ie = 10D^(-0.4 * (20D - zcalib - 2.5*alog10((scale^2D))))
        if ffer2I0 eq 1 then fer2I0 = vprior(prifer2I0,fer2I0,a=10.)
     endif ;else Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     fer2I0 = ALOG10(fer2I0)
     ;Ie = Ie / B ;NORMALISED 

     ;Length 
     fl2 = flist[0,forder2+3]
     pril2 = flist[1,forder2+3]
     if l2 eq -1 then begin
        l2 = Uprior(1, 1., 3.)
        if fl eq 1 then l = vprior(pril2,l2,a=0.1)
     endif 
     l2 = ALOG10(l2)

     ;The shape parameter
     ffer2n = flist[0,forder2+4]
     prifer2n = flist[1,forder2+4]
     if fer2n eq -1 then begin 
        fer2n = Uprior(1, 1., 6.) 
        if ffer2n eq 1 then fer2n = vprior(prifer2n,fer2n,a=0.1)
     endif
     fer2n = ALOG10(fer2n)
     
     ;The ellipticity
     ffer2E = flist[0,forder2+5]
     prifer2E = flist[1,forder2+5]
     if fer2E eq -1 then begin 
        fer2E = Uprior(1, 0.7, 1.) 
        if ffer2E eq 1 then fer2E = vprior(prifer2E,fer2E,a=0.01)
     endif
     
     
     ;The position angle 
     ffer2PA = flist[0,forder2+6]
     prifer2PA = flist[1,forder2+6]
     if fer2PA eq -1 then begin
        fer2PA = Uprior(1, 0., 90.) 
        if ffer2PA eq 1 then fer2PA = vprior(prifer2PA,fer2PA,a=0.1)
     endif

     ;c parameter
     ffer2c = flist[0,forder2+7]
     prifer2c = flist[1,forder2+7]
     if fer2c eq -1 then begin 
        fer2c = Uprior(1, 0, 1.) 
        if ffer2c eq 1 then fer2c = vprior(prifer2c,fer2c,a=0.01)
     endif
     
     parsser = [fer2X0,fer2Y0,fer2I0,l2,fer2n,fer2E,fer2PA,fer2c]
    
     pars = [pars,parsser]
     
     afer2 = [fer2x0/10.,fer2y0/10.,fer2I0/10.,l2/100.,fer2n/100.,fer2E/10.,fer2PA/10.,fer2c/10.]
     dfer2 = [fer2x0/50.,fer2y0/50.,fer2I0/50.,l2/500.,fer2n/500.,fer2E/50.,fer2PA/50.,fer2c/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D, 0.01D]
     varmat = [varmat,parscale]
  endif
  
;CENTRES
  if cennn gt 0 then begin 
     
    x0 = mean([serx0_temp,ser2x0_temp,ser3x0_temp,ser4x0_temp,expx0_temp,bexpx0_temp,ferx0_temp,fer2x0_temp],/NAN)
    y0 = mean([sery0_temp,ser2y0_temp,ser3y0_temp,ser4y0_temp,expy0_temp,bexpy0_temp,fery0_temp,fer2y0_temp],/NAN)
       
    
    parscen = [X0,Y0]
    pars = [pars,parscen]

    parscale = [0.01D, 0.01D]
    varmat = [varmat,parscale]
  endif
  
  if skynn gt 0 then begin 
     
    sky = skymean
    parssky = [sky]
    pars = [pars,parssky]

    parscale = [0.001D]
    varmat = [varmat,parscale]
 endif

  
;The Step change array 
  astep = MAKE_ARRAY(1,/FLOAT)
  adel = make_array(1,/FLOAT)
  astep[0] = 0.234
  if fsersic eq 1 then astep = [astep,aser]
  if fexp eq 1 then astep = [astep, aexp]
  if fsersic2 eq 1 then astep = [astep,aser2]
  if fsersic3 eq 1 then astep = [astep,aser3]
  if fsersic4 eq 1 then astep = [astep,aser4]
  if fbexp eq 1 then astep = [astep, abexp]
  if ffer eq 1 then astep = [astep, afer]
  if ffer2 eq 1 then astep = [astep, afer2]
  
  if fsersic eq 1 then adel = [adel,dser]
  if fexp eq 1 then adel = [adel, dexp]
  if fsersic2 eq 1 then adel = [adel,dser2]
  if fsersic3 eq 1 then adel = [adel,dser3]
  if fsersic4 eq 1 then adel = [adel,dser4]
  if fbexp eq 1 then adel = [adel, dbexp]
  if ffer eq 1 then astep = [astep, dfer]
  if ffer2 eq 1 then astep = [astep, dfer2]
  
  astep = transpose([[astep],[adel]])



  REMOVE,0,pars
  REMOVE,0,varmat

  pri_pars = pars


                                ;Check parameters
  if (sernn eq 1) or (expnn eq 1) then begin 
     flagcnt = 0D
     REPEAT BEGIN
        
        if (sernn eq 1) then begin
           if pri_pars[4] gt ALOG10(8D) then begin
              pri_pars[4] = ALOG10(3D)
           endif
           if pri_pars[4] lt AlOG10(0.4) then begin
              pri_pars[4] = ALOG10(3D)
           endif
           if pri_pars[3] lt -1D then begin
              pri_pars[3] = 0D
           endif
           if pri_pars[0] gt (nx/2D)+20D or pri_pars[0] lt (nx/2D)-20D then begin
              pri_pars[0] = nx/2D
           endif
           if pri_pars[1] gt (ny/2D)+20D or pri_pars[1] lt (ny/2D)-20D then begin
              pri_pars[1] = ny/2D
           endif
           
        endif
        
        if (expnn eq 1) then begin
           
           if pri_pars[7] gt (nx/2D)+20D or pri_pars[7] lt (nx/2D)-20D then begin
              pri_pars[7] = nx/2D
           endif
           if pri_pars[8] gt (ny/2D)+20D or pri_pars[8] lt (ny/2D)-20D then begin
              pri_pars[8] = ny/2D
           endif
        endif
        
        if (cennn eq 1) then begin
           if (expnn eq 1) and (sernn eq 1) then begin
              if pri_pars[13] gt (nx/2D)+20D or pri_pars[13] lt (nx/2D)-20D then begin
                 pri_pars[13] = nx/2D
              endif
              if pri_pars[14] gt (ny/2D)+20D or pri_pars[14] lt (ny/2D)-20D then begin
              pri_pars[14] = ny/2D
           endif
           endif else begin
              if pri_pars[7] gt (nx/2D)+20D or pri_pars[7] lt (nx/2D)-20D then begin
                 pri_pars[7] = nx/2D
              endif
              if pri_pars[8] gt (ny/2D)+20D or pri_pars[8] lt (ny/2D)-20D then begin
                 pri_pars[8] = ny/2D
              endif
           endelse
           
        endif
        
        if (expnn eq 1) and (sernn eq 1) then begin
           pflag = modify_initial(pri_pars)
           
           if pflag[0] eq 0 then begin
                                ;Bulge-to-disc ratio < Re
              pri_pars[2] = pri_pars[2] + 0.05
           endif
           if pflag[1] eq 0 then begin
                                ;Double crossing points
              pri_pars[4] = pri_pars[4] - 0.5
              pri_pars[10] = ALOG10((10D^(pri_pars[10])) + 0.5)
           endif
           if pflag[2] eq 0 and pflag[0] eq 1 then begin
                                ;Zero crossing point
              
              if pflag[4] eq 1 then pri_pars[9] = pri_pars[9] + 0.05 else pri_pars[2] = pri_pars[2] + 0.05
              
           endif
                                ;Re/h > 1
           if pflag[3] eq 0 then begin
              pri_pars[3] = ALOG10((10D^(pri_pars[3])) - 0.5)
           endif
           
        endif
        
        
        

        if (expnn eq 0) and (sernn eq 1) then pflag=[1,1,1,1,1]
     ENDREP UNTIL (pflag[0] eq 1) and (pflag[1] eq 1) and (pflag[2] eq 1) and (pflag[3] eq 1) 
  endif
  
  pars = pri_pars


END


;*************************************************************************************************


PRO bagal_setup
  common setup1, tmcmc, bagal_like, nchain, nstep, tstep, pridis, f1dfit
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common exposure, pExpo,Expo,exponames,exponame
  common mflags, fsersic,fsersic2,fsersic3,fsersic4,fexp,fbexp,ffer,ffer2
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma
  
;Read the settings file 
  readcol, './inputs/setup.txt', index, inputs, bbb, format='A,F,A',$
           COMMENT='#',STRINGSKIP='[',/SILENT
 ;GENERAL MCMC INPUTS
  bagal_like = 'likeBay'
  nchain = inputs[0]
  nstep = inputs[1]
  f1dfit = fix(inputs[2])
;IMAGE/INSTRUMENT PARAMETERS
  scale = inputs[3]
  zcalib = inputs[4]
  tpsf = fix(inputs[5])
  betaa = inputs[6]
  fFFT = inputs[7]
  seeing = inputs[8]
  seeing2 = inputs[9]
  sky = inputs[10]
  Pnoise = fix(inputs[11])
  gain = inputs[12]
  RON = inputs[13]
  pExpo = fix(inputs[14])
  Expo = DOUBLE(inputs[15])
  Bnoise = fix(inputs[16])
  skysig = inputs[17]
  tklog = inputs[18]
;COMPONENTS TO FIT
;Serśic profile / Bulge
  fsersic = fix(inputs[19])
;Exponential profile / Disc
  fexp = fix(inputs[20])
;Serśic2 profile / Bulge
  fsersic2 = fix(inputs[21])
;Serśic3 profile / Bulge
  fsersic3 = fix(inputs[22])
;Serśic4 profile / Bulge
  fsersic4 = fix(inputs[23])
;Broken Exponential / Outer disc 
  fbexp = fix(inputs[24])
;Ferrors profile / Main Bar
  ffer = fix(inputs[25])
;Ferrors profile / Secondary Bar
  ffer2 = fix(inputs[26])
;SPECIFICATIONS
  rmin = inputs[27]
  rmax = inputs[28]
  if rmin eq 0 then rmin = fix(rmin)
  sigover = inputs[29]
  fmask = fix(inputs[30])
  maskonly = fix(inputs[31])
  
  dir_out = './outputs/'
  
  readcol, './inputs/setup.txt', index, inputs, b, format='A,A,A', $
           skipline=3, COMMENT='#',STRINGSKIP='(',/SILENT
  
  psfname = strcompress('./inputs/PSF/'+inputs[0]+'.fits',/remove_all)
  exponames = inputs[1] 
  if (exponames ne 'file') and (pexpo eq 2) then exponame = strcompress('./inputs/Exposures/'+exponames+'.fits',/remove_all)

  signame = strcompress('./inputs/sigma_maps/'+inputs[2]+'.fits',/remove_all)
  inname = inputs[3]
  fcen = fix(STRSPLIT(inputs[4],'-',/EXTRACT))
  psky = fix(STRSPLIT(inputs[5],'-',/EXTRACT))
;Sersic
  fserx0 = fix(STRSPLIT(inputs[6],'-',/EXTRACT))
  fsery0 = fix(STRSPLIT(inputs[7],'-',/EXTRACT))
  fIe = fix(STRSPLIT(inputs[8],'-',/EXTRACT))
  fRe = fix(STRSPLIT(inputs[9],'-',/EXTRACT))
  fn = fix(STRSPLIT(inputs[10],'-',/EXTRACT))
  fserE = fix(STRSPLIT(inputs[11],'-',/EXTRACT))
  fserPA = fix(STRSPLIT(inputs[12],'-',/EXTRACT))
 ;Exponential
  fexpx0 = fix(STRSPLIT(inputs[13],'-',/EXTRACT))
  fexpy0 = fix(STRSPLIT(inputs[14],'-',/EXTRACT))
  fI0 = fix(STRSPLIT(inputs[15],'-',/EXTRACT))
  fh = fix(STRSPLIT(inputs[16],'-',/EXTRACT))
  fexpE = fix(STRSPLIT(inputs[17],'-',/EXTRACT))
  fexpPA = fix(STRSPLIT(inputs[18],'-',/EXTRACT))
 ;Sersic2
  fser2x0 = fix(STRSPLIT(inputs[19],'-',/EXTRACT))
  fser2y0 = fix(STRSPLIT(inputs[20],'-',/EXTRACT))
  fIe2 = fix(STRSPLIT(inputs[21],'-',/EXTRACT))
  fRe2 = fix(STRSPLIT(inputs[22],'-',/EXTRACT))
  fn2 = fix(STRSPLIT(inputs[23],'-',/EXTRACT))
  fser2E = fix(STRSPLIT(inputs[24],'-',/EXTRACT))
  fser2PA = fix(STRSPLIT(inputs[25],'-',/EXTRACT))
 ;Sersic3
  fser3x0 = fix(STRSPLIT(inputs[26],'-',/EXTRACT))
  fser3y0 = fix(STRSPLIT(inputs[27],'-',/EXTRACT))
  fIe3 = fix(STRSPLIT(inputs[28],'-',/EXTRACT))
  fRe3 = fix(STRSPLIT(inputs[29],'-',/EXTRACT))
  fn3 = fix(STRSPLIT(inputs[30],'-',/EXTRACT))
  fser3E = fix(STRSPLIT(inputs[31],'-',/EXTRACT))
  fser3PA = fix(STRSPLIT(inputs[32],'-',/EXTRACT))
 ;Sersic4
  fser4x0 = fix(STRSPLIT(inputs[33],'-',/EXTRACT))
  fser4y0 = fix(STRSPLIT(inputs[34],'-',/EXTRACT))
  fIe4 = fix(STRSPLIT(inputs[35],'-',/EXTRACT))
  fRe4 = fix(STRSPLIT(inputs[36],'-',/EXTRACT))
  fn4 = fix(STRSPLIT(inputs[37],'-',/EXTRACT))
  fser4E = fix(STRSPLIT(inputs[38],'-',/EXTRACT))
  fser4PA = fix(STRSPLIT(inputs[39],'-',/EXTRACT))
 ;Broken Exponential
  fbexpx0 = fix(STRSPLIT(inputs[40],'-',/EXTRACT))
  fbexpy0 = fix(STRSPLIT(inputs[41],'-',/EXTRACT))
  fbI0 = fix(STRSPLIT(inputs[42],'-',/EXTRACT))
  fbh1 = fix(STRSPLIT(inputs[43],'-',/EXTRACT))
  fbh2 = fix(STRSPLIT(inputs[44],'-',/EXTRACT))
  fbexpa = fix(STRSPLIT(inputs[45],'-',/EXTRACT))
  frb = fix(STRSPLIT(inputs[46],'-',/EXTRACT))
  fbexpE = fix(STRSPLIT(inputs[47],'-',/EXTRACT))
  fbexpPA = fix(STRSPLIT(inputs[48],'-',/EXTRACT))
  fbexpc = fix(STRSPLIT(inputs[49],'-',/EXTRACT))
 ;Ferrors 
  fferx0 = fix(STRSPLIT(inputs[50],'-',/EXTRACT))
  ffery0 = fix(STRSPLIT(inputs[51],'-',/EXTRACT))
  fferI0 = fix(STRSPLIT(inputs[52],'-',/EXTRACT))
  fl = fix(STRSPLIT(inputs[53],'-',/EXTRACT))
  ffern = fix(STRSPLIT(inputs[54],'-',/EXTRACT))
  fferE = fix(STRSPLIT(inputs[55],'-',/EXTRACT))
  fferPA = fix(STRSPLIT(inputs[56],'-',/EXTRACT))
  fferc = fix(STRSPLIT(inputs[57],'-',/EXTRACT))
 ;Ferrors 2
  ffer2x0 = fix(STRSPLIT(inputs[58],'-',/EXTRACT))
  ffer2y0 = fix(STRSPLIT(inputs[59],'-',/EXTRACT))
  ffer2I0 = fix(STRSPLIT(inputs[60],'-',/EXTRACT))
  fl2 = fix(STRSPLIT(inputs[61],'-',/EXTRACT))
  ffer2n = fix(STRSPLIT(inputs[62],'-',/EXTRACT))
  ffer2E = fix(STRSPLIT(inputs[63],'-',/EXTRACT))
  ffer2PA = fix(STRSPLIT(inputs[64],'-',/EXTRACT))
  ffer2c = fix(STRSPLIT(inputs[65],'-',/EXTRACT))
  
;Specifications 
  
  if fcen[0] eq 1 then fixcen = 1 else fixcen = 0
  if psky[0] eq 1 then fsky = 1 else fsky = 0
 
;Create bagal file from user input file 
  BGmkfile, './inputs/'+inname

;CREATING THE PARAMETER MATRIX AND THE MODEL LIST 
   
  if (fixcen eq 1) then begin
     fx0 = fcen
     fy0 = fcen
     fserx0 = [0,0,0]
     fsery0 = [0,0,0]
     fexpx0 = [0,0,0]
     fexpy0 = [0,0,0]
     fser2x0 = [0,0,0]
     fser2y0 = [0,0,0]
     fser3x0 = [0,0,0]
     fser3y0 = [0,0,0]
     fser4x0 = [0,0,0]
     fser4y0 = [0,0,0]

     fbexpx0 = [0,0,0]
     fbexpy0 = [0,0,0]
     fferx0 = [0,0,0]
     ffery0 = [0,0,0]
     ffer2x0 = [0,0,0]
     ffer2y0 = [0,0,0]
  endif else begin
     fx0 = 0
     fy0 = 0    
  endelse

  pcen = [[fx0],[fy0]]
  pser = [[fserx0], [fsery0], [fIE], [fRe], [fn], [fserE], [fserPA]]
  pexp = [[fexpx0], [fexpy0], [fI0], [fh], [fexpE], [fexpPA]]
  pser2 = [[fser2x0], [fser2y0], [fIE2], [fRe2], [fn2], [fser2E], [fser2PA]]
  pser3 = [[fser3x0], [fser3y0], [fIE3], [fRe3], [fn3], [fser3E], [fser3PA]]
  pser4 = [[fser4x0], [fser4y0], [fIE4], [fRe4], [fn4], [fser4E], [fser4PA]]

  pbexp = [[fbexpx0], [fbexpy0], [fbI0], [fbh1], [fbh2], [fbexpa], [frb], [fbexpE], [fbexpPA]]
  pfer = [[fserx0], [fsery0], [fIE], [fRe], [fn], [fserE], [fserPA]]
  pfer2 = [[fserx0], [fsery0], [fIE], [fRe], [fn], [fserE], [fserPA]]
 

  serlist = ['serX0','serY0','Ie','Re','n','serE','serPA']
  serlist2 = ['ser2X0','ser2Y0','Ie2','Re2','n2','ser2E','ser2PA']
  serlist3 = ['ser3X0','ser3Y0','Ie3','Re3','n3','ser3E','ser3PA']
  serlist4 = ['ser4X0','ser4Y0','Ie4','Re4','n4','ser4E','ser4PA']
  explist = ['expX0','expY0','I0','h','expE','expPA']
  bexplist = ['bexpX0','bexpY0','bI0','bh1','bh2','bexpa','rb','bexpE','bexpPA']
  ferlist = ['ferX0','ferY0','ferI0','l','fern','ferE','ferPA','ferc']
  fer2list = ['fer2X0','fer2Y0','fer2I0','l2','fer2n','fer2E','fer2PA','fer2c']
  cenlist = ['X0','Y0']
  skylist = ['sky']

  ;The flag list
 ;flist = MAKE_ARRAY(1,/STRING)
  flist = MAKE_ARRAY(3,1,/INTEGER)
  if fsersic eq 1 then flist = [[flist],[pser]]
  if fexp eq 1 then flist = [[flist], [pexp]]
  if fsersic2 eq 1 then flist = [[flist],[pser2]]
  if fsersic3 eq 1 then flist = [[flist],[pser3]]
  if fsersic4 eq 1 then flist = [[flist],[pser4]]
  if fbexp eq 1 then flist = [[flist], [pbexp]]
  if ffer eq 1 then flist = [[flist],[pfer]]
  if ffer2 eq 1 then flist = [[flist],[pfer2]]
  if fixcen eq 1 then flist = [[flist],[pcen]]
  if fsky eq 1 then flist = [[flist],[psky]]

  plength = n_elements(flist[0,*])-1
  remov = [[0,0],[1,0],[2,0]]
  REMOVE, remov, flist
  flist = reform(flist,3,plength)
  
  ;Defining groups 
  groups = flist[2,*] 
  flist = flist[0:1,*]
 
;The parameter order
  nlist = MAKE_ARRAY(1,/STRING)
  if fsersic eq 1 then nlist = [nlist, serlist]
  if fexp eq 1 then nlist = [nlist, explist]
  if fsersic2 eq 1 then nlist = [nlist, serlist2]
  if fsersic3 eq 1 then nlist = [nlist, serlist3]
  if fsersic4 eq 1 then nlist = [nlist, serlist4]
  if fbexp eq 1 then nlist = [nlist, bexplist]
  if ffer eq 1 then nlist = [nlist, ferlist]
  if ffer2 eq 1 then nlist = [nlist, fer2list]
  if fixcen eq 1 then nlist = [nlist, cenlist]
  if fsky eq 1 then nlist = [nlist,skylist]
  REMOVE, 0, nlist
 
 ;par_matrix = transpose([[flist],[nlist]])


  
  
  
;The model routine 
  mlist = MAKE_ARRAY(1,/STRING)
  ;if fsersic eq 1 then mlist = [mlist, 'sersic(spars)']
  ;if fexp eq 1 then mlist = [mlist, 'expon(epars)']
  ;REMOVE, 0 , mlist
  ;mlist = strjoin(mlist,'+')
 

;Noise 
  if (Bnoise eq 2) then sigmap = readfits(signame,h,/noscale,/SILENT) else sigmap = ''
  
 return
END

PRO ima_setup
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove 
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  ima = readfits(imaname,h,/noscale,/silent)
  im_dim = size(ima)
  nx = im_dim[1]
  ny = im_dim[2]
  lx = im_dim[4]
  ii_t = lindgen(lx)
  sigma = readfits(whtname,h,/noscale,/SILENT)
  ;sigma1 = SQRT(ima);+ (44.561605210683638D))
  ;sigma = sqrt(ima*10D)
  ;sigma = sqrt(sigma^2D + ima)
  
  ;sigma = sqrt(ima+ (sky*gain) + (ron^2D))  
  
;Check for NANs
  
  nansig = where(FINITE(sigma, /NAN))
  if (min(nansig) ne -1) then sigma[nansig] = 0D
  nanima = where(FINITE(ima, /NAN))
  if (min(nanima) ne -1) then ima[nanima] = 0D
  
  
;Normalisation of the image and sigma
  A = min(ima)
  imanew = ima - A
  B = (max(imanew))
  ima = ima ;/ B
 
  
  sigma = sigma ;/ B
  if (tklog eq 1) then sigma = (1D / ima) * sigma
  if (tklog eq 1) then ima = ALOG(ima)
  ivar = 1D/((sigma)^2D)
 

  if (fmask eq 1) then begin
     mask = readfits(maskname,h,/noscale,/SILENT)
     x0 = nx/2D
     y0 = ny/2D
     
     ;NEED TO MAKE SURE THE OBJECT YOU WANT TO FIT IS CENTRED TO THE IMAGE!!!!!!
     
     mask_value = mask[x0,y0]
     
     if mask_value ne -1 then begin
        if (maskonly eq 0) then begin 
           mask_index = where_xyz((mask ne mask_value and mask gt 0), xind=xind1, yind=yind1)        
           ima[xind1,yind1] = !Values.F_NAN
        
        endif else begin
           mask_index = where_xyz((mask ne mask_value), xind=xind1, yind=yind1)        
           ima[xind1,yind1] = !Values.F_NAN 
        endelse
     endif else begin
         if (maskonly eq 0) then begin 
           
           mask_index = where_xyz((mask eq -1), xind=xind1, yind=yind1)        
           ima[xind1,yind1] = !Values.F_NAN
        
        endif else begin
           mask_index = where_xyz((mask eq -1), xind=xind1, yind=yind1) 
           ima[xind1,yind1] = !Values.F_NAN 
        endelse
        
     endelse
  endif


;Find sky level and the good pixels to fit above 1sigma of sky
  newima = ima
  newima[(nx/2)-5:(nx/2)+5,(ny/2)-5:(ny/2)+5] = !Values.F_NAN
  meanclip, newima,skymean,skysigma, CLIPSIG=2,MAXITER=10
  
  sigabove = (skysigma+skymean)
  
  print, 'IMAGE INFORMATION: '
  print, strcompress('Sky mean: '+string(skymean))
  print, strcompress('Sky sigma: '+string(skysigma))
  
  
  xc = nx/2D
  yc = ny/2D


  radi = 2D
  REPEAT BEGIN
     circ = cir_mask(ima, xc, yc, radi)
     minima = min(ima[circ])
     radi = radi + 1D
  ENDREP UNTIL (minima le sigabove) or (radi eq ((nx-1D)/2D)) OR (radi eq ((ny-1D)/2D))
  minx = xc - radi
  maxx = xc + radi
  miny = yc - radi
  maxy = yc + radi
  subima = ima[minx:maxx,miny:maxy]

  abovesky = where(ima ge sigabove)
  Nabov = DOUBLE(n_elements(ima[abovesky]))
  abovbox = where(subima ge sigabove)
  Nbox = DOUBLE(n_elements(subima[abovbox]))
  frac = Nbox/Nabov
  newradi = radi
  if frac le 0.99 then begin
     REPEAT BEGIN
        circ = cir_mask(ima, xc, yc, newradi)
        minx = xc - newradi
        maxx = xc + newradi
        miny = yc - newradi
        maxy = yc + newradi
        subima = ima[minx:maxx,miny:maxy]
        abovbox = where(subima ge sigabove)
        Nbox = double(n_elements(subima[abovbox]))
        frac = Nbox/Nabov
        newradi = newradi + 1D
     ENDREP UNTIL (newradi eq ((nx/2D)-1)) OR (frac ge 0.9) OR (newradi eq ((ny/2D)-1))
  endif
  
  if (rmin ne 0) and (newradi lt rmin) and (rmin lt ((nx/2D)-1)) and (rmin lt ((ny/2D)-1)) then begin
     newradi = rmin
     minx = xc - newradi
     maxx = xc + newradi
     miny = yc - newradi
     maxy = yc + newradi 
  endif
  
  
  if (rmax ne 0) and (newradi gt rmax) and (rmax lt ((nx/2D)-1)) and (rmax lt ((ny/2D)-1)) then begin
     newradi = rmax
     minx = xc - newradi
     maxx = xc + newradi
     miny = yc - newradi
     maxy = yc + newradi 
  endif

  print, strcompress('Fitting radius: '+string(newradi))
  if (radi lt rmin) then radi = rmin
  print, strcompress('Maximum object radius: '+string(radi))
  print,''

  imasubb = [fix(minx),fix(maxx),fix(miny),fix(maxy)]
  
  
  ;PSF 
  case tpsf of
     0: psf = psf_construct(ndimension=31,func=tpsf,fwhm=seeing)    
     1: psf = psf_construct(ndimension=31,func=tpsf,fwhm=seeing,betaa=betaa) 
     2: psf = psf_construct(ndimension=31,func=tpsf,fwhm=[seeing,seeing2]) 
     3: psf = readfits(psfname,h,/noscale,/silent)
  ENDCASE
  
  return
END

;********************************************************************************
;PHOTMETRIC MODELS + MODEL IMAGE GENERATION

FUNCTION model_array,x0,y0,PA,q,x,y
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi  
  rx = (double(ii_t mod nx) - (x0)) ; x (micropixel) coordinate
  ry = (double(ii_t / nx) - y0)     ; y (micropixel) coordinate
  
  x = reform(rx,nx,ny)
  y = reform(ry,nx,ny)

  rxrot = (-rx * sin(pa)) + (ry * cos(pa))
  ryrot = (-rx * cos(pa)) - (ry * sin(pa))
  
  r = sqrt(rxrot^2D + (ryrot/q)^2D)
  r = reform(r,nx,ny)
  
  return, r
END

FUNCTION remodel,imasub,x0,y0,PA,q,x,y,xsize
;COORDINATES OF MATRIX TO RESAMPLE 
  
  newx = x[imasub[0]:imasub[1],imasub[2]:imasub[3]]
  newy = y[imasub[0]:imasub[1],imasub[2]:imasub[3]]
 
  xsize = size(newx)
  ysize = size(newy)

  xbig = REBIN((newx-0.5), (xsize[1]*10), (xsize[2]*10));,/INTERP)
  ybig = REBIN((newy-0.5), (ysize[1]*10), (ysize[2]*10));,/INTERP)
 
  
  xbig = xbig + 0.05
  ybig = ybig + 0.05
  
 
  rxbig = (-xbig * sin(pa)) + (ybig * cos(pa))
  rybig = (-xbig * cos(pa)) - (ybig * sin(pa))
  
  newr = sqrt(rxbig^2D + (rybig/q)^2D)

  
  return, newr
END

FUNCTION model_array_gen,x0,y0,PA,q,c,x,y
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi  
  rx = (double(ii_t mod nx) - (x0)) ; x (micropixel) coordinate
  ry = (double(ii_t / nx) - y0)     ; y (micropixel) coordinate
  
  x = reform(rx,nx,ny)
  y = reform(ry,nx,ny)

  rxrot = (-rx * sin(pa)) + (ry * cos(pa))
  ryrot = (-rx * cos(pa)) - (ry * sin(pa))

  ellExp = c + 2D
  invEllExp = 1D / ellExp
  
  r = (rxrot^ellExp + (ryrot/q)^ellExp) ^ invEllExp
  r = reform(r,nx,ny)
  
  return, r
END

FUNCTION remodel_gen,imasub,x0,y0,PA,q,c,x,y,xsize
;COORDINATES OF MATRIX TO RESAMPLE 
  
  newx = x[imasub[0]:imasub[1],imasub[2]:imasub[3]]
  newy = y[imasub[0]:imasub[1],imasub[2]:imasub[3]]
 
  xsize = size(newx)
  ysize = size(newy)

  xbig = REBIN((newx-0.5), (xsize[1]*10), (xsize[2]*10));,/INTERP)
  ybig = REBIN((newy-0.5), (ysize[1]*10), (ysize[2]*10));,/INTERP)
 
  
  xbig = xbig + 0.05
  ybig = ybig + 0.05
  
  ellExp = c + 2D
  invEllExp = 1D / ellExp
  
  rxbig = (-xbig * sin(pa)) + (ybig * cos(pa))
  rybig = (-xbig * cos(pa)) - (ybig * sin(pa))
  
  newr = (rxbig^ellExp + (rybig/q)^ellExp) ^ invEllExp

  
  return, newr
END


FUNCTION bagal_model,pars
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove 
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma
  
  Is=0 & Iexp=0 & Is2=0 & Is3=0 & Is4=0 & Ibexp=0 & Ifer=0 & Ifer2=0
  
  
  
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  sorder2 = where(nlist eq 'ser2X0', sernn2)
  sorder3 = where(nlist eq 'ser3X0', sernn3)
  sorder4 = where(nlist eq 'ser4X0', sernn4)
  beorder = where(nlist eq 'bexpX0', bexpnn)
  forder = where(nlist eq 'ferX0', fernn)
  forder2 = where(nlist eq 'fer2X0', fernn2)
  
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)
  
  excstr = make_array(1,/STRING)
  xc=[] & yc=[]
  
  ;T = SYSTIME(1)
  if sernn gt 0 then begin


     sorder = fix(sorder[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        serxc = (pars[cenorder])
        seryc = (pars[cenorder+1])
     endif else begin
        serxc = (pars[sorder])
        seryc = (pars[sorder+1])
     endelse
     xc = [xc,serxc] & yc = [yc, seryc]
     
     Ie = 10D^(pars[sorder+2])
     ;Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Re = 10D^(pars[sorder+3])
     n = 10D^(pars[sorder+4])
     serell = pars[sorder+5]
     ;serpa =  ACOS(pars[sorder+6])*180D/!DPI
     serpa =  (pars[sorder+6])*!DPI/180D
      
;DEFINE THE RADIAL POSITION VECTORS 
     rs = model_array(serxc,seryc,serpa,serell,serx,sery)
     ;writefits, dir_out+'test_r.fits', rs
;CREATE THE MODEL   
     invn = 1D / n
     bn = (1.9992 * n) - 0.3271
     Is = Ie * EXP(-bn * ( ((rs/Re)^(invn)) - 1D) )
     
     excstr = [excstr,'Is']
  endif

 ;Sersic 2
  if sernn2 gt 0 then begin


     sorder2 = fix(sorder2[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        ser2xc = (pars[cenorder])
        ser2yc = (pars[cenorder+1])
     endif else begin
        ser2xc = (pars[sorder2])
        ser2yc = (pars[sorder2+1])
     endelse
     xc = [xc,ser2xc] & yc = [yc, ser2yc]
     
     Ie2 = 10D^(pars[sorder2+2])
     ;Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Re2 = 10D^(pars[sorder2+3])
     n2 = 10D^(pars[sorder2+4])
     ser2ell = pars[sorder2+5]
     ;serpa =  ACOS(pars[sorder+6])*180D/!DPI
     ser2pa =  (pars[sorder2+6])*!DPI/180D
      
;DEFINE THE RADIAL POSITION VECTORS 
     rs2 = model_array(ser2xc,ser2yc,ser2pa,ser2ell,ser2x,ser2y)
     ;writefits, dir_out+'test_r.fits', rs
;CREATE THE MODEL   
     invn2 = 1D / n2
     bn2 = (1.9992 * n2) - 0.3271
     Is2 = Ie2 * EXP(-bn2 * ( ((rs2/Re2)^(invn2)) - 1D) )
     
     excstr = [excstr,'Is2']
  endif

  if sernn3 gt 0 then begin


     sorder3 = fix(sorder3[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        ser3xc = (pars[cenorder])
        ser3yc = (pars[cenorder+1])
     endif else begin
        ser3xc = (pars[sorder3])
        ser3yc = (pars[sorder3+1])
     endelse
     xc = [xc,ser3xc] & yc = [yc, ser3yc]
     
     Ie3 = 10D^(pars[sorder3+2])
     ;Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Re3 = 10D^(pars[sorder3+3])
     n3 = 10D^(pars[sorder3+4])
     ser3ell = pars[sorder3+5]
     ;serpa =  ACOS(pars[sorder+6])*180D/!DPI
     ser3pa =  (pars[sorder3+6])*!DPI/180D
      
;DEFINE THE RADIAL POSITION VECTORS 
     rs3 = model_array(ser3xc,ser3yc,ser3pa,ser3ell,ser3x,ser3y)
     ;writefits, dir_out+'test_r.fits', rs
;CREATE THE MODEL   
     invn3= 1D / n3
     bn3 = (1.9992 * n3) - 0.3271
     Is3 = Ie3 * EXP(-bn3 * ( ((rs3/Re3)^(invn3)) - 1D) )
     
     excstr = [excstr,'Is3']
  endif

  if sernn4 gt 0 then begin


     sorder4 = fix(sorder4[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        ser4xc = (pars[cenorder])
        ser4yc = (pars[cenorder+1])
     endif else begin
        ser4xc = (pars[sorder4])
        ser4yc = (pars[sorder4+1])
     endelse
     xc = [xc,ser4xc] & yc = [yc, ser4yc]
     
     Ie4 = 10D^(pars[sorder4+2])
     ;Ie = 10D^(-0.4 * (mue - zcalib - 2.5*alog10((scale^2D))))
     Re4 = 10D^(pars[sorder4+3])
     n4 = 10D^(pars[sorder4+4])
     ser4ell = pars[sorder4+5]
     ;serpa =  ACOS(pars[sorder+6])*180D/!DPI
     ser4pa =  (pars[sorder4+6])*!DPI/180D
      
;DEFINE THE RADIAL POSITION VECTORS 
     rs4 = model_array(ser4xc,ser4yc,ser4pa,ser4ell,ser4x,ser4y)
     ;writefits, dir_out+'test_r.fits', rs
;CREATE THE MODEL   
     invn4 = 1D / n4
     bn4 = (1.9992 * n4) - 0.3271
     Is4 = Ie4 * EXP(-bn4 * ( ((rs4/Re4)^(invn4)) - 1D) )
     
     excstr = [excstr,'Is4']
  endif
  
  if expnn gt 0 then begin
     eorder = fix(eorder[0])
   
     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        expxc = (pars[cenorder])
        expyc = (pars[cenorder+1])
     endif else begin
        expxc = (pars[eorder])
        expyc = (pars[eorder+1])
     endelse
     xc = [xc,expxc] & yc = [yc, expyc]
     
     I0 = 10D^(pars[eorder+2])
     h = 10D^(pars[eorder+3])
     expell = pars[eorder+4]
     exppa = (pars[eorder+5])*!DPI/180D
     
;DEFINE THE RADIAL POSITION VECTORS 
     
     rexp = model_array(expxc,expyc,exppa,expell,expx,expy)
     
     Iexp = I0*EXP(-rexp/h)
     
     excstr = [excstr,'Iexp']
  endif

  if bexpnn gt 0 then begin
     beorder = fix(beorder[0])
   
     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        bexpxc = (pars[cenorder])
        bexpyc = (pars[cenorder+1])
     endif else begin
        bexpxc = (pars[beorder])
        bexpyc = (pars[beorder+1])
     endelse
     xc = [xc,bexpxc] & yc = [yc, bexpyc]
     
     bI0 = 10D^(pars[beorder+2])
     h1 = 10D^(pars[beorder+3])
     h2 = 10D^(pars[beorder+4])
     bexpa = pars[beorder+5]
     rb = 10D^(pars[beorder+6])
     bexpell = pars[beorder+7]
     bexppa = (pars[beorder+8])*!DPI/180D
     
;DEFINE THE RADIAL POSITION VECTORS 
     
     rbexp = model_array(bexpxc,bexpyc,bexppa,bexpell,bexpx,bexpy)
     
     dh = (1D/h1) - (1D/h2)
     inva = 1D/bexpa
     s = (1D + EXP(-bexpa*rb) ) ^ (-inva * dh)
     z1 = bexpa * (rbexp - rb)
     z = (1 + EXP(z1)) ^ (inva * dh)

     
     Ibexp = s * bi0 * EXP(-1D * (rbexp/h1)) * z
     
     excstr = [excstr,'Ibexp']
  endif 

  if fernn gt 0 then begin


     forder = fix(forder[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        ferxc = (pars[cenorder])
        feryc = (pars[cenorder+1])
     endif else begin
        ferxc = (pars[forder])
        feryc = (pars[forder+1])
     endelse
     xc = [xc,ferxc] & yc = [yc, feryc]
     
     ferI0 = 10D^(pars[forder+2])
     l = 10D^(pars[forder+3])
     fern = 10D^(pars[forder+4])
     ferell = pars[forder+5]
     ferpa =  (pars[forder+6])*!DPI/180D
     ferc = pars[forder+7]
     
;DEFINE THE RADIAL POSITION VECTORS 
     rf = model_array_gen(ferxc,feryc,ferpa,ferell,ferc,ferx,fery)
     
;CREATE THE MODEL   
     ferindex = where_xyz(rf gt l, Xind=xfer, Yind=yfer)
     rf[xfer, yfer] = 0
     Ifer = feri0 * (rf - (rf / l)^2) ^ (fern + 0.5)
     
     excstr = [excstr,'Ifer']
  endif

  if fernn2 gt 0 then begin


     forder2 = fix(forder2[0])
     

     if cennn gt 0 then begin
        cenorder = fix(cenorder[0])
        fer2xc = (pars[cenorder])
        fer2yc = (pars[cenorder+1])
     endif else begin
        fer2xc = (pars[forder2])
        fer2yc = (pars[forder2+1])
     endelse
     xc = [xc,fer2xc] & yc = [yc, fer2yc]
     
     fer2I0 = 10D^(pars[forder2+2])
     l2 = 10D^(pars[forder2+3])
     fer2n = 10D^(pars[forder2+4])
     fer2ell = pars[forder2+5]
     fer2pa =  (pars[forder2+6])*!DPI/180D
     fer2c = pars[forder2+7]
     
;DEFINE THE RADIAL POSITION VECTORS 
     rf2 = model_array_gen(fer2xc,fer2yc,fer2pa,fer2ell,fer2c,fer2x,fer2y)
     
;CREATE THE MODEL   
     ferindex = where_xyz(rf2 gt l2, Xind=xfer2, Yind=yfer2)
     rf2[xfer2, yfer2] = 0
     Ifer2 = fer2i0 * (rf2 - (rf2 / l2)^2) ^ (fer2n + 0.5)
     
     excstr = [excstr,'Ifer2']
  endif

  
  
  Itot = Is + Iexp + Is2 + Is3 + Is4 + Ibexp + Ifer + Ifer2

  
  

  
  
  
  
;RESAMPLE FOR THE STEEP GRADIENTS     
  
;GRADIENT BETWEEN EACH PIXEL
  
 gradima = gradient(Itot)
  gradgrad = gradient(gradima)
  grad_shot = sqrt(gradima)/sqrt(Itot)
  gradthres = 0.5 ;Threshold set where the gradient noise is half that of the shot noise or greater 
  gradindex = where_xyz(grad_shot ge gradthres, Xind=xind, Yind=yind)
  
  ;stop
  xc = median(xc)
  yc = median(yc)
  gradradi = 2D
  circ = cir_mask(grad_shot, xc, yc, radi)
  
  
  REPEAT BEGIN
     circ = cir_mask(grad_shot, xc, yc, gradradi)
     mingrad = min(grad_shot[circ])
     gradradi = gradradi + 1D
  ENDREP UNTIL (mingrad le gradthres) or (gradradi ge ((nx-2D)/2D)) or (gradradi ge ((ny-2D)/2D))
  xmin = Round(xc - gradradi)
  xmax = Round(xc + gradradi)
  ymin = Round(yc - gradradi)
  ymax = Round(yc + gradradi)
  
                                ;Due to how REBIN works it duplicates
                                ;the last n/m as it interpolates and
                                ;not extrapolates so we add on 2
                                ;pixels to re-sample to accomodate for this
  xsubmin = xmin-2
  xsubmax = xmax+2 
  ysubmin = ymin-2
  ysubmax = ymax+2
  
   
  xsubdis = xsubmax - xsubmin
  ysubdis = ysubmax - ysubmin

  if xsubdis ge 20 then begin
      xmin = Round(xc - 10)
      xmax = Round(xc + 10)
      xsubmin = xmin-2
      xsubmax = xmax+2 
     
  endif
  if ysubdis ge 20 then begin
      ymin = Round(yc - 10)
      ymax = Round(yc + 10)
      ysubmin = ymin-2
      ysubmax = ymax+2
  endif
  


  xsubdis = xsubmax - xsubmin
  ysubdis = ysubmax - ysubmin
  Imasub = [xsubmin,xsubmax,ysubmin,ysubmax]
  

  
  ;Remodel/resample
  Isnew=0 & Iexpnew=0 & Is2new=0 & Is3new=0 & Is4new=0 & Ibexpnew=0 & Ifernew=0 & Ifer2new=0
  if sernn gt 0 then begin
     sernewr = remodel(Imasub,serxc,seryc,serpa,serell,serx,sery,serxsize)
     Isnew = Ie * EXP( -bn * (((sernewr/Re)^invn)-1D))
  endif
   if sernn2 gt 0 then begin
     sernewr2 = remodel(Imasub,ser2xc,ser2yc,ser2pa,ser2ell,ser2x,ser2y,ser2xsize)
     Is2new = Ie2 * EXP( -bn2 * (((sernewr2/Re2)^invn2)-1D))
  endif
    if sernn3 gt 0 then begin
     sernewr3 = remodel(Imasub,ser3xc,ser3yc,ser3pa,ser3ell,ser3x,ser3y,ser3xsize)
     Is3new = Ie3 * EXP( -bn3 * (((sernewr3/Re3)^invn3)-1D))
  endif
     if sernn4 gt 0 then begin
     sernewr4 = remodel(Imasub,ser4xc,ser4yc,ser4pa,ser4ell,ser4x,ser4y,ser4xsize)
     Is4new = Ie4 * EXP( -bn4 * (((sernewr4/Re4)^invn4)-1D))
  endif
  if expnn gt 0 then begin
     expnewr = remodel(Imasub,expxc,expyc,exppa,expell,expx,expy,expxsize)
     Iexpnew = I0 * EXP(-expnewr/h)
  endif
  if bexpnn gt 0 then begin
     bexpnewr = remodel(Imasub,bexpxc,bexpyc,bexppa,bexpell,bexpx,bexpy,bexpxsize)
     z1 = bexpa * (bexpnewr - rb)
     z = (1 + EXP(z1)) ^ (inva * dh)
     Ibexpnew = s * bi0 * EXP(-1D * (bexpnewr/h1)) * z
  endif
  if fernn gt 0 then begin
     fernewr = remodel_gen(Imasub,ferxc,feryc,ferpa,ferell,ferc,ferx,fery,ferxsize)
     ferindex = where_xyz(fernewr gt l, Xind=xfer, Yind=yfer)
     fernewr[xfer, yfer] = 0
     Ifernew = feri0 * (fernewr - (fernewr / l)^2) ^ (fern + 0.5)
  endif
  if fernn2 gt 0 then begin
     fer2newr = remodel_gen(Imasub,fer2xc,fer2yc,fer2pa,fer2ell,fer2c,fer2x,fer2y,fer2xsize)
     ferindex = where_xyz(fer2newr gt l2, Xind=xfer2, Yind=yfer2)
     fer2newr[xfer2, yfer2] = 0
     Ifer2new = fer2i0 * (fer2newr - (fer2newr / l2)^2) ^ (fer2n + 0.5)
  endif
  Itotnew = Isnew + Iexpnew + Is2new + Is3new + Is4new + Ibexpnew + Ifernew
  Itotnew = rebin(Itotnew, xsubdis+1, ysubdis+1)

  
  Itot[xmin:xmax,ymin:ymax] = Itotnew[2:n_elements(itotnew[*,0])-3,2:n_elements(itotnew[0,*])-3]

  ;PREFORM THE CONVOLUTION 
  
  xmin = fix(imasubb[0])-10
  if xmin lt 1 then xmin = fix(imasubb[0])
  if xmin lt 1 then xmin = 0

  xmax =  fix(imasubb[1])+10
  if xmax ge nx then xmax = fix(imasubb[1])
  if xmax ge nx then xmax = fix(nx-1)

  ymin =  fix(imasubb[2])-10
  if ymin lt 1 then ymin = fix(imasubb[2])
  if ymin lt 1 then ymin = 0

  ymax =  fix(imasubb[3])+10
  if ymax ge ny then ymax = fix(imasubb[3])
  if ymax ge ny then ymax = fix(ny-1)
  
  ;T = SYSTIME(1)
  if skynn gt 0 then begin
     skyval = (pars[skyorder])[0]
     Itot = Itot + skyval
  endif ;else Itot = Itot + skymean

  
  Itotnew = Itot[xmin:xmax,ymin:ymax]
  Itotal = Itot
  
 
   if fFFT eq 0 then begin
     ;Itotal1 = convol_fft(Itotnew,psf) ;,KERNEL_FFT=newpsf)
     Itotal1 = convolve_j(Itotnew,psf)

  endif else begin
     Itotal1 = convol(Itotnew,psf,/edge_wrap,/NAN)
  endelse
  
  ;PRINT, 'MCMC took ', (SYSTIME(1) - T), ' Seconds to complete'
  
  Itotal[xmin:xmax,ymin:ymax] = Itotal1
  ;Itotal = Itotal1
  
  if (Pnoise ne 0) or (Bnoise ne 0) then begin
     sigstring = make_array(1,/STRING)
     if (Pnoise eq 1) then sigstring = [sigstring,'(sqrt(Itotal*gain))']
     if Bnoise eq 1 then sigstring = [sigstring,'(make_array(nx,ny,/DOUBLE,value=skysig))']
     if Bnoise eq 2 then sigstring = [sigstring,'(sigmap)']
     REMOVE,0,sigstring
     sigstring = strcompress('sigtot = sqrt('+strjoin(sigstring,' + ')+')')
     ;exc = execute(sigstring)
     noise = Nrandom(SIG=(sigtot))
     Itotal = Itotal + noise
  endif
 
  return, Itotal
END

;********************************************************************************


;ALL LIKELIHOODS ARE => LOG[L(D|M)]
FUNCTION like,pars
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  
  ;T = SYSTIME(1) 
  model= bagal_model(pars)
  
  ;PRINT, SYSTIME(1) - T, ' Seconds'
  
  chisq = total((ima - model)*(ima - model)*ivar,/double)
  like = (-0.5*chisq)
  stop
  return, like
END

;For large data sets Poisson distribution can be estimated with
;gaussian 
Function likePos,pars
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  xmin = fix(imasubb[0])
  xmax =  fix(imasubb[1])
  ymin =  fix(imasubb[2])
  ymax =  fix(imasubb[3])
  
  
  model = bagal_model(pars)
  ;Normalise the model
  model = (model) / B
  
  ;To account for -ve fluxes we give 
  ;fixed sky level to raise all pixel
  ;model = model 
  ;newima = ima

  ;Poisson distribution
  ratio = ( ( model[xmin:xmax,ymin:ymax] ) / ( ima[xmin:xmax,ymin:ymax] ) )
  like = TOTAL(-model[xmin:xmax,ymin:ymax] + (ima[xmin:xmax,ymin:ymax]*(ALOG(ratio) + 1D)),/DOUBLE,/NAN)
  Pri = call_function('priors',Pars)
  likeli = like + Pri
  ;print, like, Pri
  if (FINITE(likeli, /NAN)) then message, 'Like is NaN'  
 
 
  return, likeli
END

Function likeBay,pars
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  common chistuff, isig, numbit, Nof, Pri
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common exposure, pExpo,Expo,exponames,exponame
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma
   
  xmin = fix(imasubb[0])
  xmax =  fix(imasubb[1])
  ymin =  fix(imasubb[2])
  ymax =  fix(imasubb[3])
  
  npars = n_elements(pars)
  
  ;T = SYSTIME(1)  
  model = bagal_model(pars)

  npars = n_elements(pars)
 
  ;Number of data points 
  N = (N_ELEMENTS(IMA[xmin:xmax,ymin:ymax]))
  NoF = N
  
  newivar = ivar;*6.25D ;newsigma^(-2D)
 
  if (tklog eq 1) then model = ALOG(model) 

  ;Chi squared
  chisq = TOTAL((ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * (ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * newivar[xmin:xmax,ymin:ymax],/DOUBLE,/NAN)
 
  ;Partition function
  isig = total(ALOG(sqrt((1D /newivar[xmin:xmax,ymin:ymax]))),/DOUBLE,/NAN)
  numbit = (N/2D) * ALOG(2D * !DPI)
  
  ;Likelihood 
  like = ((-0.5D * chisq) - isig - numbit ) 

 
  ;ires = ima - model
  ;writefits, './outputs/test_image.fits', ima
  ;writefits, './outputs/test_res.fits', Ires
  ;writefits, './outputs/test.fits', model
  ;stop
  ;Priors
  ;priors = call_function('priors',Pars)
  likee = like ;+ priors

  return, likee
END


FUNCTION likeBagal,pars
  common setup1, tmcmc, bagal_like,nchain, nstep, tstep, pridis, f1dfit
  
  like = call_function(bagal_like,pars)

  return, like
END

;********************************************************************************
;Perform 1D fit to gain initial parameter values

function buldis, xfit, valores, serfit = serfit, expfit = expfit
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common d1pdf, psf_1d
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  a=valores
  
  I_0_osb=fltarr(100D);n_elements(xfit))
  I_0_osd=fltarr(100D);n_elements(xfit))
  
;CREATE THE BULGE MODEL
  if (sernn eq 1) then begin
     rb=xfit
     a(2) = abs(a(2))
     a(1) = abs(a(1))
     invn = 1D /a(2)
     
     ;bn=0.868D *a(2)-0.142D
     bn = (1.9992 * a(2)) - 0.3271
     xre = (rb/a(1))^invn
     
     ;I_0_osb = a(0)*10D ^(-bn*(xre-1D))
     I_0_osb = a(0) * EXP(-bn * (xre - 1D))
     
     
     serfit = I_0_osb
  endif
  
;CREATE THE DISK MODELS
  if (expnn eq 1) then begin
     rd=xfit
     I_0_osd=a(3)*exp(-rd/a(4))
     expfit = I_0_osd
  endif
  
;TOTAL PROFILE
  I_0_os=I_0_osb+I_0_osd
  

;Convolve with psf
  ;newpsf_1d = REFORM(psf_1d,100)
  
  Ndim = N_elements(I_0_os)
  Npsf = N_elements(psf_1d)
  delPsf = Npsf - Ndim
  
  
  if (delPsf ge 0) then newI_0_os  = convol(I_0_os ,psf_1d[(delPsf/2D):Npsf-(delPsf/2D)-1],/EDGE_truncate,CENTER=1) else newI_0_os  = convol(I_0_os ,psf_1d,/EDGE_truncate,CENTER=1)
  
  ;newI_0_os = REBIN((newI_0_os), (Ndim))
  
  
  return,newI_0_os

end


pro bagal_1d,name,pars,varmat=varmat
  common setup1, tmcmc, bagal_like,nchain, nstep, tstep, pridis, f1dfit
  common setup2, scale, zcalib,tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  common exposure, pExpo,Expo,exponames,exponame
  common priorss,  pri_pars
  common d1pdf, psf_1d
  print,'#################################################################'
  print,' OBTAINING INITIAL ESTIMATES ....    '




;########################################################################
;DEFINE SOME VARIABLES
  chisquare=dblarr(30)
  parameters=dblarr(5,30)
  sigma1=dblarr(5,15)
  valores=dblarr(5)
  valpadisk=dblarr(1)
  valbadisk=dblarr(1)
  valpabulge=dblarr(1)
  valbabulge=dblarr(1)
  radmaxaut=dblarr(1)

;#######################################################################
;READ THE 1D PROFILE FROM ELLIPSE
  modelname=strcompress('./sample/'+name+'.tab',/REMOVE_ALL)
  
  table_ext,modelname,'sma,intens,ellip,pa',x,y,ell,pa

  ;normfactor[ttt]=mean(y[0:fix(seeing[ttt])]/10D)

;#######################################################################


  numlin=n_elements(x)

  nfijo=0
  
  difprof=y-(sigabove)
  difprof1=where(difprof lt 0.,hj)
  if (hj eq 0) then begin
     nfijo=numlin
  endif else begin
     nfijo=min(difprof1)
  endelse
  

  numlin=nfijo-1
  
  x=x[0:numlin-1]
  y=y[0:numlin-1]
  ell=ell[0:numlin-1]
  pa=pa[0:numlin-1]
  
;######################################################################
;PERFORM A LINEAR FIT FOR BOTH DISK AND BULGE TO OBTAIN A FIRST GUESS
;FOR THE 1D NON-LINEAR FIT

  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  cenorder = where(nlist eq 'X0', cennn)
  
  disk=dblarr(2)
  ndisc = 1
  disc = make_array(2,ndisc,/DOUBLE)
  
  if expnn gt 0 then begin
     
     ind=where(x[0:numlin-1] gt 0.5D *x[numlin-1])

     xdisk=x(ind)
     ydisk=y(ind)
     
     ydisk1=alog(ydisk)
     disk=poly_fit(xdisk,ydisk1,1,/DOUBLE)
     
     disk[0]=exp(disk[0])
     disk[1]=-1D /(disk[1])
     
     if (disk[0] le 0.) then disk[0]=-disk[0]
     if (disk[1] le 0.) then disk[1]=-disk[1]
     
     
     ydisk2=disk[0]*exp(-(x/disk[1]))
     
    
     
  endif else begin
     ydisk2 = 0D 
  endelse
 
;FROM THIS FIT WE OBTAIN I0->disk[0] and h->disk[1]
;NOW WE SUBTRACT THE DISK AND CALCULATE THE BULGE PARAMETERS

  bulge=dblarr(2)
  if sernn gt 0 then begin
     ;n=findgen(15)*0.5D +0.5D
     n=findgen(30)*0.25D +0.25D
     
     for qqq=0,n_elements(n)-1 do begin
        
        ybulge=y-ydisk2
        minrad=0.4D
        ind1=where(ybulge gt 0 and x lt minrad*x[numlin-1], nk)
        xbulge=x^(1D /n[qqq])
        
        
        lp=0.1
        
        while (nk le 2) do begin
           if (expnn eq 1) then begin
              disk[0]=10D^(ALOG10(disk[0])-lp)
              ;disk[1]=disk[1]+(lp*0.01D)
              disk[1]=disk[1]+(lp)
           endif
           ydisk2=disk[0]*exp(-(x/disk[1]))
           ybulge=y-ydisk2
           ind1=where(ybulge gt 0 and x lt minrad*x[numlin-1],nk)
           
           if (nk le 3 and expnn eq 0) then minrad=minrad+0.5
        
           
        endwhile
        
        ybulge1=alog10(ybulge(ind1))
        xbulge=xbulge(ind1)
        
        bulge=poly_fit(xbulge,ybulge1,1,/DOUBLE)
        
        bulge[1]=-bulge[1]
        bulge[0]=10D ^(bulge[0])
        bn=0.868D *n[qqq]-0.142D
        bulge[1]=((bn/bulge[1])^n[qqq])
        ybulge2=bulge[0]*10D ^(-bn*(x/bulge[1])^(1D /n[qqq]))
        
        
        ytotal=ybulge2+ydisk2
        
        
;THE FIT VALUES FOR THE BULGE ARE are bulge[1]->Re and bulge[0]->Ie

        ap=dblarr(5)
        ap(0)=((bulge[0])/(10D ^(bn)))
        ap(1)=bulge[1]
        ap(3)=disk[0]
        ap(4)=disk[1]
        ap(2)=n[qqq]
        
        chisquare[qqq]=total(abs(y-ytotal))
        
        parameters[*,qqq]=ap
        
     endfor

     
     ;SELECT THE BEST CHIÂ² OF THE N-GRID FITTING
     bestchi=min(chisquare,lpp,/NAN) ;select the best chi^2 of the n-grid fitting
  endif else begin
     bulge[0]=0D
     bulge[1]=0D
     n=0D
     
     ap=dblarr(5)
     ap(0)=bulge[0]
     ap(1)=bulge[1]
     ap(3)=disk[0]
     ap(4)=disk[1]
     ap(2)=0D

  endelse
  
  
  if sernn eq 0 then begin
     valores[*]=ap
  endif else begin
     
     if (n_elements(lpp) gt 1) then lpp=7
     valores[*]=parameters[*,lpp]
     
     ss1=where(finite(valores[*]) eq 0,nk)
     if (nk gt 0) then valores[ss1]=1D
  endelse
;##################################################################
; IF THE PREVIOUS FIT DO NOT HAVE SENSE; WE MODIFY THE INITIAL CONDITIONS
  if sernn eq 0 then begin
     valores[0]=valores[0]
     valores[1]=valores[1]
  endif else begin
     if (valores[0] le 0.01 or valores[0] gt 60000) then begin
        see_lim=fix(2D)
        see_lim2=where(x lt see_lim)
        see_lim3=max(see_lim2)
        valores[0]=y[see_lim3]
     endif
     
     if (valores[1] lt 1.5D or valores[1] gt x[numlin-1]) then begin
        see_lim=1.5D 
        valores[1]=see_lim
     endif
  endelse
  
  if (expnn eq 0) then begin
     valores[3]=valores[3]
  endif else begin
     if (valores[3] le 0.01 or valores[3] gt 60000) then begin
        see_lim=fix(valores[4])
        see_lim2=where(x lt see_lim)
        see_lim3=max(see_lim2)
        valores[3]=y[see_lim3]
     endif
  endelse

;####################################################################
;####################################################################
;Now, once we have a first estimation of the parameters, we perform a
;non-linear fit with this values as initial conditions, in order to
;obtain  guesses more accurate
  indfit=where(x lt x[numlin-1] and x gt 2D,jk)
;###################################################################
;#Esto lo aÃ±adimos el 15/02/06 para que no se corte en los ajustes
;#de galaxias muy debiles.
  if (indfit(0) eq -1) then indfit=where(x lt x[numlin-1],jk)
;#################################################################
  xfit=x(indfit)
  yfit=y(indfit)
  erro=dblarr(jk)+1.
  
  parinfo=replicate({fixed:0,limited:[0,0],limits:[0.D,0]},5)
  
  if (sernn eq 1) then begin
     parinfo[0].limited = [1,1]
     parinfo[0].limits=[0.01D,60000D]
     
     parinfo[1].limited = [1,1]
     parinfo[1].limits=[1D,x[numlin-1]]
     
     
     parinfo[2].limited = [1,1]
     parinfo[2].limits=[0.5D,8D]
     
  endif else begin
     parinfo[0].fixed = 1
     parinfo[1].fixed = 1
     parinfo[2].fixed = 1
  endelse
  
  if (expnn eq 1) then begin
     parinfo[3].fixed = 1
     parinfo[4].fixed = 1
  endif else begin
     parinfo[3].limited = [1,1]
     parinfo[3].limits=[0.01D,60000D]
     
     parinfo[4].limited = [1,1]
     parinfo[4].limits=[1D,x[numlin-1]]
  endelse
  
 
  
;OBTAIN THE 1-D PSF FUNCTION 
  
  tmp=where(x le valores[1],ntmp)
  indbulge=[max(tmp),max(tmp)+1]
  if (ntmp ge numlin or tmp(0) lt 0.) then indbulge=[1D,$
                                                     fix(2D)+1]
  
  babulgeprom=(ell[indbulge[0]]+ell[indbulge[1]])/2D
  pabulgeprom=(pa[indbulge[0]]+pa[indbulge[1]])/2D
  
  
  valbabulge=1D -babulgeprom
  valpabulge=pabulgeprom
  if valpabulge le 0. then valpabulge=valpabulge+180D
  if (valbabulge lt 0.2D) then valbabulge = 0.5D
  if (valbabulge ge 1D) then valbabulge = 0.99D
  
;PSF 
  dim_psf = size(psf)
  updim = ceil(dim_psf[1]/2D)
  downdim = floor(dim_psf[1]/2D)
  if (updim eq downdim) then updim = updim+1

  case tpsf of
     0: begin
        psf_1d =  ROT(psf, valpabulge, /INTERP) 
     end
     1: begin
        psf_1d =  ROT(psf, valpabulge, /INTERP)
     end
     2: begin
        psf_1d =  ROT(psf, valpabulge, /INTERP)
     end
     3: begin
        psf_1d = ROT(psf, valpabulge, /INTERP)
        psf_1d = psf_1d[*,downdim:updim]
        psf_1d = (psf_1d[*,0]+psf_1d[*,1])/2D
     end
  ENDCASE
  
  
  
;####################################################################
;           PERFORM THE NON-LINEAR 1D BULGE-DISK DECOMPOSITION
;####################################################################
  newguess=mpfitfun('buldis',xfit,yfit,erro,valores[*],ftol=1.e-7,$
                    parinfo=parinfo,perror=perror,bestnorm=bestnorm,dof=dof,/QUIET)
;#####################################################################
  
 ; valores[*]=newguess
  

  skip = 1+i
  readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
           serX0, serY0, Ie, Re, n, serE, serPA,$
           expX0, expY0, I0, h, expE, expPA,$
           FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',$
           SKIPLINE=skip,/SILENT,NUMLINE=1


  pars = MAKE_ARRAY(1,/FLOAT)
  varmat = MAKE_ARRAY(1,/FLOAT)
  serx0_temp = !Values.F_NAN
  sery0_temp = !Values.F_NAN
  expx0_temp = !Values.F_NAN
  expy0_temp = !Values.F_NAN

  aa = -23.726900
  kk = 0.095942
  airmass = 1.1977736
  
  zcal= -(aa+kk*airmass)+2.5*alog10(53.907456*0.396^2.)
  
  if (sernn eq 1) then begin
     
     mue = -2.5*ALOG10(valores[0])+zcal
     Ie = ALOG10(valores[0])
     re = ALOG10(valores[1])
     n = valores[2]
     if (10D^Re) lt 1D then re = ALOG10(1.1D)
     if (n lt 0.5) then n = 0.5D
     if (n gt 12D) then n = 12D
     
     tmp=where(x le valores[1],ntmp)
     indbulge=[max(tmp),max(tmp)+1]
     if (ntmp ge numlin or tmp(0) lt 0.) then indbulge=[1D,$
                                                        fix(2D)+1]
     
     babulgeprom=(ell[indbulge[0]]+ell[indbulge[1]])/2D
     pabulgeprom=(pa[indbulge[0]]+pa[indbulge[1]])/2D
     
     
     valbabulge=1D -babulgeprom
     valpabulge=pabulgeprom
     if valpabulge le 0. then valpabulge=valpabulge+180D
     if (valbabulge lt 0.2D) then valbabulge = 0.5D
     if (valbabulge ge 1D) then valbabulge = 0.99D
          
     serx0_temp = serx0
     sery0_temp = sery0
     
     pars = [pars,serx0,sery0,Ie,re,n,valbabulge,cos(valpabulge*!DPI/180D)]
     parscale = [0.1D,0.1D, double(0.01), 0.05D, 0.5D, 0.01D, 0.01D]
     varmat = [varmat,parscale]
  endif 

  if (expnn eq 1) then begin

     mu0 = -2.5*ALOG10(valores[3])+zcal
     I0 = ALOG10(valores[3])
     h = ALOG10(valores[4])

     if (10D^h) lt 1D then h = ALOG10(2D)

     if (numlin le 11) then begin
        valbadisk=1D -median(ell(numlin-4:numlin-1))
        valpadisk=median(pa(numlin-4:numlin-1))
     endif else begin
        valbadisk=1D -median(ell(numlin-10:numlin-1))
        valpadisk=median(pa(numlin-10:numlin-1))
        if valpadisk le 0D then valpadisk=valpadisk+180D
     endelse

     if (valbadisk lt 0.2D) then valbadisk = 0.5D
     if (valbadisk ge 1D) then valbadisk = 0.99D
     
     expx0_temp = expx0
     expy0_temp = expy0

     pars = [pars,expx0,expy0,I0,h,valbadisk,cos(valpadisk*!DPI/180D)]

     parscale = [0.5D,0.5D, double(0.01), 0.05D, 0.01D, 0.01D]
     
     varmat = [varmat,parscale]
  endif

  ;CENTRES
  if cennn gt 0 then begin 
     
    x0 = mean([serx0_temp,expx0_temp],/NAN)
    y0 = mean([sery0_temp,expy0_temp],/NAN)
       
    
    parscen = [X0,Y0]
    pars = [pars,parscen]

    parscale = [0.1D,0.1D]
    varmat = [varmat,parscale]
 endif

  REMOVE,0,pars
  REMOVE,0,varmat

  pri_pars = pars
 
   

  ;Check parameters 
  flagcnt = 0D
  REPEAT BEGIN
     
     if (sernn eq 1) then begin
        if pri_pars[4] gt 5 then begin
           pri_pars[4] = 3D
        endif
         if pri_pars[4] lt 0.4 then begin
           pri_pars[4] = 3D
        endif
        if pri_pars[0] gt (nx/2D)+5D or pri_pars[0] lt (nx/2D)-5D then begin
           pri_pars[0] = nx/2D
        endif
        if pri_pars[1] gt (ny/2D)+5D or pri_pars[1] lt (ny/2D)-5D then begin
           pri_pars[1] = ny/2D
        endif
        
     endif

     if (expnn eq 1) then begin
      
        if pri_pars[7] gt (nx/2D)+5D or pri_pars[7] lt (nx/2D)-5D then begin
           pri_pars[7] = nx/2D
        endif
        if pri_pars[8] gt (ny/2D)+5D or pri_pars[8] lt (ny/2D)-5D then begin
           pri_pars[8] = ny/2D
        endif

        if (10D^pri_pars[10] ge nx-1) or (10D^pri_pars[10] ge ny-1) then pri_pars[10] = ALOG10(3D*(10D^pri_pars[3]))
     endif
     
     if (cennn eq 1) and (expnn eq 1) then begin
        
        if pri_pars[13] gt (nx/2D)+5D or pri_pars[13] lt (nx/2D)-5D then begin
           pri_pars[13] = nx/2D
        endif
        if pri_pars[14] gt (ny/2D)+5D or pri_pars[14] lt (ny/2D)-5D then begin
           pri_pars[14] = ny/2D
        endif
     endif

      if (cennn eq 1) and (expnn eq 0) then begin
        
        if pri_pars[7] gt (nx/2D)+5D or pri_pars[7] lt (nx/2D)-5D then begin
           pri_pars[7] = nx/2D
        endif
        if pri_pars[8] gt (ny/2D)+5D or pri_pars[8] lt (ny/2D)-5D then begin
           pri_pars[8] = ny/2D
        endif
     endif

     if (expnn eq 1) and (sernn eq 1) then begin
        pflag = modify_initial(pri_pars)
        
        if pflag[0] eq 0 then begin
                                ;Bulge-to-disc ratio < Re
           pri_pars[2] = pri_pars[2] + 0.05
        endif
        if pflag[1] eq 0 then begin
                                ;Double crossing points
           pri_pars[4] = pri_pars[4] - 0.5
           pri_pars[10] = ALOG10((10D^(pri_pars[10])) + 0.5)
        endif
        if pflag[2] eq 0 and pflag[0] eq 1 then begin
                                ;Zero crossing point
           pri_pars[9] = pri_pars[9] + 0.05
        endif
                                ;Re/h > 1
        if pflag[3] eq 0 then begin
           pri_pars[3] = ALOG10((10D^(pri_pars[3])) - 0.5)
        endif
         
        
         
     endif else pflag = [1,1,1,1,1]
     
     
  ENDREP UNTIL (pflag[0] eq 1) and (pflag[1] eq 1) and (pflag[2] eq 1) and (pflag[3] eq 1)
  ;stop
  pars = pri_pars
  if (expnn eq 1) and (sernn eq 1) then valores = [10D^pars[2],10D^pars[3],pars[4],10D^pars[9],10D^pars[10]]
  if (expnn eq 0) and (sernn eq 1) then valores = [10D^pars[2],10D^pars[3],pars[4]]
  
  ytotal = buldis( xfit, valores, serfit = serfit, expfit = expfit)
  chisquare=total(((yfit-ytotal)^2D)/(sqrt(ytotal/1400D)^2D))
  

  print,''
  print,' 1D fit finished ...'
  print, ' Initial parameters are ... '
  print, ''
  if (sernn eq 1) then begin
     print, 'Sersic component:'
     print, strcompress('  Effective surface brightness / Intensity: '+string(mue)+' / '+string(10D^pri_pars[2]))
     print, strcompress('  Effective radius: '+string(10D^pri_pars[3]))
     print, strcompress('  Sersic index: '+string(pri_pars[4]))
     print, strcompress('  b/a: '+string(valbabulge))
     print, strcompress('  Position angle: '+string(valpabulge))
     print,''
  endif
    if (expnn eq 1) then begin
     print, 'Exponential component:'
     print, strcompress('  Central surface brightness / Intensity: '+string(mu0)+' / '+string(10D^pri_pars[9]))
     print, strcompress('  Scale height: '+string(10D^pri_pars[10]))
     print, strcompress('  b/a: '+string(valbadisk))
     print, strcompress('  Position angle: '+string(valpadisk))
  endif
  print,'#################################################################'
  print, ''

  
end



;********************************************************************************
PRO PHI
  common setup1, tmcmc, bagal_like,nchain, nstep, tstep, pridis, f1dfit
  common setup2, scale, zcalib,tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  common exposure, pExpo,Expo,exponames,exponame
  
  spawn,'clear'
  seed= -1L

  bagal_setup
  
  
  readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
           FORMAT='A,F,F,F',$
           SKIPLINE=1,/SILENT
  readcol, './inputs/weights/weights.txt',whtnames,$
           FORMAT='A',/SILENT

  if (pExpo eq 2) and (exponames eq 'file') then readcol,'./inputs/Exposures/Exposures.txt',exponames,$
           FORMAT='A',/SILENT 

  if (fmask eq 1) then  readcol, './inputs/Masks/Masks.txt',Masknames,$
           FORMAT='A',/SILENT
     
     
  
  
  nima = n_elements(NAME)
  
  parameters = fltarr(nima,nchain,n_elements(nlist),long(nstep))
  outpars = fltarr(nima,(n_elements(nlist)*3))
  likelihood = fltarr(nima,nchain,long(nstep))
  information = make_array(nima,(6+(nchain*2)),/DOUBLE)

 

print,'#######################'
print, ''
print,'        ooooo   '
print,"        `888'   "
print,'    .d888888888b.  '
print,'  d88P"  888  "Y88b   '
print,'  888    888    888   '
print,'  888    888    888   '
print,'  888    888    888   '
print,'  888    888    888   '
print,'  Y88b.  888  .d88P   '
print,'    "Y888888888P"    '
print,'         888     '
print,'        o888o    '
print, ''           
print,'#######################'
print, ''


  for i=0l,nima-1 do begin
  
     imaname = strcompress('./sample/'+NAME[i]+'.fits',/remove_all)
     whtname = strcompress('./inputs/weights/'+whtnames[i]+'.fits',/remove_all)
     if (fmask eq 1) then maskname = strcompress('./inputs/Masks/'+masknames[i]+'.fits',/remove_all) else maskname = ''
     if (pExpo eq 2) then begin
        if (n_elements(exponames) gt 1)then begin 
           exponame = strcompress('./inputs/Exposures/'+Exponames[i]+'.fits',/remove_all)
        endif else exponame = strcompress('./inputs/Exposures/'+Exponames+'.fits',/remove_all)
        expo = readfits(exponame,h,/noscale,/SILENT)
     endif

     print, ' '
     print, 'Fitting : '+NAME[i] 
     print, ' '

     if betaaa[i] NE -1 and tpsf eq 1 then betaa = betaaa[i]
     if seein[i] NE -1 then seeing = seein[i]
     if seein2[i] NE -1 and tpsf eq 2 then seeing2 = seein2[i]
     
     
     ima_setup
     pars = []
     
     ;cgDisplay , 1000.,1000.,/ASPECT
     if (f1dfit eq 1) then bagal_1d,NAME[i],pars,varmat=varmat else bagal_priors,pars,varmat=varmat

     
     T = SYSTIME(1) 
     PHI_AM,like, nstep, pars,'likeBay',scale=varmat, $
                      adapt=nstep, accrate=0.234, $
                      nstart=1, nrep=nchain,info=info
     PRINT,''
     PRINT, 'MCMC took ', (SYSTIME(1) - T)/60D, ' Minutes to complete'
     
     ;parameters[i,*,*,*] = pars
     ;likelihood[i,*,*] = like
     sorder = where(nlist eq 'serX0', sernn)
     sorder2 = where(nlist eq 'ser2X0', sernn2)
     sorder3 = where(nlist eq 'ser3X0', sernn3)
     sorder4 = where(nlist eq 'ser4X0', sernn4)
     eorder = where(nlist eq 'expX0', expnn)
     beorder = where(nlist eq 'bexpX0', bexpnn)
     forder = where(nlist eq 'ferX0', fernn)
     forder2 = where(nlist eq 'fer2X0', fernn2)
     
     bit = '_'
     if (sernn eq 1) then bit = bit+'s'
     if (sernn2 eq 1) then bit = bit+'s'
     if (sernn3 eq 1) then bit = bit+'s'
     if (sernn4 eq 1) then bit = bit+'s'
     if (expnn eq 1) then bit = bit+'e'
     if (bexpnn eq 1) then bit = bit+'b'
     if (fernn eq 1) then bit = bit+'f'
     if (fernn2 eq 1) then bit = bit+'f'
     
     if (sernn eq 1) then pars[sorder+4,*] = 10D^pars[sorder+4,*]
     if (sernn2 eq 1) then pars[sorder2+4,*] = 10D^pars[sorder2+4,*]
     if (sernn3 eq 1) then pars[sorder3+4,*] = 10D^pars[sorder3+4,*]
     if (sernn4 eq 1) then pars[sorder4+4,*] = 10D^pars[sorder4+4,*]


     
     
     output_data = make_array(n_elements(pars[*,0])+1,n_elements(pars[0,*]),/DOUBLE)
     for ii=0, n_elements(pars[0,*])-1 do output_data[*,ii] = [pars[*,ii],like[ii]]
     WRITEFITS, dir_out+NAME[i]+bit+'_dataout.fits', output_data
     
     information[i,*] = info
     xmin = fix(imasubb[0])
     xmax =  fix(imasubb[1])
     ymin =  fix(imasubb[2])
     ymax =  fix(imasubb[3])
     like = like                ;/(n_elements(ima[xmin:xmax,ymin:ymax]))

     pars_fit = nlist[where(flist[0,*] eq 1)]
     loglist = ['Ie', 'Re',$
                'I0', 'h',$
                'Ie2', 'Re2', $
                'Ie3', 'Re3', $
                'Ie4', 'Re4', $
                'bI0','bh1','bh2','rb',$
                'ferI0','l','fern',$
                'fer2I0','l2','fer2n'$
               ]
     
     for ii=0, n_elements(nlist)-1 do begin
        if flist[0,ii] eq 1 then begin
           info = postinfo(pars[ii,*])
           if total(strmatch(loglist, nlist[ii], /FOLD_CASE)) ge 1 then info = postinfo(10D^(pars[ii,*])) else info = postinfo(pars[ii,*])
           med = info.med
           low = info.siglow
           hi = info.sighigh
           error = ((med - low) + (hi - med)) / 2D
           print, $
              strcompress(nlist[ii] + ' = ' + string(med) + '+/-' + string(error))
        endif
     endfor

     
     stop

  endfor
  
  ;parameters[*,2,*] = (parameters[*,2,*] * B)
  
  openw,2,dir_out+'BAGAL_info.dat',width=500
     
  for i=0, nima-1 do begin   
     printf,2,name[i],transpose(information[i,*])
  endfor 
  close,2

stop
END



function plotcolor,x

xnew = round(255*(x - min(x))/(max(x)-min(x)))

col = byte(xnew)
return,col

end

function arr_round, array
  array1 = Float(Round(array*100)/100.)  
  return, array1
end

pro bagal_plots,pars,like,NAME,newparas,newlists,nrep
  common setup2, scale, zcalib,tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common mflags, fsersic,fsersic2,fsersic3,fsersic4,fexp,fbexp,ffer,ffer2
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf

  
  ;Filtering the values according to the median of the likelihood
  likemed = median(like)
  indlike = where(like gt likemed)
  ;like = like[indlike]
  ;pars = pars[*,indlike]
  
  nn = long(n_elements(pars[0,*]))
  discard = 0;nn/4l
  like = like[discard:*]
  
  newparas = MAKE_ARRAY(1,/FLOAT)
  newlists = MAKE_ARRAY(1,/STRING)
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)
  
  excstr = make_array(1,/STRING)
  if sernn gt 0 then begin

     sorder = fix(sorder[0])
     serx0 = reform(pars[sorder,discard:*])
     sxinfo = postinfo(serx0)
     sxmed = sxinfo.med
     sxlow = sxinfo.siglow
     sxhi = sxinfo.sighigh
     
     sery0 = reform(pars[sorder+1,discard:*])
     syinfo = postinfo(sery0)
     symed = syinfo.med
     sylow = syinfo.siglow
     syhi = syinfo.sighigh

     Ie = (reform(pars[sorder+2,discard:*]))
     nlist[sorder+2] = 'log(Ie)'
     mue = ie;-2.5*alog10(Ie) + zcalib + 2.5*alog10((scale^2D))
     ieinfo = postinfo(mue)
     muemed = ieinfo.med
     muelow = ieinfo.siglow
     muehi = ieinfo.sighigh

     Re = 10D^(reform(pars[sorder+3,discard:*]))
     reinfo = postinfo(Re)
     remed = reinfo.med
     relow = reinfo.siglow
     rehi = reinfo.sighigh
     
     n = (reform(pars[sorder+4,discard:*]))
     ninfo = postinfo(n)
     nmed = ninfo.med
     nlow = ninfo.siglow
     nhi = ninfo.sighigh
     
     serE = reform(pars[sorder+5,discard:*])
     bainfo = postinfo(serE)
     sbamed = bainfo.med
     sbalow = bainfo.siglow
     sbahi = bainfo.sighigh

     serpa =  reform(pars[sorder+6,discard:*]);acos(reform(pars[sorder+6,discard:*]))*180D/!DPI
     ;serpa =  acos(reform(pars[sorder+6,discard:*]))*180D/!DPI
     spainfo = postinfo(serpa)
     spamed = spainfo.med
     spalow = spainfo.siglow
     spahi = spainfo.sighigh
     newparas = [newparas,sxmed, sxlow, sxhi,symed, sylow, syhi,muemed,muelow,muehi,remed,relow,rehi,nmed,nlow,nhi,sbamed,sbalow,sbahi,spamed,spalow,spahi]
     newlists = [newlists,'sxmed', 'sxlow', 'sxhi','symed', 'sylow', 'syhi','Iemed','Ielow','Iehi','remed','relow','rehigh','nmed','nlow','nhi','sbamed','sbalow','sbahi','spamed','spalow','spahi']
  endif

  if expnn gt 0 then begin
   
     eorder = fix(eorder[0])
     expx0 = reform(pars[eorder,discard:*])
     exinfo = postinfo(expx0)
     exmed = exinfo.med
     exlow = exinfo.siglow
     exhi = exinfo.sighigh
     
     expy0 = reform(pars[eorder+1,discard:*])
     eyinfo = postinfo(expy0)
     eymed = eyinfo.med
     eylow = eyinfo.siglow
     eyhi = eyinfo.sighigh

     I0 = (reform(pars[eorder+2,discard:*]))
     nlist[eorder+2] = 'log(I0)'
     mu0 = i0;-2.5*alog10(I0) + zcalib + 2.5*alog10((scale^2D))
     i0info = postinfo(mu0)
     mu0med = i0info.med
     mu0low = i0info.siglow
     mu0hi = i0info.sighigh

     h = 10D^(reform(pars[eorder+3,discard:*]))
     hinfo = postinfo(h)
     hmed = hinfo.med
     hlow = hinfo.siglow
     hhi = hinfo.sighigh
     
     expE = reform(pars[eorder+4,discard:*])
     bainfo = postinfo(expE)
     ebamed = bainfo.med
     ebalow = bainfo.siglow
     ebahi = bainfo.sighigh

     ;exppa = acos(reform(pars[eorder+5,discard:*]))*180D/!DPI
     exppa = reform(pars[eorder+5,discard:*]); acos(reform(pars[eorder+5,discard:*]))*180D/!DPI
     epainfo = postinfo(exppa)
     epamed = epainfo.med
     epalow = epainfo.siglow
     epahi = epainfo.sighigh
     newparas = [newparas,exmed, exlow, exhi,eymed, eylow, eyhi,mu0med,mu0low,mu0hi,hmed,hlow,hhi,ebamed,ebalow,ebahi,epamed,epalow,epahi]
     newlists = [newlists,'exmed', 'exlow', 'exhi','eymed', 'eylow', 'eyhi','I0med','I0low','I0hi','hmed','hlow','hhigh','ebamed','ebalow','ebahi','epamed','epalow','epahi']
  endif

  if cennn gt 0 then begin
     cenorder = fix(cenorder[0])
     x0 = reform(pars[cenorder,discard:*])
     xinfo = postinfo(x0)
     xmed = xinfo.med
     xlow = xinfo.siglow
     xhi = xinfo.sighigh
     
     y0 = reform(pars[cenorder+1,discard:*])
     yinfo = postinfo(y0)
     ymed = yinfo.med
     ylow = yinfo.siglow
     yhi = yinfo.sighigh
     
     newparas = [newparas,xmed, xlow, xhi,ymed, ylow, yhi]
     newlists = [newlists,'xmed', 'xlow', 'xhi','ymed', 'ylow', 'yhi']
 
  endif
  if skynn gt 0 then begin
     skyorder = fix(skyorder[0])
     sky = reform(pars[skyorder,discard:*])
     skyinfo = postinfo(sky)
     skymed = skyinfo.med
     skylow = skyinfo.siglow
     skyhi = skyinfo.sighigh
     
     
     newparas = [newparas,skymed, skylow, skyhi]
     newlists = [newlists,'skymed', 'skylow', 'skyhi']
 
  endif
  REMOVE, 0 , newparas
  remove, 0, newlists
  ;if expnn gt 0 then begin
    ; eorder = fix(eorder[0])
    ; expxc = pars[eorder]
    ; expyc = pars[eorder+1]
    ; I0 = pars[eorder+2]
    ; h = pars[eorder+3]
    ; expell = pars[eorder+4]
   ;  exppa = (!DPI*pars[eorder+5])/180D
     
;DEFINE THE RADIAL POSITION VECTORS 
     
  ;   rexp = model_array(expxc,expyc,exppa,expell,expx,expy)
     
  ;   Iexp = I0*EXP(-rexp/h)
  ;   excstr = [excstr,'Iexp']
  ;endif 
  ;REMOVE, 0 , excstr
  ;excstr = strjoin(excstr,'+')
  ;exc = execute('Itot ='+excstr)
  

;PLOT THE POSTERIOR DISTRIBUTIONS
  
  if sernn gt 0 then begin
     sxrange = strcompress('['+strjoin(string(arr_round([min(serx0)-(sxlow/50D), max(serx0)+(sxhi/50D)])),',')+']')
     syrange = strcompress('['+strjoin(string(arr_round([min(sery0)-(sylow/50D), max(sery0)+(syhi/50D)])),',')+']')
     nrange = strcompress('['+strjoin(string(arr_round([min(n)-(nlow/50D), max(n)+(nhi/50D)])),',')+']')
     muerange = strcompress('['+strjoin(string(arr_round([min(mue)-(muelow-min(mue)), max(mue)+(max(mue)-muehi)])),',')+']')
     rerange = strcompress('['+strjoin(string(arr_round([min(re)-(relow/50D), max(re)+(rehi/50D)])),',')+']')
     sbarange = strcompress('['+strjoin(string(arr_round([min(serE)-(sbalow/50D), max(serE)+(sbahi/50D)])),',')+']')
     sparange = strcompress('['+strjoin(string(arr_round([min(serpa)-(spalow/50D), max(serpa)+(spahi/50D)])),',')+']')
  endif
  if expnn gt 0 then begin
     exrange = strcompress('['+strjoin(string(arr_round([min(expx0)-(exlow/50D), max(expx0)+(exhi/50D)])),',')+']')
     eyrange = strcompress('['+strjoin(string(arr_round([min(expy0)-(eylow/50D), max(expy0)+(eyhi/50D)])),',')+']')
     mu0range = strcompress('['+strjoin(string(arr_round([min(mu0)-(mu0low-min(mu0)), max(mu0)+(max(mu0)-mu0hi)])),',')+']')
     hrange = strcompress('['+strjoin(string(arr_round([min(h)-(hlow/50D), max(h)+(hhi/50D)])),',')+']')
     ebarange = strcompress('['+strjoin(string(arr_round([min(expE)-(ebalow/50D), max(expE)+(ebahi/50D)])),',')+']')
     eparange = strcompress('['+strjoin(string(arr_round([min(exppa)-(epalow/50D), max(exppa)+(epahi/50D)])),',')+']')
  endif
  if cennn gt 0 then begin
     
     xrange = strcompress('['+strjoin(string(arr_round([min(x0)-(xlow/50D), max(x0)+(xhi/50D)])),',')+']')
     yrange = strcompress('['+strjoin(string(arr_round([min(y0)-(ylow/50D), max(y0)+(yhi/50D)])),',')+']')
    
  endif
  if skynn gt 0 then begin
     
     skyrange = strcompress('['+strjoin(string(arr_round([min(sky)-(skylow/50D), max(sky)+(skyhi/50D)])),',')+']')
     
    
  endif
  
  tot = fix(total(flist[0,*]))
  
  index = where(flist[0,*] eq 1)
  plist = nlist[index]
  nnlist = ['serX0','serY0','log(I$\downe$)','R$\downe$','n','b/a','$\theta$$\downPA, ser$','expX0','expY0','log(I$\down0$)','h','b/a','$\theta$$\downPA, exp$','X$\down0$','Y$\down0$','Sky background']
  pplist = nnlist[index]
  
  data = pars
  
  if sernn gt 0 then begin
     
     sorder = fix(sorder[0])
     data[sorder+6,*] = data[sorder+6,*];(acos(data[sorder+6,*])*180D/!DPI)
     ;data[sorder+6,*] = (acos(data[sorder+6,*])*180D/!DPI)
     data[sorder+3,*] = 10D^((data[sorder+3,*]))
     ;data[sorder+4,*] = 10D^((data[sorder+4,*]))
  endif
  if expnn gt 0 then begin
     eorder = fix(eorder[0])
     data[eorder+5,*] = data[eorder+5,*];(acos(data[eorder+5,*])*180D/!DPI)
     ;data[eorder+5,*] = (acos(data[eorder+5,*])*180D/!DPI)
     data[eorder+3,*] = 10D^((data[eorder+3,*]))
     
  endif
  data = data[index,*]
  triplot,data,bins=25,titles=pplist,info=info,/boxplot,/contour

  
  !P.MULTI = 0
  
  index = where(flist[0,*] eq 1)
  npars= n_elements(Pars)
  nfpars = n_elements(Pars[index])
  
;Trace history plots 
  RESTORE, dir_out+'TraceHistory.sav'
  RESTORE, dir_out+'ChainInfo.sav'
  adpstop = sumchain[0]
  transstop = sumchain[1]
  adapt_stop = sumchain[2]
  discard = sumchain[3]
  chainstop = sumchain[4]

   xx = range(0,(chainstop-adapt_stop),(chainstop-adapt_stop))
   colors = Round(cgScaleVector(Findgen(nrep), 0, 255))
   cgLoadCT, 34
   color = plotcolor(indgen(nrep))

  if (n_elements(nlist) le 7) then nplots = 5 else nplots = 9
  !p.multi = [0, 1,nplots]
  charsize = 0.8
  if (expnn eq 0) then begin
     
     miny = min(Xveclistdim[*,2,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,2,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,2,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,2,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,2,0:chainstop-adapt_stop],yr=[miny,maxy], yt='mue',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,2,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
     miny = min(10D^Xveclistdim[*,3,0:chainstop-adapt_stop])-(0.01*min(10D^Xveclistdim[*,3,0:chainstop-adapt_stop]))
     maxy = max(10D^Xveclistdim[*,3,0:chainstop-adapt_stop])+(0.01*max(10D^Xveclistdim[*,3,0:chainstop-adapt_stop]))
     cgplot, xx, 10D^(Xveclistdim[0,3,0:chainstop-adapt_stop]), yt='Re',xt='Iteration',yr=[miny,maxy],/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, 10D^(Xveclistdim[kk,3,0:chainstop-adapt_stop]),/overplot,color=color[kk]
     endfor
     
     
     miny = min(Xveclistdim[*,4,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,4,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,4,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,4,0:chainstop-adapt_stop]))
     cgplot, xx, (Xveclistdim[0,4,0:chainstop-adapt_stop]),yt='n',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, (Xveclistdim[kk,4,0:chainstop-adapt_stop]),/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,5,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,5,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,5,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,5,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,5,0:chainstop-adapt_stop],yr=[miny,maxy], yt='b/a',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,5,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,6,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,6,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,6,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,6,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,6,0:chainstop-adapt_stop],yt='ser PA',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,6,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor           
  endif else begin
     
     miny = min(Xveclistdim[*,2,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,2,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,2,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,2,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,2,0:chainstop-adapt_stop],yr=[miny,maxy], yt='mue',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,2,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
     miny = min(10D^Xveclistdim[*,3,0:chainstop-adapt_stop])-(0.01*min(10D^Xveclistdim[*,3,0:chainstop-adapt_stop]))
     maxy = max(10D^Xveclistdim[*,3,0:chainstop-adapt_stop])+(0.01*max(10D^Xveclistdim[*,3,0:chainstop-adapt_stop]))
     cgplot, xx, 10D^(Xveclistdim[0,3,0:chainstop-adapt_stop]), yt='Re',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, 10D^(Xveclistdim[kk,3,0:chainstop-adapt_stop]),/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,4,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,4,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,4,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,4,0:chainstop-adapt_stop]))
     cgplot, xx, (Xveclistdim[0,4,0:chainstop-adapt_stop]),yt='n',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, (Xveclistdim[kk,4,0:chainstop-adapt_stop]),/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,5,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,5,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,5,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,5,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,5,0:chainstop-adapt_stop],yr=[miny,maxy], yt='b/a',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,5,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,6,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,6,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,6,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,6,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,6,0:chainstop-adapt_stop],yt='ser PA',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,6,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     

     miny = min(Xveclistdim[*,9,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,9,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,9,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,9,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,9,0:chainstop-adapt_stop],yr=[miny,maxy], yt='mu0',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,9,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
           
     miny = min(10D^Xveclistdim[*,10,0:chainstop-adapt_stop])-(0.01*min(10D^Xveclistdim[*,10,0:chainstop-adapt_stop]))
     maxy = max(10D^Xveclistdim[*,10,0:chainstop-adapt_stop])+(0.01*max(10D^Xveclistdim[*,10,0:chainstop-adapt_stop]))
     cgplot, xx, (10D^Xveclistdim[0,10,0:chainstop-adapt_stop]),yt='h',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, (10D^Xveclistdim[kk,10,0:chainstop-adapt_stop]),/overplot,color=color[kk]
     endfor         
     
     miny = min(Xveclistdim[*,11,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,11,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,11,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,11,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,11,0:chainstop-adapt_stop],yr=[miny,maxy],yt='exp b/a',xt='Iteration',/nodata,charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,11,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
     miny = min(Xveclistdim[*,12,0:chainstop-adapt_stop])-(0.01*min(Xveclistdim[*,12,0:chainstop-adapt_stop]))
     maxy = max(Xveclistdim[*,12,0:chainstop-adapt_stop])+(0.01*max(Xveclistdim[*,12,0:chainstop-adapt_stop]))
     cgplot, xx, Xveclistdim[0,12,0:chainstop-adapt_stop],yt='exp PA',xt='Iteration',/nodata,yr=[miny,maxy],charsize=charsize
     for kk=0, nrep-1 do begin
        cgplot, xx, Xveclistdim[kk,12,0:chainstop-adapt_stop],/overplot,color=color[kk]
     endfor
     
  endelse
  
   
  !p.multi=0

 
;Autocorrelation plots 
;CALCULATE THE ACF and ESS

  Ncol = FIX(ceil(nfpars/2D))
  aa = 2D*Ncol
  aa = aa - n_elements(plist)
  PLOTSYM, 0 ,0.5, /FILL
     
     
  if tot gt 1 then !P.MULTI = [0,2,Ncol]
  lag = range(0,50,51)
  
  nlags = 50l
  ACF = make_array(nrep,nfpars,nlags+1,/DOUBLE)
  acfsig = make_array(nrep,nfpars,/DOUBLE)
  ESS = make_array(nrep,nfpars,/DOUBLE)

  items = make_array(nrep,/STRING)
  for j=0,nrep-1 do items[j] = strcompress('Chain '+string(j+1))
  
  for jj=0, nfpars-1 do begin
     
     fakex = 1.
     fakey = 1.
     cgplot, fakex, fakey, yr = [-0.4,1], xr=[-0.5,50],/nodata,Title = pplist[jj],xt='Lag',yt='ACF'
     if (jj eq 0) then Al_legend, items,psym=8,box=0,/bottom,/left,colors=color,charsize=0.5
     
     for kk=0, nrep-1 do begin
        ACF[kk,jj,*] =  autocorr(Xveclistdim[kk,index[jj],discard :chainstop-adapt_stop],nlags=nlags)
        acfsig[kk,jj] = 2D*sqrt((1D/(chainstop-adapt_stop+discard))*(1D + 2D*(TOTAL(ACF[kk,jj,*]^2D))))
        essindex = where(ACF[kk,jj,*] ge (2D*acfsig[kk,jj]))
        ESS[kk,jj] = DOUBLE(chainstop-adapt_stop+discard) / (1D + (2D * TOTAL((acf[kk,jj,*])[essindex])))
        cgplot, lag, ACF[kk,jj,*],/overplot,psym=8,color=color[kk]
     endfor
  endfor
  if aa eq 1 then cgplot,fakex,fakey,XSTYLE=4,YSTYLE=4,/nodata
  

  !P.multi=0
  
  
END





 
