
;+
; NAME:
;   BaGal2D
; PURPOSE:
;   A Bayesian MCMC method for 2D photometric decompositions of
;   galaxies and other astronomical objects.    
; 
; INPUTS:
;   
; KEYWORDS:
;
;
; OUTPUTS:
;
;
; BUGS:
; REVISION HISTORY:
;   Project Start: September 2014
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
     ;gauss =  psf_gaussian(NPIXEL = imag_dim,FWHM=[FWHM,FWHM],/normal)
     ;writefits,'./outputs/psf.fits',gauss
     return,gauss
  endif 
  ;Moffat PSF
  if (func eq 1) then begin
     alpha=fwhm/(2D *sqrt(2D ^(1D /betaa) -1D))
     moffat=((betaa-1D)/(!DPI*alpha^2D))*(1D +(r/alpha)^2D)^(-1D*betaa)
     moffat=reform(moffat,nx,ny)
     ;writefits,'/Users/ja66/Documents/Work/Alpha/Pro/mogal/Outputs/psf/psf.fits',moffat
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
if (n le 0) then return,result   ; error : data not defined
 
; check if speficic percentiles requested - if not: set standard
if(not keyword_set(value)) then value = [ 0., 0.25, 0.5, 0.75, 1.0 ]
 
; create a temporary copy of the data and sort
; tmp = data
; tmp = tmp(sort(tmp))
; NO: simply save the sorted index array
  ix = sort(data)
 
; loop through percentile values, get indices and add to result
; This is all we need since computing percentiles is nothing more
; than counting in a sorted array.
for i=0,n_elements(value)-1 do begin
 
   if(value(i) lt 0. OR value(i) gt 1.) then return,-1
 
;   if(value(i) le 0.5) then ind = fix(value(i)*n)    $
;   else ind = fix(value(i)*(n+1))
   if(value(i) le 0.5) then ind = long(value(i)*n)    $
   else ind = long(value(i)*(n+1))
   if (ind ge n) then ind = n-1    ; small fix for small n
                                   ; (or value eq 1.)
 
;  if(i eq 0) then result = tmp(ind)  $
;  else result = [result, tmp(ind) ]
; ## change number 2
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
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)


  Priserx0_temp = !Values.F_NAN
  Prisery0_temp = !Values.F_NAN
  Priexpx0_temp = !Values.F_NAN
  Priexpy0_temp = !Values.F_NAN

  if (f1dfit eq 0) then begin 
     skip = 1+i
     readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
              PriserX0, PriserY0, PriIe, PriRe, Prin, PriserE, PriserPA,$
              PriexpX0, PriexpY0, PriI0, Prih, PriexpE, PriexpPA,$
              FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',$
              SKIPLINE=skip,/SILENT,NUMLINE=1

     PriserPA = cos(PriserPA*!DPI/180D)
     PriexpPA = cos(PriexpPA*!DPI/180D)
     
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
  cenorder = where(nlist eq 'X0', cennn)

  if sernn le 0 then sflag = 1
  if expnn le 0 then eflag = 1
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

  if (sflag eq 0) or (eflag eq 0) or (BnDflag eq 0) or (cenflag eq 0) then flag = 0
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
  cenorder = where(nlist eq 'X0', cennn)

  if sernn le 0 then sflag = 1
  if expnn le 0 then eflag = 1
  if cennn le 0 then cenflag = 1

  ;flag = 1
  flag = [1,1,1,1,1]
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
     ;expPA = ACOS(pars[eorder+5])*108D/!DPI 
     expPA = (pars[eorder+5]) 
     if expPa lt -360D then eflag = 0
     if expPA gt 360D then eflag = 0
     
     ;if (eflag eq 0) then stop
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
     

     ;Bulge-to-disc ratio < Re
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

  ;if sflag eq 0 or eflag eq 0 or BnDflag eq 0 or cenflag eq 0 then flag = 0
  ;if flag eq 0 then stop
  return, flag
END

Pro bagal_priors,pars,varmat=varmat
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove  
  common mflags, fsersic, fexp
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma

  skip = 1+i
  readcol, './inputs/bagal_inputs.txt',NAME,seein,seein2,betaaa,$
           serX0, serY0, Ie, Re, n, serE, serPA,$
           expX0, expY0, I0, h, expE, expPA,$
           FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',$
           SKIPLINE=skip,/SILENT,NUMLINE=1
  
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
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
     ;if serPa lt 0. then serPA = 180D + serPa
     ;if serPA gt 180. then serPA = 180D - serPa
     ;serPA = cos(serPA*!DPI/180D)
     ;serPA = alog10(serPA)

     parsser = [serX0,serY0,Ie,Re,n,serE,serPA]
    
     pars = [pars,parsser]
     
     aser = [serx0/10.,sery0/10.,Ie/10.,Re/100.,n/100.,serE/10.,serPA/10.]
     dser = [serx0/50.,sery0/50.,Ie/50.,Re/500.,n/500.,serE/50.,serPA/50.]

     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 0.01D, 5D]
     ;parscale = [0.1D,0.1D, double(0.001), 0.0001D, 0.01D, 0.0005D, 0.5D]
     
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
     ;if expPa lt 0. then expPA = 180D + expPa
     ;if expPA gt 180. then expPA = 180D - expPa
     ;expPA = cos(expPA*!DPI/180D)
    ;expPA = alog10(expPA)

     parsexp = [expX0,expY0,I0,h,expE,expPA]
     pars = [pars,parsexp]
     aexp = [expx0/10.,expy0/10.,I0/10.,h/10.,expE/10.,expPA/10.]
     dexp = [expx0/100.,expy0/100.,I0/100.,h/100.,expE/100.,expPA/100.]
     
     parscale = [0.01D,0.01D, 0.01D, 0.01D, 0.01D, 5D]
     ;parscale = [0.5D,0.5D, double(0.001), 0.05D, 0.005D, 0.5D]
     varmat = [varmat,parscale]
     
  endif

;CENTRES
  if cennn gt 0 then begin 
     
    x0 = mean([serx0_temp,expx0_temp],/NAN)
    y0 = mean([sery0_temp,expy0_temp],/NAN)
       
    
    parscen = [X0,Y0]
    pars = [pars,parscen]

    ;parscale = [2.5e-3, 2.5e-3]
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
  astep[0] = af0
  if fsersic eq 1 then astep = [astep,aser]
  if fexp eq 1 then astep = [astep, aexp]
  if fsersic eq 1 then adel = [adel,dser]
  if fexp eq 1 then adel = [adel, dexp]
  astep = transpose([[astep],[adel]])



  REMOVE,0,pars
  REMOVE,0,varmat

  pri_pars = pars


  ;Check parameters 
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
           
        
           
     ;print, pflag
     
     ;table_ext,'./sample/COSMOS/i_F125/F125/COSMOS_F125_gal0011.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/SDSS/Gad09/iband/Disc/SDSS_i_gal0268.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     
     ;aa = -23.726900
     ;kk = 0.095942
     ;airmass = 1.1977736
     ;zcal= -(aa+kk*airmass)+2.5*alog10(53.907456*0.396^2.)
     ;cgplot, x,-2.5*ALOG10(y)+zcal,yr=[max(-2.5*ALOG10(y)+zcal),min(-2.5*ALOG10(y)+zcal)],xr=[0,ceil(max(x))],position=[0.1,0.1,0.5,0.5],yt='MAG',xt='Pixel'
     ;if (expnn eq 1) then begin
     ;   I0 = 10D^(pri_pars[9]) 
     ;   cgplot, x,-2.5*ALOG10(I0*exp(-x/10D^pri_pars[10]))+zcal,/overplot,color='blue'
     ;   d1exp = (I0*exp(-x/10D^pri_pars[10]))
     ;endif else d1exp = 0
     ;if (sernn eq 1) then begin      
     ;   Ie = 10D^(pri_pars[2])
     ;   invn = 1D / pri_pars[4]
     ;   bn = (1.9992 * pri_pars[4]) - 0.3271
     ;   cgplot, x,-2.5*ALOG10(Ie * EXP(-bn * ( ((x/10D^(pri_pars[3]))^(invn)) - 1D) ))+zcal,/overplot,color='red'
     ;   d1ser = Ie * EXP(-bn * ( ((x/10D^pri_pars[3])^(invn)) - 1D) )
     ;endif else d1ser = 0
     ;d1tot = -2.5*ALOG10(d1exp + d1ser) +zcal
     ;cgplot, x, d1tot,/overplot,color='purple'
     ;print, pri_pars
     

     if (expnn eq 0) and (sernn eq 1) then pflag=[1,1,1,1,1]
  ENDREP UNTIL (pflag[0] eq 1) and (pflag[1] eq 1) and (pflag[2] eq 1) and (pflag[3] eq 1) 
 ;stop
  pars = pri_pars
;stop

END


;*************************************************************************************************


PRO bagal_setup
  common setup1, tmcmc, bagal_like, nchain, nstep, tstep, pridis, f1dfit
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common exposure, pExpo,Expo,exponames,exponame
  common mflags, fsersic, fexp
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma
  
;Read the settings file 
  readcol, './inputs/setup.txt', index, inputs, bbb, format='A,F,A',$
           COMMENT='#',STRINGSKIP='[',/SILENT
 ;GENERAL MCMC INPUTS
  tmcmc = inputs[0]
  af0 = inputs[1]
  flike = inputs[2]
  case flike of
     0: bagal_like = 'like'
     1: bagal_like = 'likeBay'
     2: bagal_like = 'likePos'
     3: bagal_like = 'likeGauss'
     4: bagal_like = 'likeScaled'
  ENDCASE
  nchain = inputs[3]
  nstep = inputs[4]
  f1dfit = fix(inputs[5])
;IMAGE/INSTRUMENT PARAMETERS
  scale = inputs[6]
  zcalib = inputs[7]
  tpsf = fix(inputs[8])
  betaa = inputs[9]
  fFFT = inputs[10]
  seeing = inputs[11]
  seeing2 = inputs[12]
  sky = inputs[13]
  Pnoise = fix(inputs[14])
  gain = inputs[15]
  RON = inputs[16]
  pExpo = fix(inputs[17])
  Expo = DOUBLE(inputs[18])
  Bnoise = fix(inputs[19])
  skysig = inputs[20]
  tklog = inputs[21]
;COMPONENTS TO FIT
;Serśic profile / Bulge
  fsersic = fix(inputs[22])
;Exponential profile / Disc
  fexp = fix(inputs[23])
;Broken Exponential / Outer disc 
  fbexp = fix(inputs[24])
;Ferrors profile / Main Bar
  ffer = fix(inputs[25])
;Ferrors profile / Secondary Bar
  ffer2 = fix(inputs[26])
;King Profile 
  fking = fix(inputs[27])
;Nuclear Source 
  fnuc = fix(inputs[28])
;Gaussian Profile 		  
  fgau = fix(inputs[29])
;Moffat Profile
  fmof = fix(inputs[30])
;SPECIFICATIONS
  rmin = inputs[31]
  rmax = inputs[32]
  if rmin eq 0 then rmin = fix(rmin)
  sigover = inputs[33]
  fmask = fix(inputs[34])
  maskonly = fix(inputs[35])
  err_fac = fix(inputs[36])
  
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
;Specifications 
  
  if fcen[0] eq 1 then fixcen = 1 else fixcen = 0
  if psky[0] eq 1 then fsky = 1 else fsky = 0
 
;Create bagal file from user input file 
  BGmkfile, './inputs/'+inname

;CREATING THE PARAMETER MATRIX AND THE MODEL LIST 
 ;pser = strjoin(strtrim(string([fserx0, fsery0, fIE, fRe, fn, fserE, fserPA]),2),',')
 ;pexp = strjoin(strtrim(string([fexpx0, fexpy0, fI0, fh, fexpE, fexpE, fexpPA]),2),',')
  
  if (fixcen eq 1) then begin
     fx0 = fcen
     fy0 = fcen
     fserx0 = [0,0,0]
     fsery0 = [0,0,0]
     fexpx0 = [0,0,0]
     fexpy0 = [0,0,0]
  endif else begin
     fx0 = 0
     fy0 = 0    
  endelse

  pcen = [[fx0],[fy0]]
  pser = [[fserx0], [fsery0], [fIE], [fRe], [fn], [fserE], [fserPA]]
  pexp = [[fexpx0], [fexpy0], [fI0], [fh], [fexpE], [fexpPA]]

 

  serlist = ['serX0','serY0','Ie','Re','n','serE','serPA']
  explist = ['expX0','expY0','I0','h','expE','expPA']
  cenlist = ['X0','Y0']
  skylist = ['sky']

  ;The flag list
 ;flist = MAKE_ARRAY(1,/STRING)
  flist = MAKE_ARRAY(3,1,/INTEGER)
  if fsersic eq 1 then flist = [[flist],[pser]]
  if fexp eq 1 then flist = [[flist], [pexp]]
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
  if fixcen eq 1 then nlist = [nlist, cenlist]
  if fsky eq 1 then nlist = [nlist,skylist]
  REMOVE, 0, nlist
 
 ;par_matrix = transpose([[flist],[nlist]])



  
  
;The model routine 
  mlist = MAKE_ARRAY(1,/STRING)
  if fsersic eq 1 then mlist = [mlist, 'sersic(spars)']
  if fexp eq 1 then mlist = [mlist, 'expon(epars)']
  REMOVE, 0 , mlist
  mlist = strjoin(mlist,'+')
 

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

  ;ima[132:142,117:126] = !Values.F_NAN
  ;ima[86:113,139:172] = !Values.F_NAN
  ;ima[79:89,103:113] = !Values.F_NAN

;Find sky level and the good pixels to fit above 1sigma of sky
  newima = ima
  newima[(nx/2)-5:(nx/2)+5,(ny/2)-5:(ny/2)+5] = !Values.F_NAN
  meanclip, newima,skymean,skysigma, CLIPSIG=2,MAXITER=10
  
  sigabove = (skysigma+skymean)
  
  ;noisemap = readfits('./inputs/weights/candels_f125_blank.fits',h,/noscale,/SILENT)
  
  
  ;meanclip, (ivar^(-0.5)),skymean2,skysigma2,MAXITER=10
  ;sigma = make_array(nx,ny,/DOUBLE,value=((skysigma^2D)))
  ;sigma2 = make_array(nx,ny,/DOUBLE,value=(skymean2^2D))
  ;sigma3 = make_array(nx,ny,/DOUBLE,value=(mean([skymean2,skysigma])^2D))
  
  ;iindex = where_xyz(ima ge sigabove, xind=xind1, yind=yind1)
  ;shot = make_array(nx,ny,/DOUBLE)
  ;shot_noise = sqrt(ABS(ima[xind1,yind1])/1900D)

  ;shot[xind1,yind1] = shot_noise
    
  ;sigma4 =  sqrt(ABS(ima)/1900D)

  ;total_error = (shot^2D + (sigma))^(0.5)
  ;total_error2 = (shot^2D + sigma2)^(0.5)
  ;total_error3 = (shot^2D + (sigma3))^(0.5)
  ;total_error3 = (sigma4)
  ;total_error = (sigma2)^(0.5)
  ;ivar = total_error^(-2D)
  ;ivar = total_error2^(-2D)
  ;ivar = total_error4^(-2D)
  
  ;stop
  
  print, 'IMAGE INFORMATION: '
  print, strcompress('Sky mean: '+string(skymean))
  print, strcompress('Sky sigma: '+string(skysigma))
  
  
  xc = nx/2D
  yc = ny/2D

 

  ;!p.multi = [0,1,3]
 
  ;imasky = ima[100,0:*]
  ;nosky = where(imasky gt (sigabove))
  ;pix = range(0.,199.,200.)
  ;cgplot, pix, imasky 
  ;newpix = pix[nosky] 
  ;newima = imasky[nosky] 
  
  ;cgplot, newpix, newima, color='red', /overplot
  ;yy = range(sigabove,sigabove,200.)
  
  ;cgplot, pix, yy, color='blue', /overplot 

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
  ;pix = range(0.,n_elements(subima[0,*]),n_elements(subima[0,*]))
  ;cgplot, pix, subima[radi,*]

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
  ;pix = range(0.,n_elements(subima[0,*]),n_elements(subima[0,*]))
  ;cgplot, pix, subima[radi,*]
  ;!p.multi = 0
  imasubb = [fix(minx),fix(maxx),fix(miny),fix(maxy)]
  
  
  ;PSF 
  case tpsf of
     0: psf = psf_construct(ndimension=31,func=tpsf,fwhm=seeing)    
     1: psf = psf_construct(ndimension=31,func=tpsf,fwhm=seeing,betaa=betaa) 
     2: psf = psf_construct(ndimension=31,func=tpsf,fwhm=[seeing,seeing2]) 
     3: psf = readfits(psfname,h,/noscale,/silent)
  ENDCASE
  
  
  ;skysigma = skymean + skysigma
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
  
   
  ;writefits, './outputs/xsmall.fits', newx
  ;writefits, './outputs/xbig.fits', xbig
  ;writefits, './outputs/ysmall.fits', newy
  ;writefits, './outputs/ybig.fits',ybig

  
  xbig = xbig + 0.05
  ybig = ybig + 0.05
  
 
  rxbig = (-xbig * sin(pa)) + (ybig * cos(pa))
  rybig = (-xbig * cos(pa)) - (ybig * sin(pa))
  
  newr = sqrt(rxbig^2D + (rybig/q)^2D)
  ;newr = REBIN(newr, (xsize[1]), (xsize[2]))
  
  return, newr
END


FUNCTION bagal_model,pars
  common setup3, flist, nlist, mlist, astep, i, af0,groups
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen,seeing2
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove 
  common noise, sigmap, Pnoise,gain,RON,Bnoise,skysig,skymean,skysigma
  
  Is=0 & Iexp=0
  
  
  
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)
  
  excstr = make_array(1,/STRING)
  
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

     I0 = 10D^(pars[eorder+2])
     ;I0 = 10D^(-0.4 * (mu0 - zcalib - 2.5*alog10((scale^2D))))
     h = 10D^(pars[eorder+3])
     expell = pars[eorder+4]
     ;exppa = ACOS(pars[eorder+5])*180D/!DPI
     exppa = (pars[eorder+5])*!DPI/180D
;DEFINE THE RADIAL POSITION VECTORS 
     
     rexp = model_array(expxc,expyc,exppa,expell,expx,expy)
     
     Iexp = I0*EXP(-rexp/h)
     
     excstr = [excstr,'Iexp']
  endif 
  ;REMOVE, 0 , excstr
  ;excstr = strjoin(excstr,'+')
  ;exc = execute('Itot ='+excstr)
  Itot = Is + Iexp

  ;PRINT, 'MCMC took ', (SYSTIME(1) - T), ' Seconds to complete'
  

  
  
  
  
;RESAMPLE FOR THE STEEP GRADIENTS     
  
;GRADIENT BETWEEN EACH PIXEL
  
 gradima = gradient(Itot)
  gradgrad = gradient(gradima)
  grad_shot = sqrt(gradima)/sqrt(Itot)
  gradthres = 0.5 ;Threshold set where the gradient noise is half that of the shot noise or greater 
  gradindex = where_xyz(grad_shot ge gradthres, Xind=xind, Yind=yind)
  ;gradindex = where_xyz(grad_shot ge median(grad_shot), Xind=xind, Yind=yind)
  ;writefits, './outputs/test_presub.fits', Itot
  ;writefits, './outputs/test_grad1.fits', gradima
  ;writefits, './outputs/test_gradgrad1.fits', sqrt(gradima)/sqrt(Itot)
  ;gradimaold = gradima
 
  
  ;meanclip, grad_shot,gmean,gsigma, CLIPSIG=2,MAXITER=10
  ;cgdisplay,700,700
  ;cgplot, grad_shot[126,*]
  ;xx = range(gmean,gmean,n_elements(grad_shot[126,*]))
  ;yy = range(gmean+gsigma,gmean+gsigma,n_elements(grad_shot[126,*]))
  ;zz = range(median(grad_shot),median(grad_shot),n_elements(grad_shot[126,*]))
  ;cgplot, xx,color='red',/overplot
  ;cgplot, yy,color='blue',/overplot
  ;cgplot, zz,color='green',/overplot
  
  ;stop
  xc = serxc
  yc = seryc
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
  ;subgrad = gradima[minx:maxx,miny:maxy]

  
  ;abovegrad = where(gradima ge gradthres)
  ;Nabov = DOUBLE(n_elements(gradima[abovegrad]))
  ;abovbox = where(subgrad ge gradthres)
  ;Nbox = DOUBLE(n_elements(subgrad[abovbox]))
  ;frac = Nbox/Nabov
  
  ;if frac le 0.8 then begin
  ;   REPEAT BEGIN
  ;      circ = cir_mask(gradima, xc, yc, gradradi)
   ;     minx = xc - gradradi
  ;      maxx = xc + gradradi
  ;      miny = yc - gradradi
  ;      maxy = yc + gradradi
  ;      subgrad = gradima[minx:maxx,miny:maxy]
  ;      abovbox = where(subgrad ge gradthres)
  ;      Nbox = double(n_elements(subgrad[abovbox]))
  ;      frac = Nbox/Nabov
  ;      gradradi = gradradi + 1D
  ;   ENDREP UNTIL (gradradi eq ((nx/2D)-1)) OR (frac ge 0.8) or (gradradi eq ((ny/2D)-1))
  ;endif

  ;xmin = min(xind)
  ;xmax = max(xind)
  ;ymin = min(yind)
  ;ymax = max(yind)
  
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
  
  ;xmin = fix((serxc)-30D)
  ;xmax = fix((serxc)+30D)
  ;ymin = fix((seryc)-30D)
  ;ymax = fix((seryc)+30D)

  xsubdis = xsubmax - xsubmin
  ysubdis = ysubmax - ysubmin
  Imasub = [xsubmin,xsubmax,ysubmin,ysubmax]
  

  

  Isnew=0 & Iexpnew=0 
  if sernn gt 0 then begin
     sernewr = remodel(Imasub,serxc,seryc,serpa,serell,serx,sery,serxsize)
     ;ser_index = where_xyz(sernewr lt 0.01, Xind=xind_ser, Yind=yind_ser)
     ;sernewr[xind_ser,yind_ser] = 0.01D
     Isnew = Ie * EXP( -bn * (((sernewr/Re)^invn)-1D))
  endif
  if expnn gt 0 then begin
     expnewr = remodel(Imasub,expxc,expyc,exppa,expell,expx,expy,expxsize)
     ;exp_index = where_xyz(expnewr lt 0.01, Xind=xind_exp, Yind=yind_exp)
     ;expnewr[xind_exp,yind_exp] = 0.01D
     Iexpnew = I0 * EXP(-expnewr/h)
  endif
  Itotnew = Isnew + Iexpnew


  
  ;gradima = gradient(Itotnew)
  ;grad_shot = sqrt(gradima)/sqrt(Itotnew)
  ;gradindex = where_xyz(grad_shot gt gradthres, Xind=xind, Yind=yind)
  ;xmin = min(xind)
  ;xmax = max(xind)
  ;ymin = min(yind)
  ;ymax = max(yind)
  
  
  
  ;writefits, './outputs/test_prosub.fits', Itotnew
  ;writefits, './outputs/test_grad2.fits', gradima
  ;writefits, './outputs/test_gradgrad2.fits', grad_shot
  ;stop
  

  ;Can do PSF before or after rebinning 
  ;if fFFT eq 0 then begin
       
     ;Itotal1 = convol_fft(Itotnew,psf) ;,KERNEL_FFT=newpsf)
  ;   Itotal1 = convolve_j(Itotnew,psf)

  ;endif else begin
  ;   Itotal1 = convol(Itotnew,psf,/edge_wrap,/NAN)
  ;endelse
  ;Itotnew = Itotal1 
  Itotnew = rebin(Itotnew, xsubdis+1, ysubdis+1)
  ;Isnew = rebin(Isnew, (xmax-xmin)+1, (ymax-ymin)+1)
  ;Iexpnew = rebin(Iexpnew, (xmax-xmin)+1, (ymax-ymin)+1)
  
  Itot[xmin:xmax,ymin:ymax] = Itotnew[2:n_elements(itotnew[*,0])-3,2:n_elements(itotnew[0,*])-3]
  ;Is[xmin:(xmax-1),ymin:(ymax-1)] = Isnew[0:n_elements(isnew[*,0])-2,0:n_elements(isnew[0,*])-2]
  ;Iexp[xmin:(xmax-1),ymin:(ymax-1)] = Iexpnew[0:n_elements(iexpnew[*,0])-2,0:n_elements(iexpnew[0,*])-2]

  
 
  
  ;gradima = gradient(Itot)
  ;grad_shot = sqrt(gradima)/sqrt(Itot)
  ;gradindex = where_xyz(grad_shot gt gradthres, Xind=xind, Yind=yind)
  ;writefits, './outputs/test_prosub.fits', Itot
  ;writefits, './outputs/test_grad2.fits', gradima
  ;writefits, './outputs/test_gradgrad2.fits', grad_shot
  ;writefits, dir_out+'test_bulge.fits', Is
  ;writefits, dir_out+'test_disc.fits', Iexp
  ;Itot = Is + Iexp

  ;writefits, dir_out+'test_res.fits', itot - itotold
  ;writefits, dir_out+'test_resgrad.fits', gradima - gradimaold
  ;stop

  
  ;Is[xmin:(xmax-1),ymin:(ymax-1)] = Isnew[0:n_elements(isnew[*,0])-2,0:n_elements(isnew[0,*])-2]
  ;Iexp[xmin:(xmax-1),ymin:(ymax-1)] = Iexpnew[0:n_elements(iexpnew[*,0])-2,0:n_elements(iexpnew[0,*])-2]
  

  ;PRINT, 'MCMC took ', (SYSTIME(1) - T), ' Seconds to complete'
  ;writefits, dir_out+'test1.fits', Itot
  ;cgplot, (itot[100,*]),color='red',/overplot
  ;gradima = gradient(Itot)
  ;cgplot, ALOG(gradima[100,*]),position=[0.55,0.1,0.9,0.45],/noerase 
  ;stop
  ;ires = Itot - ima
  ;writefits, dir_out+'test_res.fits', Ires
  
  ;writefits, './outputs/test_disc.fits', Iexp
  ;writefits, './outputs/test_bulge.fits', Is
  ;PSF 
  ;case tpsf of
  ;   0: ;psf=psf_construct(ndimension=nx,func=tpsf,fwhm=seeing)    
  ;   1: ;psf=psf_construct(ndimension=nx,func=tpsf,fwhm=seeing,betaa=betaa) 
  ;   2: ;psf = readfits(psfname,h,/noscale,/silent)
  ;ENDCASE

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
  
  ;ires = ima - Itotal
  
  ;writefits, dir_out+'test2.fits', Itotal
  ;writefits, dir_out+'test3.fits', Itot
  ;writefits, dir_out+'test_res.fits', Ires
  ;stop

  ;
  ;writefits, dir_out+'bulge_test.fits', Is
  ;cgplot, ires[100,*]
 
;SELECTING THE GOOD POINTS TO USE IN THE LIKELIHOOD ESTIMATES
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
  ;boxsize = 100
  ;boxcoors = [pars[npars-2] - (boxsize/2),pars[npars-1] - (boxsize/2),pars[npars-2] + (boxsize/2),pars[npars-1] + (boxsize/2)]
  ;model[0:boxcoors[0],*] = !Values.F_NAN
  ;model[boxcoors[2]:*,*] = !Values.F_NAN
  ;model[*,0:boxcoors[1]] = !Values.F_NAN
  ;model[*,boxcoors[3]:*] = !Values.F_NAN

  ;T = SYSTIME(1) 
  ;PRINT, 'MCMC took ', (SYSTIME(1) - T), ' Seconds to complete'

  ;stop
  ;Number of data points 
  N = (N_ELEMENTS(IMA[xmin:xmax,ymin:ymax]))
  NoF = N


  ;ivar = 1D/ (model + (sky*gain) + (RON^2D))
  
  newivar = ivar;*6.25D ;newsigma^(-2D)
  ;newivar = ivar
  ;stop
  ;Normalise the model
  ;model = model ;/ B
  if (tklog eq 1) then model = ALOG(model) 

  ;Chi squared
  chisq = TOTAL((ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * (ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * newivar[xmin:xmax,ymin:ymax],/DOUBLE,/NAN)
  
  ;chisq = TOTAL((ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * (ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]),/DOUBLE,/NAN)
 
  ;Partition function
  isig = total(ALOG(sqrt((1D /newivar[xmin:xmax,ymin:ymax]))),/DOUBLE,/NAN)
  numbit = (N/2D) * ALOG(2D * !DPI)
  
  ;Likelihood 
  like = ((-0.5D * chisq) - isig - numbit ) 

  ;Poisson distribution
  ;ratio = (  (model[xmin:xmax,ymin:ymax] /  ima[xmin:xmax,ymin:ymax])  )
  ;like = TOTAL((-model[xmin:xmax,ymin:ymax] + (ima[xmin:xmax,ymin:ymax]*(ALOG(ratio) + 1D))) *  newivar[xmin:xmax,ymin:ymax],/DOUBLE,/NAN)


  ;print, like
  ;Priors
  Pri = 0d;call_function('priors',Pars)
  likee = like ;+ Pri
  
  if (FINITE(chisq, /NAN)) then stop
  
  
  ;if (likee/(-1*10d4)) lt 0.2 then stop
        
  ;print, strcompress('Reduced chisq: '+string(chisq/N))
  ;print, like1
  
  ;print, 10D^pars[6], 10D^pars[12]
  ;righthere
  stopp = 0
  if stopp eq 1 then begin
     
    
     ;table_ext,'./sample/MOGAL_se_igal0026.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/COSMOS/i_F125/F125/COSMOS_F125_gal0355.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/COSMOS/i_F160/F160/COSMOS_F160_gal0483.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/SDSS/Mock_gad09/iband/Disc/MOGAL_se_igal0026.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     table_ext,'./sample/SDSS/Gad09/iband/Disc/SDSS_i_gal0001.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     newima = ima[fix(nx/2D),*] ;-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D)
     
     
     plot = 1
     if plot eq 1 then begin

        sorder = where(nlist eq 'serX0', sernn)
        eorder = where(nlist eq 'expX0', expnn)
        redchi = -2D*(likee - pri + isig + numbit) / Nof  
        
        print, strcompress('Reduced chisq: '+string(chisq/N))
        print, strcompress('Reduced chisq2: '+string(redchi))

        ires = ima - model
        
        writefits, './outputs/test_image.fits', ima
        writefits, './outputs/test_res.fits', Ires
        writefits, './outputs/test.fits', model
        !p.multi=0


        aa = -23.726900
        kk = 0.095942
        airmass = 1.1977736

        zcal= -(aa+kk*airmass)+2.5*alog10(53.907456*0.396^2.)

        cgplot, x,-2.5*ALOG10(y)+zcal,yr=[max(-2.5*ALOG10(y)+zcal),min(-2.5*ALOG10(y)+zcal)],xr=[0,ceil(max(x))],position=[0.1,0.1,0.5,0.5],yt='MAG',xt='Pixel'
        ;cgplot,[1:(nx/2)],-2.5*ALOG10(model[pars[0],(nx/2):*])+zcal,/overplot,color='green'

        if (expnn eq 1) then begin
           ;I0 = 10D^(-0.4 * (pars[9] - 26.2303 - 2.5*alog10((0.06^2D)))) 
           I0 = 10D^(pars[9]) 
           cgplot, x,-2.5*ALOG10(I0*exp(-x/10D^pars[10]))+zcal,/overplot,color='blue'

           d1exp = (I0*exp(-x/10D^pars[10]))

        endif else d1exp = 0

        if (sernn eq 1) then begin
           ;Ie = 10D^(-0.4 * (pars[2] - 26.2303 - 2.5*alog10((0.06^2D))))
           Ie = 10D^(pars[2])
           invn = 1D / 10D^pars[4]
           bn = (1.9992 * 10D^pars[4]) - 0.3271
           cgplot, x,-2.5*ALOG10(Ie * EXP(-bn * ( ((x/10D^(pars[3]))^(invn)) - 1D) ))+zcal,/overplot,color='red'


           d1ser = Ie * EXP(-bn * ( ((x/10D^pars[3])^(invn)) - 1D) )

           
        endif else d1ser = 0
        
    
     
        d1tot = -2.5*ALOG10(d1exp + d1ser) +zcal

        cgplot, x, d1tot,/overplot,color='purple'
        
        ;cgplot, x,newmodel,color='green',/overplot
        
        cgplot,x,l,yr=[0,1],position=[0.1,0.6,0.4,0.9],yt='Ellipticity',/noerase,xt='Pixel'
        cgplot,x,pa,yr=[-180,180],position=[0.5,0.6,0.8,0.9],/noerase,xt='Pixel',yt='$\theta$$\downPA$'
        cgplot, x0,/noerase,color='blue',position=[0.8,0.3,0.99,0.5],yr=[97,103],title='Centre positions',charsize=0.9
        cgplot, y0,/noerase,color='red',/overplot
        xx0 = range(mean(x0),mean(x0),n_elements(x0))
        yy0 = range(mean(y0),mean(y0),n_elements(y0))
        xx = range(0,n_elements(x0),n_elements(x0))
        yy = range(0,n_elements(y0),n_elements(y0))
        cgplot,xx,xx0,linestyle=2,/overplot,color='blue' 
        cgplot,yy,yy0,linestyle=2,/overplot,color='red'
        print, mean(x0),mean(y0)
        
     
  ;print, chisq, like
  
        cgplot, ima[pars[0],*],/noerase,position=[0.55,0.1,0.75,0.3],xt='pixel',charsize=1
        cgplot, model[pars[0],*],color='red',/overplot
        item = 'Along y'
        al_legend,item,/top,/right,box=0
        cgplot, ima[*,pars[1]],/noerase,position=[0.55,0.3,0.75,0.5],charsize=1,xtickformat='(A1)'
        cgplot, model[*,pars[1]],color='red',/overplot
        item = 'Along x'
        al_legend,item,/top,/right,box=0
        stop
     endif
  endif

  return, likee
END

FUNCTION likeScaled, pars
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  
  npars = n_elements(pars)
  f = pars[npars]
  
  xmin = fix(imasubb[0])
  xmax =  fix(imasubb[1])
  ymin =  fix(imasubb[2])
  ymax =  fix(imasubb[3])
  
  
  model = bagal_model(pars)

  ;Normalise the model
  model = model / B
  
  ;Number of data points 
  N = (N_ELEMENTS(IMA[xmin:xmax,ymin:ymax]))

  ;Chi squared
  chisq = ((ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * (ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]) * ivar[xmin:xmax,ymin:ymax])
  
  ;Partition function
  isig = total(ALOG(sqrt((1./ivar[xmin:xmax,ymin:ymax]))),/DOUBLE,/NAN) 
  numbit = (N/2D) * ALOG(2D * !DPI)
  lnf = N * ALOG(f)

  ;Likelihood 
  invf = 1D/f
  like = (-0.5*invf*(TOTAL(chisq,/DOUBLE,/NAN)) - isig - numbit - lnf)

  ;Priors
  Pri = call_function('priors',Pars)

  like = like + Pri
  
 
  if (FINITE(chisq, /NAN)) then stop

  
  return, likee
END
Function likeGauss,pars
  common image_data, ima, ivar, nx, ny, ii_t, A, B, imasubb, radi,sigabove
  
  
  xmin = fix(imasubb[0])
  xmax =  fix(imasubb[1])
  ymin =  fix(imasubb[2])
  ymax =  fix(imasubb[3])
  
  model = bagal_model(pars)

  ;Normalise the model
  model = (model - A) / B
  
  chisq = ((ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax])*(ima[xmin:xmax,ymin:ymax] - model[xmin:xmax,ymin:ymax]))/(2D*model[xmin:xmax,ymin:ymax])
  
  like = TOTAL(( (-1D*chisq) - (0.5 * ALOG(2D*!DPI*model[xmin:xmax,ymin:ymax]))),/DOUBLE,/NAN)
  
  Pri = call_function('priors',Pars)
  
  likee = like + Pri
  
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
  print,' OBTAINING GUESSES .....................please wait    '




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
     
    
     
     
     ;disc[0,*] = Nprior(ndisc, disk[0], 1D, seed)
     ;disc[1,*] = Nprior(ndisc, disk[1], 1D, seed)
     
     

     
     ydisk2=disk[0]*exp(-(x/disk[1]))
     
     ;cgplot, x,-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D),yr=[max(-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D)),min(-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D))]
     ;cgplot, x,-2.5*ALOG10(ydisk2)+26.2303+2.5*ALOG10(0.06^2D),color='blue',/overplot
     ;cgplot, x,-2.5*ALOG10(ybulge)+26.2303+2.5*ALOG10(0.06^2D),color='red',/overplot
     
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
           ;lp=lp+0.1D
           ;cgplot, x,-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D),yr=[max(-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D)),min(-2.5*ALOG10(y)+26.2303+2.5*ALOG10(0.06^2D))]
           ;cgplot, x,-2.5*ALOG10(ydisk2)+26.2303+2.5*ALOG10(0.06^2D),color='blue',/overplot
           ;cgplot, x,-2.5*ALOG10(ybulge)+26.2303+2.5*ALOG10(0.06^2D),color='red',/overplot
           
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
     
        

     ;print, pflag
     ;table_ext,'./sample/COSMOS/i_F125/F125/COSMOS_F125_gal0251.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/SDSS/Gad09/iband/Disc/SDSS_i_gal0353.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0
     ;table_ext,'./sample/COSMOS/i_F160/F160/COSMOS_F160_gal0101.tab','sma,intens,ellip,pa,x0,y0',x,y,l,pa,x0,y0

     ;aa = -23.726900
     ;kk = 0.095942
     ;airmass = 1.1977736
     ;zcal= -(aa+kk*airmass)+2.5*alog10(53.907456*0.396^2.)
     ;   cgplot, x,-2.5*ALOG10(y)+zcal,yr=[max(-2.5*ALOG10(y)+zcal),min(-2.5*ALOG10(y)+zcal)],xr=[0,ceil(max(x))],position=[0.1,0.1,0.5,0.5],yt='MAG',xt='Pixel'
     ;   if (expnn eq 1) then begin
     ;      I0 = 10D^(pri_pars[9]) 
     ;      cgplot, x,-2.5*ALOG10(I0*exp(-x/10D^pri_pars[10]))+zcal,/overplot,color='blue'
     ;      d1exp = (I0*exp(-x/10D^pri_pars[10]))
     ;   endif else d1exp = 0
     ;   if (sernn eq 1) then begin      
     ;      Ie = 10D^(pri_pars[2])
     ;      invn = 1D / pri_pars[4]
     ;      bn = (1.9992 * pri_pars[4]) - 0.3271
     ;      cgplot, x,-2.5*ALOG10(Ie * EXP(-bn * ( ((x/10D^(pri_pars[3]))^(invn)) - 1D) ))+zcal,/overplot,color='red'
     ;      d1ser = Ie * EXP(-bn * ( ((x/10D^pri_pars[3])^(invn)) - 1D) )
     ;   endif else d1ser = 0
     ;   d1tot = -2.5*ALOG10(d1exp + d1ser) +zcal
     ;  cgplot, x, d1tot,/overplot,color='purple'
     ;print, pri_pars
     ;stop
     
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
     CASE tmcmc of
       ; 0: bagal_mcmc, seed,'bagal_step','likeBagal',nstep,pars,like,/log
       ; 1: bagal_ram, like, nstep, pars,'likeBagal',scale=varmat, $
       ;               adapt=nstep, accrate=af0, $
       ;               gamma=0.5, nstart=1
        2: PHI_AM,like, nstep, pars,'likeBay',scale=varmat, $
                      adapt=nstep, accrate=af0, $
                      nstart=1, nrep=nchain,info=info
     ENDCASE
     PRINT,''
     PRINT, 'MCMC took ', (SYSTIME(1) - T)/60D, ' Minutes to complete'
     
     ;parameters[i,*,*,*] = pars
     ;likelihood[i,*,*] = like
     sorder = where(nlist eq 'serX0', sernn)
     eorder = where(nlist eq 'expX0', expnn)

     if (sernn eq 1) and (expnn eq 1) then bit = '_se'
     if (sernn eq 0) and (expnn eq 1) then bit = '_e'
     if (sernn eq 1) and (expnn eq 0) then bit = '_s'

     ;model_pars = median(pars,DIMENSION = 2)
     ;model = bagal_model(model_pars)
     ;ires = ima - model
        
     ;writefits, './outputs/test_image.fits', ima
     ;writefits, './outputs/test_res.fits', Ires
     ;writefits, './outputs/test.fits', model
     
     if (sernn eq 1) then pars[4,*] = 10D^pars[4,*]


     output_data = make_array(n_elements(pars[*,0])+1,n_elements(pars[0,*]),/DOUBLE)
     for ii=0, n_elements(pars[0,*])-1 do output_data[*,ii] = [pars[*,ii],like[ii]]
     WRITEFITS, dir_out+NAME[i]+bit+'_dataout.fits', output_data
     
     information[i,*] = info
     xmin = fix(imasubb[0])
     xmax =  fix(imasubb[1])
     ymin =  fix(imasubb[2])
     ymax =  fix(imasubb[3])
     like = like;/(n_elements(ima[xmin:xmax,ymin:ymax]))
     ;stop
     oldlike = like
     sorder = where(nlist eq 'serX0', sernn)
     eorder = where(nlist eq 'expX0', expnn)
     cenorder = where(nlist eq 'X0', cennn)
     skyorder = where(nlist eq 'sky', skynn)

     discard = 0
     newparas = MAKE_ARRAY(1,/FLOAT)
     newlists = MAKE_ARRAY(1,/STRING)
     sorder = where(nlist eq 'serX0', sernn)
     eorder = where(nlist eq 'expX0', expnn)
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
        
        mue = reform(pars[sorder+2,discard:*])
        nlist[sorder+2] = 'mue'
                                ;mue = -2.5*alog10(Ie) + zcalib + 2.5*alog10((scale^2D))
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
        
        serpa =  reform(pars[sorder+6,discard:*])
        spainfo = postinfo(serpa)
        spamed = spainfo.med
        spalow = spainfo.siglow
        spahi = spainfo.sighigh
        newparas = [newparas,sxmed, sxlow, sxhi,symed, sylow, syhi,muemed,muelow,muehi,remed,relow,rehi,nmed,nlow,nhi,sbamed,sbalow,sbahi,spamed,spalow,spahi]
        newlists = [newlists,'sxmed', 'sxlow', 'sxhi','symed', 'sylow', 'syhi','iemed','ielow','iehi','remed','relow','rehigh','nmed','nlow','nhi','sbamed','sbalow','sbahi','spamed','spalow','spahi']
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
        
        mu0 = reform(pars[eorder+2,discard:*])
        nlist[eorder+2] = 'mu0'
                                ;mu0 = -2.5*alog10(I0) + zcalib + 2.5*alog10((scale^2D))
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

        exppa =  reform(pars[eorder+5,discard:*])
        epainfo = postinfo(exppa)
        epamed = epainfo.med
        epalow = epainfo.siglow
        epahi = epainfo.sighigh
        newparas = [newparas,exmed, exlow, exhi,eymed, eylow, eyhi,mu0med,mu0low,mu0hi,hmed,hlow,hhi,ebamed,ebalow,ebahi,epamed,epalow,epahi]
        newlists = [newlists,'exmed', 'exlow', 'exhi','eymed', 'eylow', 'eyhi','i0med','i0low','i0hi','hmed','hlow','hhigh','ebamed','ebalow','ebahi','epamed','epalow','epahi']
     endif
     REMOVE, 0 , newparas


     if sernn gt 0 then begin
        print, 'SERSIC COMPONENTS'
        sorder = where(newlists eq 'sxmed', sernn)
        print, $
           strcompress('X0 = ' + string(newparas[sorder]) + ' +/- ' + string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'symed', sernn)
        print, $
           strcompress('Y0 = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'iemed', sernn)
        print, $
           strcompress('Ie = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'remed', sernn)
        print, $
           strcompress('Re = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'nmed', sernn)
        print, $
           strcompress('n = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'sbamed', sernn)
        print, $
           strcompress('b/a = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
        sorder = where(newlists eq 'spamed', sernn)
        print, $
           strcompress('PA = ' + string(newparas[sorder]) + ' +/- ' +  string(((newparas[sorder] - newparas[sorder+1]) + (newparas[sorder+2] - newparas[sorder]))/2D))
     endif

if expnn gt 0 then begin
        print, 'EXPONENTIAL COMPONENTS'
        eorder = where(newlists eq 'exmed', expnn)
        print, $
           strcompress('X0 = ' + string(newparas[eorder]) + ' +/- ' + string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
        eorder = where(newlists eq 'eymed', expnn)
        print, $
           strcompress('Y0 = ' + string(newparas[eorder]) + ' +/- ' +  string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
        eorder = where(newlists eq 'i0med', expnn)
        print, $
           strcompress('mu0 = ' + string(newparas[eorder]) + ' +/- ' +  string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
        eorder = where(newlists eq 'hmed', expnn)
        print, $
           strcompress('h = ' + string(newparas[eorder]) + ' +/- ' +  string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
        eorder = where(newlists eq 'ebamed', expnn)
        print, $
           strcompress('b/a = ' + string(newparas[eorder]) + ' +/- ' +  string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
        eorder = where(newlists eq 'epamed', expnn)
        print, $
           strcompress('PA = ' + string(newparas[eorder]) + ' +/- ' +  string(((newparas[eorder] - newparas[eorder+1]) + (newparas[eorder+2] - newparas[eorder]))/2D))
     endif
;stop
  endfor
  ;parameters[*,2,*] = (parameters[*,2,*] * B)
  
  openw,2,dir_out+'PHI_info.dat',width=500
     
  for i=0, nima-1 do begin   
     printf,2,name[i],transpose(information[i,*])
  endfor 
  close,2


END

