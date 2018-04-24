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
end

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

function mrandomn, seed, covar, nrand, STATUS = status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;  NAME:
;     MRANDOMN
; PURPOSE:
; Function to draw NRAND random deviates from a multivariate normal
; distribution with zero mean and covariance matrix COVAR.
; 
; AUTHOR : Brandon C. Kelly, Steward Obs., Sept. 2004
;
; INPUTS : 
;
;    SEED - The random number generator seed, the default is IDL's
;           default in RANDOMN()
;    COVAR - The covariance matrix of the multivariate normal
;            distribution.    
; OPTIONAL INPUTS :
;
;    NRAND - The number of randomn deviates to draw. The default is
;            one.
; OUTPUT :
;
;    The random deviates, an [NRAND, NP] array where NP is the
;    dimension of the covariance matrix, i.e., the number of
;    parameters.
;
; OPTIONAL OUTPUT:
;     STATUS - status of the Cholesky decomposition.   If STATUS = 0 then 
;         the computation was successful.  If STATUS > 0 then the 
;         input covariance matrix is not positive definite  (see LA_CHOLDC),
;         and MRANDOMN
;         Note that if a STATUS keyword is supplied then no error message 
;         will be printed.
; REVISION HISTORY:
;     Oct. 2013  -- Use LA_CHOLDC instead of CHOLDC to enable use of STATUS
;           keyword.    W. Landsman
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() lt 2 then begin
    print, 'Syntax- Result = mrandomn( seed, covar, [nrand] , STATUS = )'
    return, 0
endif

printerr = ~arg_present(errmsg)
errmsg = '' 


;check inputs and set up defaults
if n_elements(nrand) eq 0 then nrand = 1
if size(covar, /n_dim) ne 2 then begin
    print, 'COVAR must be a matrix.'
    return, 0
endif

np = (size(covar))[1]
if (size(covar))[2] ne np then begin
    print, 'COVAR must be a square matrix.'
    return, 0
endif

epsilon = randomn(seed, nrand, np) ;standard normal random deviates (NP x NRAND matrix)

A = covar  ;store covariance into dummy variable for input into TRIRED

  la_choldc, A, /double, status=status        ;do Cholesky decomposition
  if status NE 0 then begin
     message,'Array is not positive definite, STATUS = ' + strtrim(status,2),/CON 
     return,-1
  endif   

for i = 0, np - 2  do A[i+1:*,i] = 0d        ;Zero out upper triangular portion

;transform standard normal deviates so they have covariance matrix COVAR
epsilon = A ## epsilon

return, epsilon
end

function Fractile, x, frac

;+
; FRACTILE
;	Return the requested fractile of the input data.
;
; Usage:
;	fr = fractile(x, frac)
;
; Return:
;	fr	<input>	The requested fractile.
;
; Arguments:
;	x	most	input	The array whose fractile(s) are to be
;				returned 
;	frac	float	input	The fractile(s) to return.
;
; Restrictions:
;	The input data must be a SORTable array (i.e. not complex,
;	string or structure).
;
; Example:
;	To find the interquartile range of a data set, try:
;	q = fractile(data, [.25,.75])
;	iqr = q(1)-q(0)
;	
; History:
;	Original: 26/9/95; SJT
;	Modify to interpolate: 4/2/97; SJT
;-

  if (n_params() ne 2) then message, 'Incorrect number of arguments'

  n = n_elements(x)-1
  i = sort(x)
  
  f0 = floor(frac*n)
  w1 = frac*n - f0
  f1 = floor(frac*n+1)
  w0 = f1 - frac*n

  return, x(i(f0))*w0 + x(i(f1))*w1
  
end

;+
; NAME:
;    SAMPLE_PDF
;
; PURPOSE:
;    This function randomly samples a given probability density function.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    Result = SAMPLE_PDF( Xpdf, Pdf, Nsample )
;
; INPUT:
;    Pdf:  A vector of type floating point containing the probability density 
;        function values.
;    Xpdf:  A vector of type floating point containing the location of the 
;        values in Pdf.
;    Nsample:  A scalar of type integer containing the sampling size.
;
; KEYWORD PARAMETERS:
;    -
;
; OUTPUT:
;    Result:  A vector of type floating point containing the randomly sampled 
;        locations.
;
; USES:
;    pdf_to_cdf.pro
;    var_type.pro
;
; PROCEDURE:
;    This function randomly samples the quantiles of a given probability 
;    density function.
;
; EXAMPLE:
;    Define the PDF values and where they are.
;      pdf = [ 0.1, 0.25, 0.3, 0.25, 0.1 ]
;      xpdf = [ 0., 1., 2., 4., 5. ]
;    Sample from the PDF 1000 times.
;      result = sample_pdf( xpdf, pdf, 1000 )
;    The result should have sampled 0 about 0.1*1000 times, etc.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@atm.ox.ac.uk), 2003-09-25
;    Modified:  DAS, 2011-04-13 (fixed bug by adding NORMALISE setting to call 
;        to pdf_to_cdf.pro;  modified documentation format)
;-

;***********************************************************************

FUNCTION SAMPLE_PDF, $
    Xpdf, Pdf, Nsample

;***********************************************************************
; Constants and Variables

; The size of the input PDF
npdf = n_elements( pdf )

; The location vector of the values in PDF
if n_elements( xpdf ) ne npdf then stop

;***********************************************************************
; Sample Randomly from the CDF

; Estimate the Cumulative Distribution Function
cdf = pdf_to_cdf( xpdf, pdf, xcdf=xcdf, normalise=1 )

; Initialise an output vector of random sampling from the uniform distribution
xsample = randomu( seed, nsample )

; Iterate through samples.
; Check if the counter has to be a long integer.
i0 = 0
if var_type( nsample ) ne 2 then i0 = long( i0 )
for i = i0, nsample - 1 do begin
  ; Find randomly sampled point on the CDF
  id = max( where( cdf le xsample[i] ) )
  ; Take the location of that sampled point.
  ; Note that we take the location of that segment defined for the PDF, not 
  ; for the CDF, as this is the more representative location of the segment.
  xsample[i] = xpdf[id+1]
endfor

;***********************************************************************
; The End

return, xsample
END

Function twoD_PDF,xs, bins=bins, range=range,weights=weights, $
                           maxins=maxins,minins=minins,locx=locx,locy=locy,info=info

  
  
  x = xs[0,*]
  y = xs[1,*]
  
  if (n_elements(range) eq 0) then begin
     
     ran1 = [min(x),max(x)]
     ran2 = [min(y),max(y)]
     
     if n_elements(maxins) eq 0 then xmaxin = max(x) else xmaxin = maxins[ii]
     if n_elements(minins) eq 0 then xminin = min(x) else xminin = minins[ii]
     
     rangex = Max(x < xmaxin, /NAN) - Min(x > xminin, /NAN)
     
     if n_elements(maxins) eq 0 then ymaxin = max(y) else ymaxin = maxins[jj]
     if n_elements(minins) eq 0 then yminin = min(y) else yminin = minins[jj]
     
     rangey = Max(y < ymaxin, /NAN) - Min(y > yminin, /NAN)
     
     
  endif else begin
     ran1 = range[0,*]
     ran2 = range[1,*]
     
     if n_elements(maxins) eq 0 then xmaxin = ran1[1] else xmaxin = maxins[ii]
     if n_elements(minins) eq 0 then xminin = ran1[0] else xminin = minins[ii]
     
     ;rangex = Max(x < xmaxin, /NAN) - Min(x > xminin, /NAN)
     rangex = xmaxin - xminin

     if n_elements(maxins) eq 0 then ymaxin = ran2[1] else ymaxin = maxins[jj]
     if n_elements(minins) eq 0 then yminin = ran2[0] else yminin = minins[jj]
     
     ;rangey = Max(y < ymaxin, /NAN) - Min(y > yminin, /NAN)
     rangey = ymaxin - yminin
  endelse
  
  if n_elements(bins) gt 0 then (ynbin = bins) else (ynbin = [])
  if n_elements(bins) gt 0 then (xnbin = bins) else (xnbin = [])
  
  IF N_Elements(xnbin) EQ 0 THEN BEGIN ; Scott's Choice
     xbinsize = (3.5D * StdDev(xminin > x < xmaxIn, /NAN))/N_Elements(x)^(1./3.0D)
     xnbin = (rangex / xbinsize ) + 1
  ENDIF ELSE BEGIN
     xbinsize = rangex / (xnbin -1D)
  ENDELSE
  
 
  
  
  IF N_Elements(ynbin) EQ 0 THEN BEGIN ; Scott's Choice
     ybinsize = (3.5D * StdDev(xminIn > x < xmaxIn, /NAN))/N_Elements(x)^(1./3.0D)
     ynbin = (rangey / ybinsize ) + 1
  ENDIF ELSE BEGIN
     ybinsize = rangey / (ynbin -1D)
  ENDELSE
  
 
  
  density = Hist_2D(x, y, Min1=xminin, Max1=xmaxin, Bin1=xbinsize, $
                    Min2=yminin, Max2=ymaxin, Bin2=ybinsize)
  
 
  info = {min1:xminin,max1:xmaxin,bin1:xbinsize,min2:yminin,max2:ymaxin,bin2:ybinsize}
  

  return, density
end

Function SAMPLE_PDF_ND,data,nsample=nsample
  
  
  
  ;Obtain the 2D PDF
  D2_Pdf = twoD_pdf(data,bins=25,info=info)
  
  ;Calculate the CDF grid
  histdata = Histogram(data[0,*], $
                       BINSIZE=info.bin1, $
                       MAX=info.max1, $
                       MIN=info.min1, $
                       LOCATIONS= xloc)
  histdata = Histogram(data[1,*], $
                          BINSIZE=info.bin2, $
                          MAX=info.max2, $
                          MIN=info.min2, $
                          LOCATIONS= yloc)

  xcdf = pdf_to_cdf( xloc, d2_pdf[*,0], xcdf=xcdf, normalise=1 )

  ; Initialise an output vector of random sampling from the uniform distribution
  xsample = randomu( seed, nsample )

; Iterate through samples.
; Check if the counter has to be a long integer.
  i0 = 0
  if var_type( nsample ) ne 2 then i0 = long( i0 )
  for i = i0, nsample - 1 do begin
  ; Find randomly sampled point on the CDF
     id = max( where( cdf le xsample[i] ) )
  ; Take the location of that sampled point.
  ; Note that we take the location of that segment defined for the PDF, not 
  ; for the CDF, as this is the more representative location of the segment.
     xsample[i] = xpdf[id+1]
     
     ycdf = pdf_to_cdf( xloc, d2_pdf[*,*], xcdf=ycdf, normalise=1 )
     
  endfor
  
  
  
  stop
  
  
END


PRO new_sample


  ;Read in Gadotti '09 results
  readcol, 'gad09.txt',name, xc, yc, sky, seeing_Gadotti,seeing_lower, seeing_higher, seeing_mean, seeing_1_G2, seeing_1_G2,seeing_Effwidth ,gain, Dark_current, zero_point, airmass, zcal,mue, Ie_gad, re_gad, n_gad, mu0, I0_gad, h_gad, z,format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',SKIPLINE=1

  ;Read in GASP2D results 
  readcol, '2Dfit_g09.txt',GALAXY,Ie_g2d,Re_g2d, n_g2d,ba_g2d,PA_g2d, BT_g2d,I0_g2d, h_g2d,ba_g2d,PA_g2d,DT_g2d,Xcenter_g2d,Ycenter_g2d,Chi_g2d,format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F',skipline=1
 

  RESOLVE_ROUTINE, 'cosmology_routines'
  red,omega0=0.3,omegalambda=0.7,h100=0.7D

   
  
 ;-----------------------
  dist=DANGULAR(z)
  physcale=dist/206265D
  scale=0.3961D
  band='i'
 ;----------------------
  
 


  ;Read in posterior distribution   
  path = './outputs/'
  results = file_search(path+'SDSS*', count=nima)
  ngals = n_elements(results)
  
  res0 = readfits(results[0],/SILENT)
  ressize = size(res0, /DIMENSION)
  npars = ressize[0] 
  
  n_runs = 10000 
  
  meds = make_array(ngals,npars,/DOUBLE)
  q = make_array(ngals,npars,6,/DOUBLE)
  low_3s = make_array(ngals,npars,/DOUBLE)
  low_2s = make_array(ngals,npars,/DOUBLE)
  low_1s = make_array(ngals,npars,/DOUBLE)
  
  high_1s = make_array(ngals,npars,/DOUBLE)
  high_2s = make_array(ngals,npars,/DOUBLE)
  high_3s = make_array(ngals,npars,/DOUBLE)

  err = make_array(ngals,npars,/DOUBLE)

  
  
  ie = make_array(ngals,/DOUBLE)
  re = make_array(ngals,/DOUBLE)
  n = make_array(ngals,/DOUBLE)
  serE = make_array(ngals,/DOUBLE)
  serPA = make_array(ngals,/DOUBLE)
  i0 = make_array(ngals,/DOUBLE)
  h = make_array(ngals,/DOUBLE)
  expE = make_array(ngals,/DOUBLE)
  expPA = make_array(ngals,/DOUBLE)
  x0 = make_array(ngals,/DOUBLE)
  y0 = make_array(ngals,/DOUBLE)
  

  
  newRe = make_array(n_elements(results),n_runs,/DOUBLE)
  newh = make_array(n_elements(results),n_runs,/DOUBLE)
  
  
  for i=0l, n_elements(results)-1 do begin
     ;Read in distributions 
        
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     

        ;data = 10D^res[[3,10],*]
     ;Random points from PDF
     ir = floor(Uprior(n_runs, 0, n_elements(Ie), seed))
     
     newRe[i,*] = re[ir]
     newh[i,*] = h[ir]
     
     ;Inverse transforming sampling 
     ;newsample = SAMPLE_PDF_ND(data,nsample=1000)
        
       
     for jj=0, npars-1 do begin
        quartile = fractile(res[jj,*], [0.0015,0.025,0.15865,0.84135,0.97635,0.9985])
        q[i,jj,*] = quartile
        meds[i,jj] = median(res[jj,*],/EVEN,/DOUBLE)  
     endfor
     
        
  endfor
  
  index = [3,10]
  meds[*,index] = 10D^meds[*,index]
  q[*,index,*] = 10D^(q[*,index,*])

  ;Converting Re and h into physical units
  for i=0l, n_elements(results)-1 do begin 
     newre[i,*] = (newre[i,*]*scale)*physcale[i]/1000.
  endfor
  rePhys_gad = (re_gad*scale)*physcale/1000.
  rePhys_g2d = (re_g2d*scale)*physcale/1000.
  rePhys = (meds[*,3]*scale)*physcale/1000.
  rePhys_q1 = (q[*,3,2]*scale)*physcale/1000.
  rePhys_q2 = (q[*,3,3]*scale)*physcale/1000.

  for i=0l, n_elements(results)-1 do begin
     newh[i,*] = (newh[i,*]*scale)*physcale[i]/1000.
  endfor
  hPhys_gad = (h_gad*scale)*physcale/1000.
  hPhys_g2d = (h_g2d*scale)*physcale/1000.
  hPhys = (meds[*,10]*scale)*physcale/1000.
  hPhys_q1 = (q[*,10,2]*scale)*physcale/1000.
  hPhys_q2 = (q[*,10,3]*scale)*physcale/1000.
  
  n_meds = meds[*,4] 
  nind = where(n_meds ge 2)

  lin_fit = make_array(n_runs,2,/double)
  pears = make_array(n_runs,/double)
  spear = make_array(n_runs,/double)
  kendal = make_array(n_runs,/double)

  for j=0l, n_runs-1 do begin
     
     result = LINFIT(newre[*,j], newh[*,j])
     lin_fit[j,*] = result
     pears[j] = CORRELATE(newre[*,j], newh[*,j])
     result = R_CORRELATE(newre[*,j], newh[*,j])
     spear[j] = result[0]
     result = R_CORRELATE(newre[*,j], newh[*,j], /KENDALL)
     kendal[j] = result[0]

  endfor
  
  
  ;Same again but generating random numbers from gaussian
  ;Using the sigmas from the posterior distributions 

  ;n_runs = 10000
  GauRe = make_array(n_elements(results),n_runs,/DOUBLE)
  Gauh = make_array(n_elements(results),n_runs,/DOUBLE)
  
  for i=0l, n_elements(results)-1 do begin
     sig_re = median([(rephys[i] - rephys_q1[i]),(rephys_q2[i] - rephys[i])])
     Gaure[i,*] = Nprior(n_runs, rephys[i], sig_re)
     sig_h = median([(hphys[i] - hphys_q1[i]),(hphys_q2[i] - hphys[i])])
     Gauh[i,*] = Nprior(n_runs, hphys[i], SIG_h)
     
  endfor


  Gau_fit = make_array(n_runs,2,/double)
  Gau_pears = make_array(n_runs,/double)
  Gau_spear = make_array(n_runs,/double)
  Gau_kendal = make_array(n_runs,/double)

  for j=0l, n_runs-1 do begin
     
     result = LINFIT(Gaure[*,j], Gauh[*,j])
     Gau_fit[j,*] = result
     Gau_pears[j] = CORRELATE(Gaure[*,j], Gauh[*,j])
     result = R_CORRELATE(Gaure[*,j], Gauh[*,j])
     Gau_spear[j] = result[0]
     result = R_CORRELATE(Gaure[*,j], Gauh[*,j], /KENDALL)
     Gau_kendal[j] = result[0]

  endfor
  
  cgLoadCT, 33
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  ;2D hist for sdss sample normalised 
  data = make_array(2,N_ELEMENTS(newre),/DOUBLE)
  data[0,*] = (newre[*,*])
  data[1,*] = (newh[*,*])
  density = twoD_pdf(data,bins=500,info=info)
  good_info = info
  master_density = density
  master_density[*,*] = 0D

  

  for i=0, n_elements(results)-1 do begin
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     ran = make_array(2,2)
     ran[0,*] = [good_info.min1,good_info.max1]
     ran[1,*] = [good_info.min2,good_info.max2]

     re = (re*scale)*physcale[i]/1000.
     h = (h*scale)*physcale[i]/1000.

     new_data = make_array(2,N_ELEMENTS(re),/DOUBLE)
     new_data[0,*] = (re)
     new_data[1,*] = (h)
     new_dist = twoD_pdf(new_data,bins=500, range=ran, info=newinfo)
     ;new_dist = (new_dist/total(new_dist))* ((newinfo.bin1)*(newinfo.bin2))
     
     maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
     scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
     new_dist = scaledDensity
     master_density = master_density + new_dist
     
  endfor
  
  ;maxDensity = Ceil(Max(Master_density)/1e2) * 1e2
  ;scaledDensity = BytScl(Master_density, Min=0, Max=maxDensity)
  
  
  
  cgDisplay , 700.,700.,/aspect
  cgPS_Open,strcompress('./plots/SDSS_corr.ps')

  ;cgplot, newre, newh, psym=3, /xlog, /ylog,xrange=[0.1,5],yrange=[0.5,20], position=[0.1,0.5,0.5,0.9],xt='R$\downe$ (/Kpc)', yt='h (/Kpc)',title='$\Phi$-PDF',charsize=0.9
  
  cgImage, master_Density, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
           XTitle='R$\downe$ (/Kpc)', YTitle='h (/Kpc)',title='$\Phi$-PDF', $
           Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9
  
  ;for i=0, n_runs-1 do begin
     x = range(min(newre[*,*]), max(newre[*,*]), 1000)
     y = median(lin_fit[*,0]) + (median(lin_fit[*,1]) * x)
     cgplot, (x),(y),color='red',/overplot
  ;endfor
  
  data = make_array(2,N_ELEMENTS(Gaure),/DOUBLE)
  data[0,*] = (Gaure[*,*])
  data[1,*] = (Gauh[*,*])
  density = twoD_pdf(data,bins=500, range=ran, info=good_info)

  master_density = density
  master_density[*,*] = 0D

  for i=0, n_elements(results)-1 do begin
     
     Re = Gaure[i,*]
        
     h= Gauh[i,*]
     
     ran = make_array(2,2)
     ran[0,*] = [good_info.min1,good_info.max1]
     ran[1,*] = [good_info.min2,good_info.max2]

     new_data = make_array(2,N_ELEMENTS(re),/DOUBLE)
     new_data[0,*] = (re)
     new_data[1,*] = (h)
     new_dist = twoD_pdf(new_data,bins=500, range=ran, info=newinfo)
     ;new_dist = (new_dist/total(new_dist))* ((newinfo.bin1)*(newinfo.bin2))
     
     maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
     scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
     new_dist = scaledDensity
     master_density = master_density + new_dist
     
  endfor
  
  ;maxDensity = Ceil(Max(Master_density)/1e2) * 1e2
  ;scaledDensity = BytScl(Master_density, Min=0, Max=maxDensity)
  
  
  
 
  ;cgplot, newre, newh, psym=3, /xlog, /ylog,xrange=[0.1,5],yrange=[0.5,20], position=[0.1,0.5,0.5,0.9],xt='R$\downe$ (/Kpc)', yt='h (/Kpc)',title='$\Phi$-PDF',charsize=0.9
  
  
  cgImage, master_Density, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
           XTitle='R$\downe$ (/Kpc)', YTitle='', title='N-PDF',$
           Position=[0.55,0.5,0.95,0.9],$
           charsize=0.9,/noerase
  ;for i=0, n_runs-1 do begin
     x = range(min(Gaure[*,*]), max(Gaure[*,*]), 1000)
     y = median(Gau_fit[*,0]) + (median(Gau_fit[*,1]) * x)
     cgplot, (x),(y),color='red',/overplot
  ;endfor
  
  cghistoplot, lin_fit[*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='a ($\Phi$-PDF)'

  cghistoplot, lin_fit[*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='b ($\Phi$-PDF)'

  cghistoplot, Gau_fit[*,0], /outline, position=[0.55,0.3,0.7,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='a (N-PDF)'

  cghistoplot, Gau_fit[*,1], /outline, position=[0.8,0.3,0.95,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='b (N-PDF)'
  
  cghistoplot, pears,/outline, position=[0.225,0.1,0.375,0.25],xrange=[min([min(pears),min(spear),min(kendal)]), max([max(pears),max(spear),max(kendal)])],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='r,$\rho$,$\tau$ ($\Phi$-PDF)'
  cghistoplot, spear,/outline,/oplot,color='blue'
  cghistoplot, kendal,/outline,/oplot,color='dark green'

   cghistoplot, Gau_pears,/outline, position=[0.675,0.1,0.825,0.25],xrange=[min([min(Gau_pears),min(Gau_spear),min(Gau_kendal)]), max([max(Gau_pears),max(Gau_spear),max(Gau_kendal)])],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='r,$\rho$,$\tau$ (N-PDF)'
  cghistoplot, Gau_spear,/outline,/oplot,color='blue'
  cghistoplot, Gau_kendal,/outline,/oplot,color='dark green'


  
;*********************************************************************************************
;Looking at the difference between n<2 and n>2 galaxies
;C - classical bulge
;P - pseudo bulgle 


  n_meds = meds[*,4] 
  Cind = where(n_meds ge 2)
  Pind = where(n_meds lt 2)
  
  lin_fit_C = make_array(n_runs,2,/double)
  pears_C  = make_array(n_runs,/double)
  spear_C  = make_array(n_runs,/double)
  kendal_C  = make_array(n_runs,/double)
  
  lin_fit_P = make_array(n_runs,2,/double)
  pears_P  = make_array(n_runs,/double)
  spear_P  = make_array(n_runs,/double)
  kendal_P  = make_array(n_runs,/double)

  for j=0l, n_runs-1 do begin
     
     result = LINFIT(newre[Cind,j], newh[Cind,j])
     lin_fit_C[j,*] = result
     pears_C[j] = CORRELATE(newre[Cind,j], newh[Cind,j])
     result = R_CORRELATE(newre[Cind,j], newh[Cind,j])
     spear_C[j] = result[0]
     result = R_CORRELATE(newre[Cind,j], newh[Cind,j], /KENDALL)
     kendal_C[j] = result[0]

     result = LINFIT(newre[Pind,j], newh[Pind,j])
     lin_fit_P[j,*] = result
     pears_p[j] = CORRELATE(newre[pind,j], newh[pind,j])
     result = R_CORRELATE(newre[pind,j], newh[pind,j])
     spear_p[j] = result[0]
     result = R_CORRELATE(newre[pind,j], newh[pind,j], /KENDALL)
     kendal_p[j] = result[0]

  endfor
  
  
  ;Same again but generating random numbers from gaussian
  ;Using the sigmas from the posterior distributions 
  
  cgLoadCT, 33
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  ;2D hist for sdss sample normalised 
  data = make_array(2,N_ELEMENTS(newre),/DOUBLE)
  data[0,*] = (newre[*,*])
  data[1,*] = (newh[*,*])
  density = twoD_pdf(data,bins=500,info=info)
  good_info = info
  master_density_C = density
  master_density_C[*,*] = 0D
  master_density_P = master_density_C
  nTot_C = []
  nTot_P = []

  for i=0, n_elements(results)-1 do begin
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     ran = make_array(2,2)
     ran[0,*] = [good_info.min1,good_info.max1]
     ran[1,*] = [good_info.min2,good_info.max2]

     re = (re*scale)*physcale[i]/1000.
     h = (h*scale)*physcale[i]/1000.

     new_data = make_array(2,N_ELEMENTS(re),/DOUBLE)
     new_data[0,*] = (re)
     new_data[1,*] = (h)
     new_dist = twoD_pdf(new_data,bins=500, range=ran, info=newinfo)
     ;new_dist = (new_dist/total(new_dist))* ((newinfo.bin1)*(newinfo.bin2))
     
     maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
     scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
     new_dist = scaledDensity
     
     if (n_meds[i] ge 2) then begin
        master_density_C = master_density_C + new_dist
        nTot_C = [nTot_C, n] 
     endif else begin
        master_density_P = master_density_P + new_dist
        nTot_P = [nTot_P, n] 
     endelse
  endfor
  
  ;maxDensity = Ceil(Max(Master_density)/1e2) * 1e2
  ;scaledDensity = BytScl(Master_density, Min=0, Max=maxDensity)
  
  
  
  ;cgplot, newre, newh, psym=3, /xlog, /ylog,xrange=[0.1,5],yrange=[0.5,20], position=[0.1,0.5,0.5,0.9],xt='R$\downe$ (/Kpc)', yt='h (/Kpc)',title='$\Phi$-PDF',charsize=0.9
  
  cgImage, master_Density_C, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
           XTitle='R$\downe$ (/Kpc)', YTitle='h (/Kpc)',title='n $\geq$ 2', $
           Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9
  
  ;for i=0, n_runs-1 do begin
     x = range(min(newre[cind,*]), max(newre[cind,*]), 1000)
     y = median(lin_fit_C[*,0]) + (median(lin_fit_C[*,1]) * x)
     cgplot, (x),(y),color='red',/overplot
  ;endfor
  
  
  
  
  cgImage, master_Density_P, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
           XTitle='R$\downe$ (/Kpc)', YTitle='', title='n < 2',$
           Position=[0.55,0.5,0.95,0.9],$
           charsize=0.9,/noerase
  ;for i=0, n_runs-1 do begin
  x = range(min(newre[pind,*]), max(newre[pind,*]), 1000)
  y = median(lin_fit_p[*,0]) + (median(lin_fit_p[*,1]) * x)
  cgplot, (x),(y),color='red',/overplot
  ;endfor
  
  cghistoplot, lin_fit_C[*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='a '

  cghistoplot, lin_fit_C[*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='b '

  cghistoplot, lin_fit_P[*,0], /outline, position=[0.55,0.3,0.7,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='a '

  cghistoplot, lin_fit_P[*,1], /outline, position=[0.8,0.3,0.95,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='b '
  
  cghistoplot, pears_C,/outline, position=[0.225,0.1,0.375,0.25],xrange=[min([min(pears_C),min(spear_C),min(kendal_C)]), max([max(pears_C),max(spear_C),max(kendal_C)])],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='r,$\rho$,$\tau$'
  cghistoplot, spear_C,/outline,/oplot,color='blue'
  cghistoplot, kendal_C,/outline,/oplot,color='dark green'

   cghistoplot, pears_P,/outline, position=[0.675,0.1,0.825,0.25],xrange=[min([min(pears_P),min(spear_P),min(kendal_P)]), max([max(pears_P),max(spear_P),max(kendal_P)])],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='r,$\rho$,$\tau$ '
  cghistoplot, spear_P,/outline,/oplot,color='blue'
  cghistoplot, kendal_P,/outline,/oplot,color='dark green',binsize=0.008


  cghistoplot, nTot_c, /outline
  cghistoplot, nTot_p, /outline, /oplot, color='blue'
 

  cgPS_close
  cgPS2PDF, strcompress('./plots/SDSS_corr.ps')


  stop
END

function plotcolor,x

  xnew = round(255*(x - min(x))/(max(x)-min(x)))

  col = byte(xnew)
  return,col
  
end

pro test_corr
  
  covar = make_array(2,2,/double)
  covar[0,0] = 1D
  covar[1,0] = 0.9D
  covar[0,1] = 0.9D
  covar[1,1] = 1D
  
  
  ran = mrandomn(seed, covar, 100000)
  new_dist = twoD_pdf(transpose(ran),bins=100,info=good_info)
    
  maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
  scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
  
  cgLoadCT, 33
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  ;cgImage, scaledDensity, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
  ;         XTitle='a', YTitle='b', title='N-PDF',$
  ;         charsize=0.9
  
  n = 100
  n_runs = 1000
  a = 2D
  b = 1D
  x = range(0,10,n)
  y = a + (b * x)
  sig_x = [0.1,1D]
  sig_y = [0.1D,1D]

  ;sig_x = [range(0.1,0.9,10),range(1,5,5)]
  ;sig_y = [range(0.1,0.9,10),range(1,5,5)]

  sig_x = [0.1D,0.25D,0.5D,0.75D,1D,2D,5D,10D]
  sig_y = [0.1D,0.25D,0.5D,0.75D,1D,2D,5D,10D]

  ;varr = [range(0.01,0.09,10),range(0.1,0.9,10),range(1,10,10)]
  varr = [0.1,0.25,0.5,0.75,1,5]
  covarr = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
  ;varr = 0.1D
  ;covarr = 0.9

  total_data = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n,2,n_runs,/double)


  for j=0, n_elements(sig_x)-1 do begin
     xscat = Nprior(n, 0, sig_x[j])
     yscat = Nprior(n, 0, sig_y[j])
     
     
     newx = x + xscat
     newy = y + yscat
     
     
     for k=0, n_elements(varr)-1 do begin
        for kk=0, n_elements(covarr)-1 do begin
           covar = make_array(2,2,/double)
           covar[0,0] = varr[k] ;x varrience
           covar[1,0] = varr[k] * covarr[kk] ;covarience x,y
           covar[0,1] = varr[k] * covarr[kk] ;covarience x,y
           covar[1,1] = varr[k]   ;y varrience
        
           
           data = make_array(2,n,/double)
           data[0,*] = ABS(newx)
           data[1,*] = ABS(newy)
           dist = twoD_pdf(data,bins=500,info=info)
           
           
           ran = make_array(2,2)
           ran[0,*] = [info.min1,info.max1]
           ran[1,*] = [info.min2,info.max2]
           
           master_dist = dist
           master_dist[*,*] = 0D
           
       
           
           for i=0, n-1 do begin
              
              rand = mrandomn(seed, covar, n_runs)
              rand = transpose(rand)
              randx = data[0,i] + rand[0,*]
              randy = data[1,i] + rand[1,*]
              
              
              
              
              new_data = make_array(2,n_runs,/double)
              new_data[0,*] = ABS(randx)
              new_data[1,*] = ABS(randy)
              
              total_data[j,k,kk,i,0,*] = new_data[0,*]
              total_data[j,k,kk,i,1,*] = new_data[1,*]
              
              new_dist = twoD_pdf(new_data,bins=500,range=ran,info=info)
              
              maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
              scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
              new_dist = scaledDensity
              
              master_dist = master_dist + new_dist
           endfor


        endfor
     endfor
  endfor

 ; maxDensity = Ceil(Max(master_dist)/1e2) * 1e2
 ; scaledDensity = BytScl(master_dist, Min=0, Max=maxDensity)
  

  ;cgImage, master_dist, XRange=[info.min1,info.max1], YRange=[info.min2,info.max2], /Axes, Palette=palette, $
  ;         XTitle='x', YTitle='y',$
  ;         charsize=0.9

  ;pears = CORRELATE(total_data[*,0,*], total_data[*,1,*])  
  ;Gau_fit = make_array(n_runs,2,/double)
  ;Gau_pears = make_array(n_runs,/double)
  ;Gau_spear = make_array(n_runs,/double)
  ;Gau_kendal = make_array(n_runs,/double)

  

  lin_fit = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n_runs,2,/double)
  pears = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n_runs,/double)
  spear = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n_runs,/double)
  kendal = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n_runs,/double)

  cordata_x = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n,/double)
  cordata_y = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr),n,/double)
              

  for i=0, n_elements(sig_x)-1 do begin
     for j=0, n_elements(varr)-1 do begin
        for k=0, n_elements(covarr)-1 do begin
           for l=0l, n_runs-1 do begin
              

              
              cordata_x[i,j,k,*] = total_data[i,j,k,*,0,l]
              cordata_y[i,j,k,*] = total_data[i,j,k,*,1,l]
             
              
              result = LINFIT(cordata_x[i,j,k,*], cordata_y[i,j,k,*])
              lin_fit[i,j,k,l,*] = result
     
              pears[i,j,k,l] = CORRELATE(cordata_x[i,j,k,*], cordata_y[i,j,k,*]) 
              result = R_CORRELATE(cordata_x[i,j,k,*], cordata_y[i,j,k,*])
              spear[i,j,k,l] = result[0]
              result = R_CORRELATE(cordata_x[i,j,k,*], cordata_y[i,j,k,*], /KENDALL)
              kendal[i,j,k,l] = result[0]
              
           endfor
        endfor
     endfor
  endfor

  cgLoadCT, 34
  
  position = [[0.1,0.76,0.3,0.96],[0.35,0.76,0.55,0.96],[0.1,0.54,0.3,0.74],[0.35,0.54,0.55,0.74],[0.1,0.32,0.3,0.52],[0.35,0.32,0.55,0.52],[0.1,0.1,0.3,0.3],[0.35,0.1,0.55,0.3]]

  
  ;bb = (lin_fit[*,*,*,*,1])
  
  
  mean_a = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr))
  mean_b = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr))
  mean_pears = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr))
  mean_spear = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr))
  mean_kendal = make_array(n_elements(sig_x),n_elements(varr),n_elements(covarr))
  
  
  for i=0, n_elements(sig_x)-1 do begin
     for j=0, n_elements(varr)-1 do begin
        for k=0, n_elements(covarr)-1 do begin
           mean_a[i,j,k] = mean(lin_fit[i,j,k,*,0])
           mean_b[i,j,k] = mean(lin_fit[i,j,k,*,1])
           mean_pears[i,j,k] = mean(pears[i,j,k,*])
           mean_spear[i,j,k] = mean(spear[i,j,k,*])          
        endfor
     endfor
  endfor
  
  cgDisplay , 700.,700.,/aspect
  cgPS_Open,strcompress('./plots/SDSS_corrtest.ps')

  color = plotcolor(mean_a)

  cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,0],/nodata,charsize=1
  
  for i=0, n_elements(sig_x)-1 do begin     
     cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,i],charsize=1,/noerase
    
     for j=0, n_elements(varr)-1 do begin     
        for k=0, n_elements(covarr)-1 do begin
           cgplot, varr[j], covarr[k], color=color[i,j,k],/overplot,psym=16
        endfor
     endfor
     
  endfor
  
  cgCOLORBAR, range=[min(mean_a),max(mean_a)],/vertical,POSITION=[0.6, 0.32, 0.62, 0.74],color='black',/right

  cghistoplot, lin_fit[*,*,*,*,0],/outline 


  color = plotcolor(mean_b)

  cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,0],/nodata,charsize=1
  
  for i=0, n_elements(sig_x)-1 do begin     
     cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,i],charsize=1,/noerase
    
     for j=0, n_elements(varr)-1 do begin     
        for k=0, n_elements(covarr)-1 do begin
           cgplot, varr[j], covarr[k], color=color[i,j,k],/overplot,psym=16
        endfor
     endfor
     
  endfor
  
  cgCOLORBAR, range=[min(mean_b),max(mean_b)],/vertical,POSITION=[0.6, 0.32, 0.62, 0.74],color='black',/right

  cghistoplot, lin_fit[*,*,*,*,1],/outline


 color = plotcolor(mean_pears)

  cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,0],/nodata,charsize=1
  
  for i=0, n_elements(sig_x)-1 do begin     
     cgplot, varr[0],covarr[0],psym=16,yr=[0.09,1],xr=[0.05,10],/xlog,/ylog,position=position[*,i],charsize=1,/noerase
    
     for j=0, n_elements(varr)-1 do begin     
        for k=0, n_elements(covarr)-1 do begin
           cgplot, varr[j], covarr[k], color=color[i,j,k],/overplot,psym=16
        endfor
     endfor
     
  endfor
  
  cgCOLORBAR, range=[min(mean_pears),max(mean_pears)],/vertical,POSITION=[0.6, 0.32, 0.62, 0.74],color='black',/right

  cghistoplot, pears[*,*,*,*],/outline

  cgPS_close
  cgPS2PDF, strcompress('./plots/SDSS_corrtest.ps')

  stop
end

pro danm_corr
  
  n = 100
  n_runs = 1000
  

  x = Uprior(n, 0D, 10D, seed)
  y = Uprior(n, 0D, 10D, seed)
  
  ;varr = [range(0.01,0.09,10),range(0.1,0.9,10),range(1,10,10)]
  varr = [0.01,0.3]
  covarr = [0.1,0.9]
  ;varr = 0.1D
  ;covarr = 0.1

  total_data = make_array(n_elements(varr),n_elements(covarr),n,2,n_runs,/double)

  cgLoadCT, 33
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]
   
  cgDisplay , 700.,700.,/aspect
  cgPS_Open,strcompress('./plots/corrtest_again.ps')
     
  position = [[0.1,0.1,0.5,0.5],[0.6,0.1,0.95,0.5],[0.1,0.6,0.5,0.95],[0.6,0.6,0.95,0.95]]

  for k=0, n_elements(varr)-1 do begin
     for kk=0, n_elements(covarr)-1 do begin
        covar = make_array(2,2,/double)
        covar[0,0] = varr[k]                 ;x varrience
        covar[1,0] = varr[k] * covarr[kk]    ;covarience x,y
        covar[0,1] = varr[k] * covarr[kk]    ;covarience x,y
        covar[1,1] = varr[k]                 ;y varrience
        
           
        data = make_array(2,n,/double)
        data[0,*] = ABS(x)
        data[1,*] = ABS(y)
        dist = twoD_pdf(data,bins=250,info=info)
        
           
        ran = make_array(2,2)
        ran[0,*] = [info.min1,info.max1]
        ran[1,*] = [info.min2,info.max2]
        
        master_dist = dist
        master_dist[*,*] = 0D
        
       
           
        for i=0, n-1 do begin
           
           rand = mrandomn(seed, covar, n_runs)
           rand = transpose(rand)
           randx = data[0,i] + rand[0,*]
           randy = data[1,i] + rand[1,*]
           
           
              
              
           new_data = make_array(2,n_runs,/double)
           new_data[0,*] = ABS(randx)
           new_data[1,*] = ABS(randy)
           
           total_data[k,kk,i,0,*] = new_data[0,*]
           total_data[k,kk,i,1,*] = new_data[1,*]
           
           new_dist = twoD_pdf(new_data,bins=500,range=ran,info=info)
           
           maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
           scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
           new_dist = scaledDensity
           
           master_dist = master_dist + new_dist
        endfor

        if k eq 0 and kk eq 0 then cgImage, master_dist, XRange=[info.min1,info.max1], YRange=[info.min2,info.max2], /Axes, Palette=palette, $
           XTitle='x', YTitle='y', $
           Position=position[*,k],$
           charsize=0.9
        if k gt 0 and kk gt 0 then cgImage, master_dist, XRange=[info.min1,info.max1], YRange=[info.min2,info.max2], /Axes, Palette=palette, $
           XTitle='x', YTitle='y', $
           Position=position[*,k],$
           charsize=0.9,/noerase
        
        
     endfor
  endfor
  


  

  lin_fit = make_array(n_elements(varr),n_elements(covarr),n_runs,2,/double)
  pears = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  spear = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  s_prob = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  kendal = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)

  cordata_x = make_array(n_elements(varr),n_elements(covarr),n,/double)
  cordata_y = make_array(n_elements(varr),n_elements(covarr),n,/double)
              
  
  
  for j=0, n_elements(varr)-1 do begin
     for k=0, n_elements(covarr)-1 do begin
        for l=0l, n_runs-1 do begin
           
           
           
           cordata_x[j,k,*] = total_data[j,k,*,0,l]
           cordata_y[j,k,*] = total_data[j,k,*,1,l]
           
           
           result = LINFIT(cordata_x[j,k,*], cordata_y[j,k,*])
           lin_fit[j,k,l,*] = result
           
           pears[j,k,l] = CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*]) 
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*])
           spear[j,k,l] = result[0]
           s_prob[j,k,l] = result[1]
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*], /KENDALL)
           kendal[j,k,l] = result[0]
           
        endfor
     endfor
  endfor
  
  x = range(0,10,100)
  y1 = x
  y2 = x/0.05D

  cgplot, total_data[0,0,*,0,*],total_data[0,0,*,1,*],psym=3,Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9,xt='x',yt='y',color='green' 
  cgplot, total_data[0,1,*,0,*],total_data[0,1,*,1,*],psym=3,/overplot,color='purple'
  
  
  
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cgplot, total_data[1,0,*,0,*],total_data[1,0,*,1,*],psym=3,Position=[0.55,0.5,0.95,0.9],$
           charsize=0.9,xt='x',yt='y',color='orange',/noerase
  cgplot, total_data[1,1,*,0,*],total_data[1,1,*,1,*],psym=3,/overplot,color='yellow'
  
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cghistoplot, lin_fit[0,0,*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='a ',color='green',xrange=[min(lin_fit[*,*,*,0]), max(lin_fit[*,*,*,0])]
  cghistoplot, lin_fit[0,1,*,0], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,0], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,0], /outline,/oplot,color='yellow'

    
    
  cghistoplot, lin_fit[0,0,*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='b ',color='green',xrange=[min(lin_fit[*,*,*,1]), max(lin_fit[*,*,*,1])]
  cghistoplot, lin_fit[0,1,*,1], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,1], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,1], /outline,/oplot,color='yellow'
 
  cghistoplot, spear[0,0,*],/outline, position=[0.1,0.1,0.25,0.25],xrange=[min(spear), max(spear)],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$',color='green'
  cghistoplot, spear[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, spear[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, spear[1,1,*], /outline,/oplot,color='yellow'
  
  cghistoplot, s_prob,/outline, position=[0.35,0.1,0.5,0.25],xrange=[min(s_prob), max(s_prob)],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$ probability',color='green'
  cghistoplot, s_prob[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, s_prob[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, s_prob[1,1,*], /outline,/oplot,color='yellow'
  
  
;**********************************************************************************
;Same again but now with the 0.05 < x/y < 1 prior

  
  x = Uprior(150, 0D, 6D, seed)
  ;x_y = Uprior(n, 0.05D, 1.D, seed)
  ;y = x / x_y
  y = Uprior(150, 0D, 10D, seed)
  xy = x/y
  index = where((xy le 1D) and (xy ge 0.05))
  x = x[index]
  y = y[index]
  
  n = n_elements(x)
  ;varr = [range(0.01,0.09,10),range(0.1,0.9,10),range(1,10,10)]
  ;varr = [0.1,0.25,0.5,0.75,1,5]
  ;covarr = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
  ;varr = 0.1D
  ;covarr = 0.9

  total_data = make_array(n_elements(varr),n_elements(covarr),n,2,n_runs,/double)


  xx = range(0,10,100)
  yy1 = xx
  yy2 = xx/0.05D
     
  for k=0, n_elements(varr)-1 do begin
     for kk=0, n_elements(covarr)-1 do begin
        covar = make_array(2,2,/double)
        covar[0,0] = varr[k]                 ;x varrience
        covar[1,0] = varr[k] * covarr[kk]    ;covarience x,y
        covar[0,1] = varr[k] * covarr[kk]    ;covarience x,y
        covar[1,1] = varr[k]                 ;y varrience
        
           
        data = make_array(2,n,/double)
        data[0,*] = ABS(x)
        data[1,*] = ABS(y)
        dist = twoD_pdf(data,bins=250,info=info)
        
           
        ran = make_array(2,2)
        ran[0,*] = [info.min1,info.max1]
        ran[1,*] = [info.min2,info.max2]
        
        master_dist = dist
        master_dist[*,*] = 0D
        
       
           
        for i=0, n-1 do begin
           
           rand = mrandomn(seed, covar, n_runs)
           rand = transpose(rand)
           randx = data[0,i] + rand[0,*]
           randy = data[1,i] + rand[1,*]
           
           
              
              
           new_data = make_array(2,n_runs,/double)
           new_data[0,*] = ABS(randx)
           new_data[1,*] = ABS(randy)
           
           total_data[k,kk,i,0,*] = new_data[0,*]
           total_data[k,kk,i,1,*] = new_data[1,*]
           
           new_dist = twoD_pdf(new_data,bins=500,range=ran,info=info)
           
           maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
           scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
           new_dist = scaledDensity
           
           master_dist = master_dist + new_dist
        endfor
        
        if k eq 0 and kk gt 0 then cgImage, master_dist, XRange=[info.min1,info.max1], YRange=[info.min2,info.max2], /Axes, Palette=palette, $
           XTitle='x', YTitle='y', $
           Position=position[*,k],$
           charsize=0.9
        if k gt 0 and kk gt 0 then cgImage, master_dist, XRange=[info.min1,info.max1], YRange=[info.min2,info.max2], /Axes, Palette=palette, $
           XTitle='x', YTitle='y', $
           Position=position[*,k],$
           charsize=0.9,/noerase
        cgplot, xx,yy1, color='red',/overplot
        cgplot, xx, yy2, color='red',/overplot
        
     endfor
  endfor
  


  

  lin_fit = make_array(n_elements(varr),n_elements(covarr),n_runs,2,/double)
  pears = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  spear = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  kendal = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  s_prob = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  cordata_x = make_array(n_elements(varr),n_elements(covarr),n,/double)
  cordata_y = make_array(n_elements(varr),n_elements(covarr),n,/double)
              

  
  for j=0, n_elements(varr)-1 do begin
     for k=0, n_elements(covarr)-1 do begin
        for l=0l, n_runs-1 do begin
           
           
           
           cordata_x[j,k,*] = total_data[j,k,*,0,l]
           cordata_y[j,k,*] = total_data[j,k,*,1,l]
           
           
           result = LINFIT(cordata_x[j,k,*], cordata_y[j,k,*])
           lin_fit[j,k,l,*] = result
           
           pears[j,k,l] = CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*]) 
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*])
           spear[j,k,l] = result[0]
           s_prob[j,k,l] = result[1]
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*], /KENDALL)
           kendal[j,k,l] = result[0]
           
        endfor
     endfor
  endfor
  
  x = range(0,10,100)
  y1 = x
  y2 = x/0.05D

  cgplot, total_data[0,0,*,0,*],total_data[0,0,*,1,*],psym=3,Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9,xt='x',yt='y',color='green' 
  cgplot, total_data[0,1,*,0,*],total_data[0,1,*,1,*],psym=3,/overplot,color='purple'
  cgplot, x,y1, color='red',/overplot
  cgplot, x, y2, color='red',/overplot
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cgplot, total_data[1,0,*,0,*],total_data[1,0,*,1,*],psym=3,Position=[0.55,0.5,0.95,0.9],$
           charsize=0.9,xt='x',yt='y',color='orange',/noerase
  cgplot, total_data[1,1,*,0,*],total_data[1,1,*,1,*],psym=3,/overplot,color='yellow'
  cgplot, x,y1, color='red',/overplot
  cgplot, x, y2, color='red',/overplot
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cghistoplot, lin_fit[0,0,*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='a ',color='green',xrange=[min(lin_fit[*,*,*,0]), max(lin_fit[*,*,*,0])]
  cghistoplot, lin_fit[0,1,*,0], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,0], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,0], /outline,/oplot,color='yellow'

    
    
  cghistoplot, lin_fit[0,0,*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='b ',color='green',xrange=[min(lin_fit[*,*,*,1]), max(lin_fit[*,*,*,1])]
  cghistoplot, lin_fit[0,1,*,1], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,1], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,1], /outline,/oplot,color='yellow'
 
  cghistoplot, spear[0,0,*],/outline, position=[0.1,0.1,0.25,0.25],xrange=[min(spear), max(spear)],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$',color='green'
  cghistoplot, spear[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, spear[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, spear[1,1,*], /outline,/oplot,color='yellow'
  
  cghistoplot, s_prob,/outline, position=[0.35,0.1,0.5,0.25],xrange=[min(s_prob), 1e-11],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$ probability',color='green'
  cghistoplot, s_prob[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, s_prob[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, s_prob[1,1,*], /outline,/oplot,color='yellow'
  

cgPS_close
  cgPS2PDF, strcompress('./plots/corrtest_again.ps')

  stop
;**********************************************************************************
;Same again but now with the 0.001< x/y < 2. prior

  n=100
  x = Uprior(n, 0D, 10D, seed)
  x_y = Uprior(n, 0.001D, 2D, seed)
  y = x / x_y
  
  ;varr = [range(0.01,0.09,10),range(0.1,0.9,10),range(1,10,10)]
  ;varr = [0.1,0.25,0.5,0.75,1,5]
  ;covarr = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
  ;varr = 0.1D
  ;covarr = 0.9

  total_data = make_array(n_elements(varr),n_elements(covarr),n,2,n_runs,/double)


   
     
  for k=0, n_elements(varr)-1 do begin
     for kk=0, n_elements(covarr)-1 do begin
        covar = make_array(2,2,/double)
        covar[0,0] = varr[k]                 ;x varrience
        covar[1,0] = varr[k] * covarr[kk]    ;covarience x,y
        covar[0,1] = varr[k] * covarr[kk]    ;covarience x,y
        covar[1,1] = varr[k]                 ;y varrience
        
           
        data = make_array(2,n,/double)
        data[0,*] = ABS(x)
        data[1,*] = ABS(y)
        dist = twoD_pdf(data,bins=500,info=info)
        
           
        ran = make_array(2,2)
        ran[0,*] = [info.min1,info.max1]
        ran[1,*] = [info.min2,info.max2]
        
        master_dist = dist
        master_dist[*,*] = 0D
        
       
           
        for i=0, n-1 do begin
           
           rand = mrandomn(seed, covar, n_runs)
           rand = transpose(rand)
           randx = data[0,i] + rand[0,*]
           randy = data[1,i] + rand[1,*]
           
           
              
              
           new_data = make_array(2,n_runs,/double)
           new_data[0,*] = ABS(randx)
           new_data[1,*] = ABS(randy)
           
           total_data[k,kk,i,0,*] = new_data[0,*]
           total_data[k,kk,i,1,*] = new_data[1,*]
           
           new_dist = twoD_pdf(new_data,bins=500,range=ran,info=info)
           
           maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
           scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
           new_dist = scaledDensity
           
           master_dist = master_dist + new_dist
        endfor
        
        
     endfor
  endfor
  


  

  lin_fit = make_array(n_elements(varr),n_elements(covarr),n_runs,2,/double)
  pears = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  spear = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  kendal = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  s_prob = make_array(n_elements(varr),n_elements(covarr),n_runs,/double)
  cordata_x = make_array(n_elements(varr),n_elements(covarr),n,/double)
  cordata_y = make_array(n_elements(varr),n_elements(covarr),n,/double)
              

  
  for j=0, n_elements(varr)-1 do begin
     for k=0, n_elements(covarr)-1 do begin
        for l=0l, n_runs-1 do begin
           
           
           
           cordata_x[j,k,*] = total_data[j,k,*,0,l]
           cordata_y[j,k,*] = total_data[j,k,*,1,l]
           
           
           result = LINFIT(cordata_x[j,k,*], cordata_y[j,k,*])
           lin_fit[j,k,l,*] = result
           
           pears[j,k,l] = CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*]) 
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*])
           spear[j,k,l] = result[0]
           s_prob[j,k,l] = result[1]
           result = R_CORRELATE(cordata_x[j,k,*], cordata_y[j,k,*], /KENDALL)
           kendal[j,k,l] = result[0]
           
        endfor
     endfor
  endfor
  
  x = range(0,10,100)
  y1 = x/2D
  y2 = x/0.001D

  cgplot, total_data[0,0,*,0,*],total_data[0,0,*,1,*],psym=3,Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9,xt='x',yt='y',color='green' 
  cgplot, total_data[0,1,*,0,*],total_data[0,1,*,1,*],psym=3,/overplot,color='purple'
  cgplot, x,y1, color='red',/overplot
  cgplot, x, y2, color='red',/overplot
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cgplot, total_data[1,0,*,0,*],total_data[1,0,*,1,*],psym=3,Position=[0.55,0.5,0.95,0.9],$
           charsize=0.9,xt='x',yt='y',color='orange',/noerase
  cgplot, total_data[1,1,*,0,*],total_data[1,1,*,1,*],psym=3,/overplot,color='yellow'
  cgplot, x,y1, color='red',/overplot
  cgplot, x, y2, color='red',/overplot
  y = median(lin_fit[*,*,*,0]) + (median(lin_fit[*,*,*,1]) * x)
  cgplot, (x),(y),color='blue',/overplot
  
  
  cghistoplot, lin_fit[0,0,*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='a ',color='green',xrange=[min(lin_fit[*,*,*,0]), max(lin_fit[*,*,*,0])]
  cghistoplot, lin_fit[0,1,*,0], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,0], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,0], /outline,/oplot,color='yellow'

    
    
  cghistoplot, lin_fit[0,0,*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='b ',color='green',xrange=[min(lin_fit[*,*,*,1]), max(lin_fit[*,*,*,1])]
  cghistoplot, lin_fit[0,1,*,1], /outline,/oplot,color='purple'
  cghistoplot, lin_fit[1,0,*,1], /outline,/oplot,color='orange'
  cghistoplot, lin_fit[1,1,*,1], /outline,/oplot,color='yellow'
 
  cghistoplot, spear[0,0,*],/outline, position=[0.1,0.1,0.25,0.25],xrange=[min(spear), max(spear)],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$',color='green'
  cghistoplot, spear[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, spear[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, spear[1,1,*], /outline,/oplot,color='yellow'
  
  cghistoplot, s_prob,/outline, position=[0.35,0.1,0.5,0.25],xrange=[min(s_prob), 0.0002],/noerase,charsize=.3,ytitle='',ytickformat='(A1)',xtitle='$\rho$ probability',color='green'
  cghistoplot, s_prob[0,1,*],/outline,/oplot,color='purple'
  cghistoplot, s_prob[1,0,*], /outline,/oplot,color='orange'
  cghistoplot, s_prob[1,1,*], /outline,/oplot,color='yellow'


  cgPS_close
  cgPS2PDF, strcompress('./plots/corrtest_again.ps')

  stop
end



pro newcorr
  
   ;Read in Gadotti '09 results
  readcol, 'gad09.txt',name, xc, yc, sky, seeing_Gadotti,seeing_lower, seeing_higher, seeing_mean, seeing_1_G2, seeing_1_G2,seeing_Effwidth ,gain, Dark_current, zero_point, airmass, zcal,mue, Ie_gad, re_gad, n_gad, mu0, I0_gad, h_gad, z,format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',SKIPLINE=1

  ;Read in GASP2D results 
  readcol, '2Dfit_g09.txt',GALAXY,Ie_g2d,Re_g2d, n_g2d,ba_g2d,PA_g2d, BT_g2d,I0_g2d, h_g2d,ba_g2d,PA_g2d,DT_g2d,Xcenter_g2d,Ycenter_g2d,Chi_g2d,format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F',skipline=1
 
  RESOLVE_ROUTINE, 'cosmology_routines'
  red,omega0=0.3,omegalambda=0.7,h100=0.7D

   
  
 ;-----------------------
  dist=DANGULAR(z)
  physcale=dist/206265D
  scale=0.3961D
  band='i'
 ;----------------------
  
 


  ;Read in posterior distribution   
  path = './outputs/SDSS/Gad09/iband/Disc/'
  results = file_search(path+'SDSS*', count=nima)
  ngals = n_elements(results)
  
  res0 = readfits(results[0],/SILENT)
  ressize = size(res0, /DIMENSION)
  npars = ressize[0] 
  
  n_runs = 10000 
  
  meds = make_array(ngals,npars,/DOUBLE)
  q = make_array(ngals,npars,6,/DOUBLE)
  low_3s = make_array(ngals,npars,/DOUBLE)
  low_2s = make_array(ngals,npars,/DOUBLE)
  low_1s = make_array(ngals,npars,/DOUBLE)
  
  high_1s = make_array(ngals,npars,/DOUBLE)
  high_2s = make_array(ngals,npars,/DOUBLE)
  high_3s = make_array(ngals,npars,/DOUBLE)

  err = make_array(ngals,npars,/DOUBLE)

  
  
  ie = make_array(ngals,/DOUBLE)
  re = make_array(ngals,/DOUBLE)
  n = make_array(ngals,/DOUBLE)
  serE = make_array(ngals,/DOUBLE)
  serPA = make_array(ngals,/DOUBLE)
  i0 = make_array(ngals,/DOUBLE)
  h = make_array(ngals,/DOUBLE)
  expE = make_array(ngals,/DOUBLE)
  expPA = make_array(ngals,/DOUBLE)
  x0 = make_array(ngals,/DOUBLE)
  y0 = make_array(ngals,/DOUBLE)
  

  
  newRe = make_array(n_elements(results),n_runs,/DOUBLE)
  newh = make_array(n_elements(results),n_runs,/DOUBLE)
  
  
  for i=0l, n_elements(results)-1 do begin
     ;Read in distributions 
        
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     

        ;data = 10D^res[[3,10],*]
     ;Random points from PDF
     ir = floor(Uprior(n_runs, 0, n_elements(Ie), seed))
     
     newRe[i,*] = re[ir]
     newh[i,*] = h[ir]
     
     ;Inverse transforming sampling 
     ;newsample = SAMPLE_PDF_ND(data,nsample=1000)
        
       
     for jj=0, npars-1 do begin
        quartile = fractile(res[jj,*], [0.0015,0.025,0.15865,0.84135,0.97635,0.9985])
        q[i,jj,*] = quartile
        meds[i,jj] = median(res[jj,*],/EVEN,/DOUBLE)  
     endfor
     
        
  endfor
  
  index = [3,10]
  meds[*,index] = 10D^meds[*,index]
  q[*,index,*] = 10D^(q[*,index,*])

  ;Converting Re and h into physical units
  for i=0l, n_elements(results)-1 do begin 
     newre[i,*] = (newre[i,*]*scale)*physcale[i]/1000.
  endfor
  rePhys_gad = (re_gad*scale)*physcale/1000.
  rePhys_g2d = (re_g2d*scale)*physcale/1000.
  rePhys = (meds[*,3]*scale)*physcale/1000.
  rePhys_q1 = (q[*,3,2]*scale)*physcale/1000.
  rePhys_q2 = (q[*,3,3]*scale)*physcale/1000.

  for i=0l, n_elements(results)-1 do begin
     newh[i,*] = (newh[i,*]*scale)*physcale[i]/1000.
  endfor
  hPhys_gad = (h_gad*scale)*physcale/1000.
  hPhys_g2d = (h_g2d*scale)*physcale/1000.
  hPhys = (meds[*,10]*scale)*physcale/1000.
  hPhys_q1 = (q[*,10,2]*scale)*physcale/1000.
  hPhys_q2 = (q[*,10,3]*scale)*physcale/1000.
  

cgLoadCT, 33
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  ;2D hist for sdss sample normalised 
  data = make_array(2,N_ELEMENTS(newre),/DOUBLE)
  data[0,*] = (newre[*,*])
  data[1,*] = (newh[*,*])
  density = twoD_pdf(data,bins=500,info=info)
  good_info = info
  master_density = density
  master_density[*,*] = 0D

  

  for i=0, n_elements(results)-1 do begin
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     ran = make_array(2,2)
     ran[0,*] = [good_info.min1,good_info.max1]
     ran[1,*] = [good_info.min2,good_info.max2]

     re = (re*scale)*physcale[i]/1000.
     h = (h*scale)*physcale[i]/1000.

     new_data = make_array(2,N_ELEMENTS(re),/DOUBLE)
     new_data[0,*] = (re)
     new_data[1,*] = (h)
     new_dist = twoD_pdf(new_data,bins=500, range=ran, info=newinfo)
     ;new_dist = (new_dist/total(new_dist))* ((newinfo.bin1)*(newinfo.bin2))
     
     maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
     scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
     new_dist = scaledDensity
     master_density = master_density + new_dist
     
  endfor
  
  
;*****************************************************
;Adding mock galxies from gadotti 09 fitting results

  ;Read in the information about gaalxies

  path = './sample/SDSS/Mock_Gad09/iband/Disc/'
  readcol,path+'MOGAL_se_ivals.txt',NAME,z,serX0,serY0,mue,Re_gad,n,serE,serPA,serM,ser2tot,serSFH,serSF,expX0,expY0,mu0,h_gad,expE,expPA,expM,exp2tot,expSFH,expSF,zz,seeing,SNR,Nsam,Seein,M_ser,m_ser,M_exp,m_exp, FORMAT='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F',skipline=1
  
  ;-----------------------
  dist=DANGULAR(z)
  physcale=dist/206265D
  scale=0.3961D
  band='i'
 ;----------------------
  
  
  
  ;Read in posterior distribution   
  path = './outputs/SDSS/Mock_Gad09/iband/Disc/2component/'
  results = file_search(path+'MOGAL*', count=nima)
  ngals = n_elements(results)
  
  res0 = readfits(results[0],/SILENT)
  ressize = size(res0, /DIMENSION)
  npars = ressize[0] 
  
  n_runs = 10000 
  
  meds = make_array(ngals,npars,/DOUBLE)
  q = make_array(ngals,npars,6,/DOUBLE)
  low_3s = make_array(ngals,npars,/DOUBLE)
  low_2s = make_array(ngals,npars,/DOUBLE)
  low_1s = make_array(ngals,npars,/DOUBLE)
  
  high_1s = make_array(ngals,npars,/DOUBLE)
  high_2s = make_array(ngals,npars,/DOUBLE)
  high_3s = make_array(ngals,npars,/DOUBLE)

  err = make_array(ngals,npars,/DOUBLE)

  
  
  ie = make_array(ngals,/DOUBLE)
  re = make_array(ngals,/DOUBLE)
  n = make_array(ngals,/DOUBLE)
  serE = make_array(ngals,/DOUBLE)
  serPA = make_array(ngals,/DOUBLE)
  i0 = make_array(ngals,/DOUBLE)
  h = make_array(ngals,/DOUBLE)
  expE = make_array(ngals,/DOUBLE)
  expPA = make_array(ngals,/DOUBLE)
  x0 = make_array(ngals,/DOUBLE)
  y0 = make_array(ngals,/DOUBLE)
  

  
  newRe = make_array(n_elements(results),n_runs,/DOUBLE)
  newh = make_array(n_elements(results),n_runs,/DOUBLE)
  
  
  for i=0l, n_elements(results)-1 do begin
     ;Read in distributions 
        
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     

        ;data = 10D^res[[3,10],*]
     ;Random points from PDF
     ir = floor(Uprior(n_runs, 0, n_elements(Ie), seed))
     
     newRe[i,*] = re[ir]
     newh[i,*] = h[ir]
     
     ;Inverse transforming sampling 
     ;newsample = SAMPLE_PDF_ND(data,nsample=1000)
        
       
     for jj=0, npars-1 do begin
        quartile = fractile(res[jj,*], [0.0015,0.025,0.15865,0.84135,0.97635,0.9985])
        q[i,jj,*] = quartile
        meds[i,jj] = median(res[jj,*],/EVEN,/DOUBLE)  
     endfor
     
        
  endfor
  
  index = [3,10]
  meds[*,index] = 10D^meds[*,index]
  q[*,index,*] = 10D^(q[*,index,*])

  ;Converting Re and h into physical units
  for i=0l, n_elements(results)-1 do begin 
     newre[i,*] = (newre[i,*]*scale)*physcale[i]/1000.
  endfor
  rePhys = (meds[*,3]*scale)*physcale/1000.
  rePhys_q1 = (q[*,3,2]*scale)*physcale/1000.
  rePhys_q2 = (q[*,3,3]*scale)*physcale/1000.

  for i=0l, n_elements(results)-1 do begin
     newh[i,*] = (newh[i,*]*scale)*physcale[i]/1000.
  endfor
  hPhys = (meds[*,10]*scale)*physcale/1000.
  hPhys_q1 = (q[*,10,2]*scale)*physcale/1000.
  hPhys_q2 = (q[*,10,3]*scale)*physcale/1000.

  rePhys_gad = (re_gad*scale)*physcale/1000.
  hPhys_gad = (h_gad*scale)*physcale/1000.
 

  ;2D hist for sdss sample normalised 
  data = make_array(2,N_ELEMENTS(newre),/DOUBLE)
  data[0,*] = (newre[*,*])
  data[1,*] = (newh[*,*])
  density = twoD_pdf(data,bins=500,info=info)
  good_info = info
  master_density = density
  master_density[*,*] = 0D

  

  for i=0, n_elements(results)-1 do begin
     res = readfits(results[i],/SILENT)
     Ie = reform(10D^res[2,*])
     Re = reform(10D^res[3,*])
     n= reform(res[4,*])
     serE = reform(res[5,*])        
     serPA = reform(res[6,*])
     I0= reform(10D^res[9,*])      
     h= reform(10D^res[10,*])
     expE = reform(res[11,*])
     expPA = reform(res[12,*])
     x0 = reform(res[13,*])
     y0 = reform(res[14,*])
     ran = make_array(2,2)
     ran[0,*] = [good_info.min1,good_info.max1]
     ran[1,*] = [good_info.min2,good_info.max2]

     re = (re*scale)*physcale[i]/1000.
     h = (h*scale)*physcale[i]/1000.

     new_data = make_array(2,N_ELEMENTS(re),/DOUBLE)
     new_data[0,*] = (re)
     new_data[1,*] = (h)
     new_dist = twoD_pdf(new_data,bins=500, range=ran, info=newinfo)
     ;new_dist = (new_dist/total(new_dist))* ((newinfo.bin1)*(newinfo.bin2))
     
     maxDensity = Ceil(Max(new_dist)/1e2) * 1e2
     scaledDensity = BytScl(new_dist, Min=0, Max=maxDensity)
     new_dist = scaledDensity
     master_density = master_density + new_dist
     
  endfor
stop

;**************************************
  n_meds = meds[*,4] 
  nind = where(n_meds ge 2)

  lin_fit = make_array(n_runs,2,/double)
  pears = make_array(n_runs,/double)
  spear = make_array(n_runs,/double)
  kendal = make_array(n_runs,/double)

  for j=0l, n_runs-1 do begin
     
     result = LINFIT(newre[*,j], newh[*,j])
     lin_fit[j,*] = result
     pears[j] = CORRELATE(newre[*,j], newh[*,j])
     result = R_CORRELATE(newre[*,j], newh[*,j])
     spear[j] = result[0]
     result = R_CORRELATE(newre[*,j], newh[*,j], /KENDALL)
     kendal[j] = result[0]

  endfor

  

  


  
  
  cgDisplay , 700.,700.,/aspect
  cgPS_Open,strcompress('./plots/SDSS_newcorr.ps')

  ;cgplot, newre, newh, psym=3, /xlog, /ylog,xrange=[0.1,5],yrange=[0.5,20], position=[0.1,0.5,0.5,0.9],xt='R$\downe$ (/Kpc)', yt='h (/Kpc)',title='$\Phi$-PDF',charsize=0.9
  
  cgImage, master_Density, XRange=[good_info.min1,good_info.max1], YRange=[good_info.min2,good_info.max2], /Axes, Palette=palette, $
           XTitle='R$\downe$ (/Kpc)', YTitle='h (/Kpc)',title='$\Phi$-PDF', $
           Position=[0.1,0.5,0.5,0.9],$
           charsize=0.9
  
  ;for i=0, n_runs-1 do begin
     x = range(min(newre[*,*]), max(newre[*,*]), 1000)
     y = median(lin_fit[*,0]) + (median(lin_fit[*,1]) * x)
     cgplot, (x),(y),color='red',/overplot
  ;endfor


  
  cghistoplot, lin_fit[*,0], /outline, position=[0.1,0.3,0.25,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='a ($\Phi$-PDF)'

  cghistoplot, lin_fit[*,1], /outline, position=[0.35,0.3,0.5,0.45],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='b ($\Phi$-PDF)'


  cghistoplot, pears,/outline, position=[0.225,0.1,0.375,0.25],xrange=[min([min(pears),min(spear),min(kendal)]), max([max(pears),max(spear),max(kendal)])],/noerase,charsize=.5,ytitle='',ytickformat='(A1)',xtitle='r,$\rho$,$\tau$ ($\Phi$-PDF)'
  cghistoplot, spear,/outline,/oplot,color='blue'
  cghistoplot, kendal,/outline,/oplot,color='dark green'


  
  cgPS_close
  cgPS2PDF, strcompress('./plots/SDSS_newcorr.ps')


  
  
  
  
  stop
end
