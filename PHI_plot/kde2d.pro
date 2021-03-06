;+
; NAME:
;    KDE
;
; PURPOSE:
;    Estimate the probability density underlying a set of discrete
;    samples (measurements) using the kernel density estimator method.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    d = kde(x, t)
;
; INPUTS:
;    x: discrete samples of the desired distribution.
;
;    t: values for which the probability density is required.
;
; KEYWORD PARAMETERS:
;    scale: smoothing factor, also called the bandwidth
;        used to compute the kernel density estimate
;
;    weight: weighting for sampled points.
;
; KEYWORD FLAGS:
;    By default, KDE uses the Epanechnikov kernel to compute
;    the kernel density estimate, this can be overridden by
;    setting one of the following flags:
;    GAUSSIAN: use Gaussian kernel
;    TRIANGULAR: use triangular kernel
;    BIWEIGHT: use biweight kernel
;
; OUTPUTS:
;    d: probability density estimated at each value specified by t
;

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

;+
; NAME:
;    QUANTILE_THRESHOLD
;
; PURPOSE:
;    This function estimates the thresholds in the input data array 
;    corresponding to the specified quantiles.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    result = quantile_threshold( data, quantile )
;
; INPUTS:
;    DATA:  A numerical array of data, of size N_DATA.
;    QUANTILE:  A floating point vector listing quantiles to find in DATA.  
;        Values must be within [0,1].  Of length N_QUANTILE
;
; KEYWORD PARAMETERS:
;    PRESORTED:  If set then the function assumes that the values in DATA have 
;        already been sorted in ascending order, thus running more efficiently. 
;        The default is to assume that they have not been sorted and thus to 
;        spend some time sorting them.
;
; OUTPUTS:
;    RESULT:  A numerical array containing the N_QUANTILE threshold values 
;        corresponding to the N_QUANTILE quantiles in DATA as specified in 
;        QUANTILE.
;
; USES:
;    -
;
; PROCEDURE:
;    This function estimates the thresholds corresponding to the quantiles 
;    using method #5 as documented in R.
;
; EXAMPLE:
;    data = randomn( 2, 100 )
;    quantile = [ 0.05, 0.95 ]
;    result = quantile_threshold( data, quantile )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dstone@lbl.gov), 2012-05-31
;    Modified:  DAS, 2012-07-05 (Added PRESORTED keyword)
;-

FUNCTION QUANTILE_THRESHOLD, $
    DATA, QUANTILE, $
    PRESORTED=presorted_opt

;***********************************************************************
; Constants and checks

; Ensure valid quantile values
if ( min( quantile, max=temp ) lt 0. ) or ( temp gt 1. ) then stop

; Count the quantiles
n_quantile = n_elements( quantile )

;***********************************************************************
; Estimate the thresholds

; Initialise the output vector of thresholds
threshold = fltarr( n_quantile )

; Identify valid values in DATA
id_data = where( finite( data ) eq 1, n_id_data )
if n_id_data eq 0 then stop
; Sort good values in data
data_sort = data[id_data]
if not( keyword_set( presorted_opt ) ) then begin
  id = sort( data_sort )
  data_sort = data_sort[id]
endif

; Define the index of the lower values closest to the thresholds
index = quantile * n_id_data - 0.5
; Deal with unresolved tail quantiles
id = where( index lt 0., n_id )
if n_id gt 0 then threshold[id] = data_sort[0]
id = where( index gt n_id_data-1., n_id )
if n_id gt 0 then threshold[id] = data_sort[n_id_data-1]
id = where( ( index ge 0. ) and ( index le n_id_data - 1. ), n_id )
index = index[id]
index_floor = floor( index )
threshold[id] = data_sort[index_floor] $
    + ( index - index_floor ) $
    * ( data_sort[index_floor+1] - data_sort[index_floor] )

;***********************************************************************
; The end

return, threshold
END


FUNCTION bandwidth,x
  ;r = percentiles(x,value=[0.25,0.75])
  r = quantile_threshold(x,[0.25D,0.75D]) 
  h = (r[1] - r[0])/1.34D
  return, 4D * 1.06D * min(sqrt(VARIANCE(x)), h) * N_ELEMENTS(x) ^ (-1D/5D)
END

FUNCTION outer_sub, X, Y
  
  xsize = size(X, /DIMENSIONS)
  ysize = size(Y, /DIMENSIONS)
  
  
  z = make_array(ysize[0],xsize[0],/DOUBLE)
  
  for iy = 0l, ysize[0]-1 do begin
     for ix = 0l, xsize[0]-1 do z[iy,ix] = X[ix] - Y[iy]
  endfor
  
  
  return, z
END

FUNCTION kde2d, x, y, h = h, n = n, lims = lims

  x = double(x)
  y = double(y)
  if N_ELEMENTS(lims) eq 0 then lims = [min(x), max(x), min(y), max(y)]
  if N_ELEMENTS(n) eq 0 then n = 25
  
  nx = N_ELEMENTS(x)
  if (N_elements(y) ne nx) then message,"data vectors must be the same length"
  if ( ARRAY_EQUAL(finite(x), 0) || ARRAY_EQUAL(finite(y), 0) ) then message,"missing or infinite values in the data are not allowed"
  if(ARRAY_EQUAL(finite(lims), 0)) then message,"only finite values are allowed in 'lims'"
  
  
  
  gx = range(lims[0], lims[1], n)
  gy = range(lims[2], lims[3], n)
  if N_ELEMENTS(h) eq 0 then h = [bandwidth(x), bandwidth(y)] else h = REPLICATE(h,2)
  
  if ARRAY_EQUAL(h lt 0, 1) then message,"bandwidths must be strictly positive"
  h = h/4                       ; for S's bandwidth scale
  ax = outer_sub(gx, x)/h[0]
  ay = outer_sub(gy, y)/h[1]

  
  z = ( transpose(Gprior(ax)) # Gprior(ay) ) / (nx * h[0] * h[1])
  return, z
END

