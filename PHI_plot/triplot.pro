;data must be organized as [ndims,npoints]
function plotcolor,x

xnew = round(255*(x - min(x))/(max(x)-min(x)))

col = byte(xnew)
return,col

end
Function RNORM, N, AVG=AVG, SIG=SIG, seed=seed
  
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(SIG) then SIG = 1D
 
  R = make_array(N,/DOUBLE)
  for i=0, N-1 do begin
     R[i] = RANDOMN(seed)
  endfor
  
  Result = (R * SIG) + AVG
  return, result
END
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
Function Gprior, x, AVG=AVG, SIG=SIG, seed=seed
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(SIG) then SIG = 1D
  
  Px = (1D/(sqrt(2D*!DPI)*SIG)) * EXP(-0.5*(((x - AVG)/sig)^2D))
  
  return,Px
END

pro triplot,xs, bins=bins, range=range,weights=weights,color=color,$
            kdeon=kdeon, smooth1d=smooth1d,titles=titles,info=info,$
            maxins=maxins,minins=minins,boxplot=boxplot,posplots=posplots,$
            cross=cross,crxcolor=crxcolor,crxlable=crxlable,contour=contour,GiveBIC=GiveBIC,Npix=Npix,dof=dof

  cgdisplay,700,700
  !p.multi=0
  ;xs : data array of the samples. [ndim, nsamples]
  ;     Array is either a one-D or a 2d array 
   
  ;range : array of ranges for each dimension [ndim, range]

  
  sx = size(xs)
  
  if (sx[0] lt 1) then message, 'The input sample array must be 1- or 2-D'
  if (sx[0] gt sx[1]) then message, "How cute! You want more dimensions than samples."

  if (keyword_set(stats) eq 1) then stats = make_array(1,/DOUBLE)

  
  
  
  
  ndims = (sx[0] eq 2) ? sx[1] : 1
  if (n_elements(titles) gt 0) and (n_elements(titles) ne ndims) then message,$
     "Title array needs to have same number of dimensions as the data array i.e. 9 parameters = 9 titles"
 
  if KEYWORD_Set(GIVEBIC) eq 0 then info = make_array(ndims,11,/FLOAT) else info = make_array(ndims,12,/FLOAT)
  
  ;Plot position array 
  posgrid = make_array(ndims,ndims,/STRING)
  
  
  
  xmin = 0.1
  xmax = 0.9
  ymin = 0.1
  ymax = 0.9
  len = xmax - xmin
  xarr = fltarr(ndims+1)
  yarr = fltarr(ndims+1)
  
  for k=0, ndims do begin
     ii = float(k)
     x1 = xmin + (ii * len/ndims)
     xarr[k] = x1
     y1 = ymin + (ii * len/ndims)
     yarr[k] = y1
  endfor

  
  for ii=1l, ndims do begin
     for jj=1l, ndims do begin
        posarr = string([xarr[ii-1],yarr[jj-1],xarr[ii]-0.005,yarr[jj]-0.005])
        posgrid[ii-1,jj-1] = strcompress('['+strjoin(posarr,',')+']',/REMOVE_ALL)
     endfor
  endfor
  posgridstring = ROTATE(transpose(posgrid),3)
 
  
  
  
  ;Triangle plot
  
  
  if (ndims gt 1) then begin
     s = make_array(ndims,ndims)
     
     for jj=0, ndims-2 do s[jj+1:*,jj] = 0
     for jj=0, ndims-1 do s[jj,jj] = 1         ;1-D Histogram
     for jj=0, ndims-2 do s[jj,jj+1:*] = 2     ;2-D Histogram

     charsize = 1D - (ndims*0.07)
     ticks = 6
     for ii=0,ndims-1 do begin
        for jj=0, ndims-1 do begin
           exc = execute('position = ' + posgridstring[ii,jj])
           
           if (ii eq 0) then begin
              xtitle = ''
              ytitle = 'Posterior probability'
              xtickformat='(A1)'
              xticks = ticks
           endif else if (ii eq ndims-1) then begin
              if (n_elements(titles) eq 0) then begin
              endif else begin
                 xtitle = titles[ii]
              endelse
              xtickformat='(F6.1)'
              ytitle = ''
              xticks = ticks
           endif 
           if (ii ne 0) and (ii ne (ndims-1)) then begin
              xtitle = ''
              ytitle = ''
              xtickformat='(A1)'
              xticks = ticks
           endif
           if (s[ii,jj] eq 1) then begin
               x = xs[ii,*]
               if (n_elements(range) eq 0) then begin 
                  ran = [min(x),max(x)]
               endif else ran = range[ii,*]
               
               if keyword_set(kdeon) eq 0 then begin
                  if (ii eq 0) and (jj eq 0) then begin
                     cghistoplot, x, xrange=ran, mininput=ran[0],maxinput=ran[1],position=position,$
                                  /outline,$
                                  xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat='(A1)',$
                                  charsize=charsize,nbins=30,/FILLPOLYGON,POLYCOLOR='Light Gray',thick=2 ;,xticks=xticks
                  endif else begin
                     cghistoplot, x, xrange=ran, mininput=ran[0],maxinput=ran[1],position=position,/noerase,$
                                  /outline,$
                                  xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat='(A1)',$
                                  charsize=charsize,nbins=30,/FILLPOLYGON,POLYCOLOR='Light Gray',thick=2 ;,xticks=xticks
                  endelse
               endif else begin
                  t = range(ran[0],ran[1],1000)
                  d = akde(x[*],t)
                  if (ii eq 0) and (jj eq 0) then begin
                     cgplot, t,d,position=position,xrange=ran,$
                             xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat='(A1)',$
                             charsize=charsize,color='indian red' ;,xticks=xticks
                  endif else begin
                      cgplot, t,d,position=position,/noerase,xrange=ran,$
                             xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat='(A1)',$
                             charsize=charsize,color='indian red' ;,xticks=xticks
                    
                  endelse
                  
               endelse 
               
            endif $
            else if (s[ii,jj] eq 2) then begin
               x = xs[ii,*]
               y = xs[jj,*]
               
               if (ii eq 0) then begin
                  xtitle = ''
                  if (n_elements(titles) eq 0) then begin
                     ytitle=''
                  endif else begin
                     ytitle = titles[jj]
                  endelse 
                  xtickformat='(A1)'
                  ytickformat='(F6.2)'
                  xticks=ticks
                  yticks=ticks
               endif else if (jj eq (ndims-1) ) then begin
                  if (n_elements(titles) eq 0) then begin
                     xtitle = ''
                  endif else begin
                     xtitle = titles[ii]
                  endelse
                  xtickformat='(F6.1)'
                  ytickformat='(A1)'
                  ytitle = ''
                  xticks=ticks
                  yticks=ticks
               endif
               if (ii eq 0) and (jj eq (ndims-1)) then begin
                  if (n_elements(titles) eq 0) then begin
                     ytitle=''
                     xtitle=''
                  endif else begin
                     ytitle = titles[jj]
                     xtitle = titles[ii]
                  endelse 
                  xtickformat='(F6.1)'
                  ytickformat='(F6.1)'
                  xticks=ticks
                  yticks=ticks
               endif
               
                if (ii ne 0) and (jj ne (ndims-1)) then begin
                 
                     ytitle=''
                     xtitle=''
                 
                  xtickformat='(A1)'
                  ytickformat='(A1)'
                  xticks=ticks
                  yticks=ticks
               endif
               ;axis_format = {XTICKFORMAT:xtickformat,YTICKFORMAT:ytickformat};,XTICKS:xticks,YTICKS:yticks}

               if (n_elements(range) eq 0) then begin
                 
                  ran1 = [min(x),max(x)]
                  ran2 = [min(y),max(y)]
                  
                  if n_elements(maxins) eq 0 then xmaxin = max(x) else xmaxin = maxins[ii]
                  if n_elements(minins) eq 0 then xminin = min(x) else xminin = minins[ii]
                  
                  rangex = Max(x < xmaxin, /NAN) - Min(x > xminin, /NAN)

                  if n_elements(maxins) eq 0 then ymaxin = max(y) else ymaxin = maxins[jj]
                  if n_elements(minins) eq 0 then yminin = min(y) else yminin = minins[jj]
                  
                  rangey = Max(y < ymaxin, /NAN) - Min(y > yminin, /NAN)
                  
                  ran = make_Array(2,2,/DOUBLE)
                  ran[0,*] = ran1
                  ran[1,*] = ran2
                  
               endif else begin
                  ran1 = range[ii,*]
                  ran2 = range[jj,*]

                  ran = make_Array(2,2,/DOUBLE)
                  ran[0,*] = ran1
                  ran[1,*] = ran2

                  if n_elements(maxins) eq 0 then xmaxin = ran1[1] else xmaxin = maxins[ii]
                  if n_elements(minins) eq 0 then xminin = ran1[0] else xminin = minins[ii]
                  
                  rangex = Max(x < xmaxin, /NAN) - Min(x > xminin, /NAN)

                  if n_elements(maxins) eq 0 then ymaxin = ran2[1] else ymaxin = maxins[jj]
                  if n_elements(minins) eq 0 then yminin = ran2[0] else yminin = minins[jj]
                  
                  rangey = Max(y < ymaxin, /NAN) - Min(y > yminin, /NAN)

               endelse
               
               ;if n_elements(bins) gt 0 then (ynbin = bins) else (ynbin = [])
               ;if n_elements(bins) gt 0 then (xnbin = bins) else (xnbin = [])
               
               ;IF N_Elements(xnbin) EQ 0 THEN BEGIN ; Scott's Choice
               ;   xbinsize = (3.5D * StdDev(xminin > x < xmaxIn, /NAN))/N_Elements(x)^(1./3.0D)
               ;   xnbin = (rangex / xbinsize ) + 1
               ;ENDIF ELSE BEGIN
               ;   xbinsize = rangex / (xnbin -1)
               ;ENDELSE
               
               
               
               
               ;IF N_Elements(ynbin) EQ 0 THEN BEGIN ; Scott's Choice
               ;   ybinsize = (3.5D * StdDev(xminIn > x < xmaxIn, /NAN))/N_Elements(x)^(1./3.0D)
               ;   ynbin = (rangey / ybinsize ) + 1
               ;ENDIF ELSE BEGIN
               ;   ybinsize = rangey / (ynbin -1)
               ;ENDELSE
               


               ;density = Hist_2D(x, y, Min1=xminin, Max1=xmaxin, Bin1=xbinsize, $
                                 ;Min2=yminin, Max2=ymaxin, Bin2=ybinsize)
   
               ;maxDensity = Ceil(Max(density)/1e2) * 1e2
               ;scaledDensity = BytScl(density, Min=0, Max=maxDensity)
               
               ;cgLoadCT, 33
               ;TVLCT, cgColor('white',/Triple), 0
               ;TVLCT, r, g, b, /Get
               ;palette = [ [r], [g], [b] ]
   
               
               if (ii eq 0) and (jj eq 0) then begin               
                  ;cgImage, scaledDensity, XRange=ran1, YRange=ran2, /Axes, Palette=palette, $
                  ;         XTitle=xtitle, YTitle=ytitle, $
                  ;         Position=position,$
                  ;         charsize=charsize,$
                  ;         AXKEYWORDS=axis_format
                  
                  post_plots,reform(x),refor(y),nbins=nbins,range=ran,position=position,xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat=ytickformat,/med,charsize=charsize,/contour,/smooth
               endif else begin
                  ;cgImage, scaledDensity, XRange=ran1, YRange=ran2, /Axes, Palette=palette, $
                  ;         XTitle=xtitle, YTitle=ytitle, $
                  ;         Position=position,/noerase,$
                  ;         charsize=charsize,$
                  ;         AXKEYWORDS=axis_format
                  
                  post_plots,reform(x),reform(y),nbins=nbins,range=ran,position=position,xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,ytickformat=ytickformat,/med,charsize=charsize,/contour,/smooth,/noerase
               endelse
               ;Over plot input crosshairs 
               ;if (n_elements(cross) gt 0) then begin

               ;   crx_size = size(cross)
                  
               ;   if (crx_size[0] gt 1) then begin
               ;   crox1 = cross[ii,*]
               ;   crox2 = cross[jj,*]
                  
               ;   if n_elements(crxcolor) eq 0 then crxcolor = make_array(crx_size[0],/string, VALUE = 'red')
               ;   for iii=0, crx_size[0]-1 do begin
               ;      cgplot, cross[ii,iii], cross[jj,iii], /overplot, psym=7, symsize=1.2,color=crxcolor[iii]
               ;   endfor
               ;endif else begin
               ;   crox1 = cross[ii]
               ;   crox2 = cross[jj]
               ;   cgplot, crox1, crox2, /overplot, psym=7, symsize=1.2,color='red'
               ;endelse
               ;endif

               
               ;thick = (!D.Name EQ 'PS') ? 6 : 2
               ;sdensity = SMOOTH(density,2)
               ;levels = cgConLevels( sdensity, Factor=0.01,nlevels=4)
               ;LEVELS=maxDensity*[0.68, 0.95, 0.997]
               
               ;cgContour, sdensity, LEVELS=maxDensity*[0.25, 0.5, 0.75], /OnImage, $
               ;           C_Colors='black',LABEL=0,$
               ;           C_Thick=0.5
               ;if (keyword_set(contour) eq 1) then begin
               ;   cgContour,  sdensity,levels=levels,/OnImage, $
               ;               C_Colors='white',LABEL=0,$
               ;               C_Thick=0.5
               ;endif
            endif
            
            
         endfor
 
        
     endfor
     
  endif

  


  ;Posterioir probability plots
  ;!p.multi=0
  ;if (keyword_set(posplots) eq 1) then begin
  ;Ncol = FIX(round(ndims/2D))
  
  ;if (ndims gt 1) then !P.MULTI = [0,2,Ncol]
  
  ;for ii=0, ndims-1 do begin
     
  ;   x = xs[ii,*]
  ;   
  ;   if (n_elements(range) eq 0) then begin 
  ;      ran = [min(x),max(x)]
  ;   endif else ran = range[ii,*]
  ;   if (n_elements(titles) eq 0) then begin
  ;      xtitle = ''
  ;   endif else begin
  ;      xtitle = titles[ii]
  ;   endelse
     
  ;   t = range(ran[0],ran[1],1000)
  ;   d = kde(x[*],t)
     
  ;   cgplot, t,d,xrange=ran,$
  ;           xtitle=xtitle,ytitle='Posterior probability',$
  ;           charsize=1.,color='indian red'
     
  ;endfor
;endif 
  ;!p.multi=0
  
  ;Box plots
  ;if (keyword_set(boxplot) eq 1) then begin
  ;   Ncol = FIX(round(ndims/2D))
     
  ;   if (ndims gt 1) then !P.MULTI = [0,2,Ncol]
     
  ;   for ii=0, ndims-1 do begin
        
  ;      x = xs[ii,*]
        
  ;      if (n_elements(range) eq 0) then begin 
  ;         ran = [min(x),max(x)]
  ;      endif else ran = range[ii,*]
  ;      if (n_elements(titles) eq 0) then begin
  ;         xtitle = ''
  ;      endif else begin
  ;         xtitle = titles[ii]
  ;      endelse
        
        
  ;      cgBoxplot,x,stats=stats,ytitle=xtitle
  ;   endfor
  ;endif
  
;Statistical information
  if (keyword_set(info) eq 1) then begin 
     for ii=0, ndims-1 do begin
        percen = percentiles(x,value=[0.16,0.82])
        
        info[ii,0] = median(x);stats.median
        info[ii,1] = mean(x);stats.mean
        info[ii,2] = percen[0]
        info[ii,3] = percen[1]
        info[ii,4] = min(x);stats.min
        info[ii,5] = max(x);stats.max

        percen = percentiles(x,value=[0.25,0.75])

        info[ii,6] = percen[0];stats.Q25
        info[ii,7] = percen[1];stats.Q75
        info[ii,8] = percen[1] - percen[0];stats.IQR
        info[ii,9] = STDDEV(x);stats.SDEV
        info[ii,10] = n_elements(x);stats.N

        if KEYWORD_Set(GIVEBIC) eq 1 then begin
           sdata = size(data,/DIMENSIONS)
           BIC = -1D * 2D * ALOG(MAX(data[sdata[0]])) - (dof*ALOG(NPIX))
           info[ii,11] = BIC
        endif
        
     endfor
     
     !p.multi=0
  endif
end


pro tri_example
  
  x1 = RNORM( 10000, AVG=0., SIG=0.5) 
  x2 = RNORM( 10000, AVG=2., SIG=1.) 
  x = [x1,x2]
  
  y = randomn(seed, 20000)
  
  z = RNORM( 20000, AVG=20., SIG=1.)

  k = RNORM( 20000, AVG=2., SIG=1.)
  
  l = RNORM( 20000, AVG=0., SIG=10.)

  l2 = RNORM( 20000, AVG=1., SIG=3.)
  
  l3 = RNORM( 20000, AVG=2., SIG=1.)

  l4 = RNORM( 20000, AVG=1.5, SIG=1.)

  l5 = RNORM( 20000, AVG=0.5, SIG=0.1)

  l6 = RNORM( 20000, AVG=45., SIG=10.)

  l7 = RNORM( 20000, AVG=22., SIG=2.)


  data = transpose([[x],[y],[z]]);,[k],[l],[l2],[l3],[l4],[l5],[l6],[l7]])
  data = transpose([[x],[y],[z],[k]]);,[l],[l2],[l3],[l4],[l5],[l6],[l7]])
  
  titles = ['x','y','z','k']
  cgPS_Open,strcompress('./Plots/tri_test.ps')
  triplot,data,bins=25,titles=titles,info=info

  
  ;triplot,data,bins=25,titles=titles,info=info,/GiveBIC,Npix=Npix,dof=dof


  cgPS_close
  cgPS2PDF, './Plots/tri_test.ps'
  stop
end
