function contourlevels_lc, image, enclosedfrac
	
	imlevelhist = histogram(image, binsize=0.000001, min=0.0000001, omax=maxhist, loc=imlevelloc)
	cumulativepixels = total(reverse(imlevelhist*imlevelloc), /cumulative)
	reverseloc = reverse(imlevelloc)
	fraccumulativepixels = cumulativepixels
	result = reverseloc[value_locate(fraccumulativepixels, enclosedfrac)]
	result = result[sort(result)]

return, result

end

pro post_plots,xx,yy,nbins=nbins,noerase=noerase,position=position,xtitle=xtitle,ytitle=ytitle,charsize=charsize,range=range,ytickformat=ytickformat,xtickformat=xtickformat,contour=contour,yflip=yflip,smooth=smooth,trues=trues,med=med,kdeimage=kdeimage,conover=conover,sigfi=sigfi,title=title

  
  ;cgLoadCT, 39
  cgLoadCT, 20, /REVERSE
  ;cgLoadCT, 2, NColors=3, /REVERSE
  ;cgLoadCT, 33, NColors=12, Bottom=3
     
  TVLCT, cgColor('white',/Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  if keyword_set(xtickformat) eq 0 then xtickformat = '(A1)' 
  if keyword_set(ytickformat) eq 0 then ytickformat = '(A1)' 
  axis_format = {XTICKFORMAT:xtickformat,YTICKFORMAT:ytickformat}

  if N_elements(sigfi) eq 0 then sigfi=2
  
  
  data = make_array(2,N_ELEMENTS(xx),/DOUBLE)
  data[0,*] = xx
  data[1,*] = yy

  
  density = twoD_pdf(data,bins=nbins,info=info,range=range)
  
  
 
  
  xrange = [info.min1,info.max1]
  if keyword_set(yflip) eq 0 then yrange = [info.min2,info.max2] else begin
     yrange = [info.max2,info.min2]
     density = REVERSE(density, 2)
  endelse
  densize = size(density,/DIMENSION)
  xbsize = (xrange[1] - xrange[0]) / (densize[0]-1)
  ybsize = (yrange[1] - yrange[0]) / (densize[1]-1)
  
  ;xtickV = [xrange[0]:xrange[1]]
  
  if keyword_set(kdeimage) eq 0 then begin 
     maxDensity = Ceil(Max(density)/1e2) * 1e2
     scaledDensity = BytScl(density, Min=0, Max=maxDensity)
  endif else begin
     index = where(xx ge xrange[0] and xx le xrange[1] and yy ge yrange[0] and yy le yrange[1])
     density_kde = kde2d(xx[index],yy[index],n=100, lims=[xrange,yrange])
     scaledDensity = BytScl(density_kde)
  endelse 

  if (keyword_set(contour) eq 0) then begin
     if keyword_set(noerase) eq 0 then begin               
        cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
                 XTitle=xtitle, YTitle=ytitle, $
                 Position=position,$
                 charsize=charsize,$
                 AXKEYWORDS=axis_format
     endif else begin
        cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
                 XTitle=xtitle, YTitle=ytitle, $
                 Position=position,/noerase,$
                 charsize=charsize,$
                 AXKEYWORDS=axis_format
     endelse

      if n_elements(trues) gt 0 then begin
        xtruex = range(trues[0],trues[0],100)
        xtruey = range(-100,100,100)
        cgplot, xtruex,xtruey,/overplot,thick=1 
        
        ytruey = range(trues[1],trues[1],100)
        ytruex = range(-100,100,100)
        cgplot, ytruex,ytruey,/overplot,thick=1  
     endif
     
     if keyword_set(med) eq 1 then begin
        xmedx = range(median(xx),median(xx),100)
        xmedy = range(-100,100,100)
        cgplot, xmedx,xmedy,/overplot,linestyle=1,color='Blue',thick=1 
        ymedy = range(median(yy),median(yy),100)
        ymedx = range(-100,100,100)
        cgplot, ymedx,ymedy,/overplot,linestyle=1,color='Blue',thick=1  
        
        q = fractile(xx, [.16,.84])
        xmedx = range(q[0],q[0],100)
        xmedy = range(-100,100,100)
        cgplot, xmedx,xmedy,/overplot,linestyle=2,color='Blue',thick=1 
        xmedx = range(q[1],q[1],100)
        xmedy = range(-100,100,100)
        cgplot, xmedx,xmedy,/overplot,linestyle=2,color='Blue',thick=1 
        
        q = fractile(yy, [.16,.84])
        ymedy = range(q[0],q[0],100)
        ymedx = range(-100,100,100)
        cgplot, ymedx,ymedy,/overplot,linestyle=2,color='Blue',thick=1 
        ymedy = range(q[1],q[1],100)
        ymedx = range(-100,100,100)
        cgplot, ymedx,ymedy,/overplot,linestyle=2,color='Blue',thick=1 
        
        
        
     endif
  endif else begin 
    
     if keyword_set(kdeimage) eq 1 then begin
        cgLoadCT, 51
        SetDecomposedState, 0, CurrentState=currentState
        TVLCT, cgColor('white',/Triple), 0
        TVLCT, r, g, b, /Get
        palette = [ [r], [g], [b] ]

        if keyword_set(noerase) eq 0 then begin               
           cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
                 XTitle=xtitle, YTitle=ytitle, $
                    Position=position,$
                    charsize=charsize,$
                    AXKEYWORDS=axis_format
        endif else begin
           
           cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
                    XTitle=xtitle, YTitle=ytitle, $
                    Position=position,/noerase,$
                    charsize=charsize,$
                    AXKEYWORDS=axis_format
        endelse
        SetDecomposedState, currentState
     endif 
     
   
     
     if keyword_set(smooth) eq 1 then begin 
        index = where(xx ge xrange[0] and xx le xrange[1] and yy ge yrange[0] and yy le yrange[1])
        density = kde2d(xx[index],yy[index],n=nbins, lims=[xrange,yrange])
     endif
     ;xvec = range(min(xx),max(xx),nbins)
     ;yvec = range(min(yy),max(yy),nbins)
     densize = size(density,/DIMENSION)
     xvec = range(xrange[0],xrange[1],densize[0])
     yvec = range(yrange[0],yrange[1],densize[1])
     
     ;xvec = [xrange[0],xrange[1],xbsize]
     ;yvec = [yrange[0],yrange[1],ybsize]

     ;; Normalize the total are of the 2D probability density function (pdf)
     pdf = density/TOTAL(density)
     
;;	Contour levels containing 68% and 99.7% of points, respectively
     pdfcont = contourlevels_lc(pdf, [.682,.954,.996])
    
     
   
     nlevels = n_elements(pdfcont)
     
     if Keyword_set(conover) eq 1 then begin
       
           
           
           cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,$
                      Background=cgColor('white'), $
                      Color=cgColor('black'),$
                      C_ANNOTATION=' ',C_CHARSIZE=0.01, /overplot
       
           
     endif else begin
     


     
     
     cgLoadCT, 51, NColors=5,Bottom=2 ;,/reverse
     SetDecomposedState, 0, CurrentState=currentState
     
     
     
     if keyword_set(noerase) eq 0 then begin
        
        if keyword_set(smooth) eq 0  then begin
           cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,$
                      Position=position, $
                      Background=cgColor('white'), $
                      Color=cgColor('black'), XStyle=1, YStyle=1, $
                      C_Colors=cgColor('black'),xtitle=xtitle,ytitle=ytitle,$
                      xtickformat=xtickformat,ytickformat=ytickformat,charsize=charsize,C_ANNOTATION=' ',C_CHARSIZE=0.01,$
                      xticks=3,title=title,thick=1;,xtickv=[xrange[0],median(xx),xrange[1]],yticks=2,ytickv=[yrange[0],median(yy),yrange[1]] 
        endif else begin
           cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,$
                      Position=position, $
                      /Fill, C_Colors=Indgen(nlevels)+3,$
                      Background=cgColor('white'), $
                      Color=cgColor('black'), XStyle=1, YStyle=1, $
                      xtitle='!C'+xtitle,ytitle=ytitle,$
                      xtickformat=xtickformat,ytickformat=ytickformat,charsize=charsize,C_ANNOTATION=' ',C_CHARSIZE=0.01,$
                      xticks=2,title=title,thick=1;,xtickv=[xrange[0],median(xx),xrange[1]],yticks=2,ytickv=[yrange[0],median(yy),yrange[1]] 
           Contour, pdf, xvec, yvec, /Overplot, Levels=pdfcont,NLevels=n_elements(pdfcont), /Follow, Color=cgColor('black'),C_ANNOTATION=' ',C_CHARSIZE=0.01,thick=1
        endelse
        ;cgContour, pdf, xvec,yvec,levels=pdfcont[4],/overplot
        ;cgContour, pdf, xvec,yvec,levels=pdfcont,/overplot


     endif else begin
        if keyword_set(smooth) eq 0 then begin

           cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,/noerase,$
                      Position=position, $
                      Background=cgColor('white'), $
                      Color=cgColor('black'), XStyle=1, YStyle=1,C_ANNOTATION=' ',C_CHARSIZE=0.01,$
                      C_Colors=cgColor('black'),xtitle=xtitle,ytitle=ytitle,xtickformat=xtickformat,$
                      ytickformat=ytickformat,charsize=charsize,title=title,XTHICK=1,YTHICK=1,thick=1 ;,$
                      ;xticks=3;,xtickv=[xrange[0],median(xx),xrange[1]],yticks=2,ytickv=[yrange[0],median(yy),yrange[1]]
        endif else begin
           
           if xtickformat eq '(A1)' then begin
              cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,/noerase,$
                         Position=position, $
                         /Fill, C_Colors=Indgen(nlevels)+3,$
                         Background=cgColor('white'), $
                         Color=cgColor('black'), XStyle=1, YStyle=1, $
                         xtitle='!C'+xtitle,ytitle=ytitle,$
                         xtickformat=xtickformat,ytickformat=ytickformat,charsize=charsize,C_ANNOTATION=' ',C_CHARSIZE=0.01,$
                         XTick_Get=tickValues,title=title,XTHICK=1,YTHICK=1,thick=1 ;xticks=3;,xtickv=[xrange[0],median(xx),xrange[1]],yticks=2,ytickv=[yrange[0],median(yy),yrange[1]] 
              
              Contour, pdf, xvec, yvec, /Overplot, Levels=pdfcont,NLevels=n_elements(pdfcont), /Follow, Color=cgColor('black'),C_ANNOTATION=' ',C_CHARSIZE=0.01,thick=1
              
              
           endif else begin
              
              cgContour, pdf, xvec,yvec,NLevels=n_elements(pdfcont), levels=pdfcont,/noerase,$
                         Position=position, $
                         /Fill, C_Colors=Indgen(nlevels)+3,$
                         Background=cgColor('white'), $
                         Color=cgColor('black'), XStyle=1, YStyle=1, $
                         xtitle='!C'+xtitle,ytitle=ytitle,$
                         xtickformat=xtickformat,ytickformat=ytickformat,charsize=charsize,C_ANNOTATION=' ',C_CHARSIZE=0.01,$
                         XTick_Get=tickValues,title=title,XTHICK=1,YTHICK=1,thick=1 ;xticks=3;,xtickv=[xrange[0],median(xx),xrange[1]],yticks=2,ytickv=[yrange[0],median(yy),yrange[1]] 
              
              Contour, pdf, xvec, yvec, /Overplot, Levels=pdfcont,NLevels=n_elements(pdfcont), /Follow, Color=cgColor('black'),C_ANNOTATION=' ',C_CHARSIZE=0.01,thick=1
              
              ;numticks = n_elements(tickValues)
              ;ypos = Replicate(!Y.Window[0] - 0.03, numticks+1)
              ;if numticks lt 5 then xsh = 0.04
              ;if numticks ge 5 then xsh = 0.07
 
              ;xpos = !X.Window[0] + (!X.Window[1] + xsh - !X.Window[0]) *  Findgen(numticks + 1) / numticks
              ;FOR j=0, numticks-2 DO BEGIN 
              ;   if tickvalues[j] ge -0.0001 and tickvalues[j] le 0.0001 then str = '0.00' else str = (sigfig(tickvalues[j], sigfi))
              ;XYOutS, xpos[j], ypos[j],str, Alignment=0.5, Orientation=45, /Normal, color=cgcolor('black'),CHARSIZE=charsize
              ;ENDFOR
           endelse
           
        endelse
                                ;cgContour, pdf, xvec,yvec,levels=pdfcont[4],/overplot
     endelse
     SetDecomposedState, state
     
  endelse
    if n_elements(trues) gt 0 then begin
        xtruex = range(trues[0],trues[0],100)
        xtruey = range(-100,100,100)
        cgplot, xtruex,xtruey,/overplot,color='dark grey',thick=1,linestyle=3
        
        ytruey = range(trues[1],trues[1],100)
        ytruex = range(-100,100,100)
        cgplot, ytruex,ytruey,/overplot,color='dark grey',thick=1,linestyle=3 
     endif
     
     if keyword_set(med) eq 1 then begin
        xmedx = range(median(xx),median(xx),100)
        xmedy = range(-100,100,100)
        cgplot, xmedx,xmedy,/overplot,color='dark grey',thick=1 
        ymedy = range(median(yy),median(yy),100)
        ymedx = range(-100,100,100)
        cgplot, ymedx,ymedy,/overplot,color='dark grey',thick=1  
        
     endif
                                ;Contour, density,xvec,yvec,/Overplot, $
                                ;         Color=cgColor('black'), Levels=levels, $
                                ;         C_Labels=everyOther
  
     
     SetDecomposedState, currentState
     
                                ;cgContour, density, C_Colors=Indgen(nlevels)+3, Background=cgColor('white'),NLevels=userLevels
     
  endelse




  
end
