;;+
; NAME:
;   bagal_AM
; PURPOSE:
;  Implementation of the Adaptive Metropolis sampler of Roberts and Rosenthal (2009). 
; INPUTS:
;   like       - Returns the log probability density to sample from 
;   like_func  - function that computes the likelihood
;   nstep      - number of iterations 
;   pars       - initial parameters (can be an array or structure)
; KEYWORDS:
;   Scale      - Vector with variences or covarience matricies of jump
;                distribution 
;   Adapt      - Adaptor, number of steps to adapt to 
;   Acc_rate   - Desired acceptance rate 
;   list       - List matrix of stuff
;   nstart     - Iteration where adaption should start
; OUTPUTS:
;   pars       - array of parameters, sorted from most to least likely
;   like       - array of likelihoods for the pars
; BUGS:
; REVISION HISTORY:
;   
;-
pro FnY,a
  n = n_elements(a) 
  if n lt 1 then begin
    ;print,' usage: randomize,SomeArray'
    ;print,'   Elements of SomeArray are placed in random order'
    return
  endif
  for i=n-1,1,-1 do begin
    ;; pick at random one of the first i elements
    indy = floor(randomu(seed)*(i+1))
    ;; move this element to the ith position
    tmp = a[i]
    a[i] = a[indy]
    a[indy] = tmp
  endfor
end

FUNCTION geweke, x, burn = burn, zstat = zstat
  
  x = REFORM(x)
  frac1 = 0.1D 
  frac2 = 0.5D
  
  ndims = n_elements(x)
  
  xmin = 0l
  
  REPEAT BEGIN
     xstart = [xmin, FLOOR(ndims - 1 - frac2*(ndims - 1 - xmin))]
     xend = [ceil(xmin + frac1*(ndims - 1 - xmin)),ndims-1]
     
     xvar = make_array(2,/DOUBLE)
     xmean = make_array(2,/DOUBLE)
     for ii=0,1 do begin
        xvar[ii] = VARIANCE(x[xstart[ii]:xend[ii]]) 
        xmean[ii] = mean(x[xstart[ii]:xend[ii]]) 
     endfor
     zstat = (xmean[0] - xmean[1]) / (sqrt(xvar[0] + xvar[1])) 
     
     if (ABS(zstat) ge 3) then frac1 = frac1 + 0.1D

  ENDREP UNTIL (ABS(zstat) le 3) OR (frac1 eq frac2)
  
  ;gewstat = 0; bad, gewstat = 1; good
  gewstat = 1
  if (frac1 eq frac2) or (ABS(zstat) ge 3) then gewstat = 0 
  
  burn = xend[0]
  return, gewstat
END
Function RCAUY, N, AVG=AVG, GAM=GAM, seed=seed
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(GAM) then GAM = 1D
  
  
  R = make_array(N,/DOUBLE)
  for i=0, N-1 do begin
     R[i] = RANDOMU(seed)
  endfor
  
  
  result = GAM * (TAN(!DPI * (R - 0.5D))) + AVG
  return, result 
end

Function RNORM, N, AVG=AVG, SIG=SIG, seed=seed
  
  if NOT keyword_set(avg) then AVG = 0D
  if NOT keyword_set(SIG) then SIG = 1D
 
  
  
  ;R = fltarr(N)
  ;rx = fltarr(N)
  ;ry = fltarr(N)
  
  ;for i=0, N-1 do begin 
  ;   repeat begin 
  ;      x = 2*RANDOMU(seed, 1)-1
  ;      y = 2*RANDOMU(seed, 1)-1
  ;      R_s = x^2 + y^2
  ;   endrep until R_s LT 1.0 
  ;   R[i] = sqrt(R_s)          
  ;   rx[i] = x
  ;   ry[i] = y
  ;endfor
  
  ;G1 = sqrt(-ALOG(R))*2D*rx/R
  ;G2 = AVG+(G1*SIG)
  R = make_array(N,/DOUBLE)
  for i=0, N-1 do begin
     R[i] = RANDOMN(seed)
  endfor
  
  Result = (R * SIG) + AVG
  return, result
END

Function nearPSD, c
  n = (size(c))[2]
  eval = LA_EIGENPROBLEM(c,eigenvectors=evec)
  val = real_part(eval)
  
  order = REVERSE(sort(val))
  val = val[order]
  index = where(val lt 0,count)
  if count ne 0 then val[index] = 0.
  
  vec = transpose(real_part(evec))
  vec = vec[order,*]
  aaa = 1D / ((vec^2D) ## val)
  T = sqrt(DIAG_MATRIX(aaa))
  bb = DIAG_MATRIX(sqrt(val))
  B = (T ## vec ## bb)
  out = B ## (transpose(B))
  return, out
end

FUNCTION MVTNORM, N=N,k=k,means=means,sig=sig
  
  if keyword_set(N) eq 0 then N = 1
  if keyword_set(k) eq 0 then k = 2
  if keyword_set(means) eq 0 then means = make_array(k,/FLOAT)
  if keyword_set(sig) eq 0 then begin
     kvec = range(1.,1.,k)
     sig = diag_matrix(kvec)
  endif

  if (n_elements(sig) eq 1) then begin
     if (sig ge -1.) and (sig le 1.) then begin
        siga = [1.,sig,sig,1.]
        sig = make_array(2,2,/FLOAT)
        sig[*] = siga
     endif
  endif
  
  tsig = sig ne transpose(sig)
  if (total(tsig) gt 0) then message, 'Incorrect covarience matrix specified'
  sigsize = size(sig,/DIMENSIONS)
  if (n_elements(means) ne k) OR (sigsize[0] ne k) OR (sigsize[1] ne k) then message, 'Incorrect covariance matrix dimensions.'
  
  
  ;eval = (HQR(ELMHES(sig), /DOUBLE))
  ;evec = transpose(real_part(EIGENVEC(sig, eval)))
  ;eval = real_part(eval)
  eval = EIGENQL(sig, EIGENVECTORS = evec)
  eval = round(eval*(1.e10))/(1.e10)
  
  if (where(eval lt 0)) gt -1 then message, 'Inappropriate covariance matrix specified.' ;begin 
  
     ;sig = nearPSD(sig)          
     ;;eval = EIGENQL(sig, EIGENVECTORS = evec)
     ;eval = round(eval*(1.e10))/(1.e10)
  ;endif
  
  
  A = evec ##  (diag_matrix(sqrt(eval)))
  x = make_array(k,n,/float)
  ran = RNORM(n * n_elements(means))
  x[*] = ran
  output = means + A ## x

  return, output
END

FUNCTION rmvtnorm, n, mu, sig
  d = n_elements(mu)
  temp = TRANSPOSE(Nrandom(SIG=(make_array(d,/DOUBLE,value=1D))))
  S = sig
  ;LA_CHOLDC, S 
  CHOLDC, S, Pin, /DOUBLE
  for jj = 0,d-2 Do s[jj+1:*,jj] = 0
  for jj = 0,d-1 DO s[jj,jj] = Pin[jj]
  rvector = mu + (S##temp)
  return, rvector
END


PRO PHI_AM,  like, nstep, pars,like_func,scale=scale,$
               adapt=adapt, accrate=accrate,$
               nstart=nstart,nrep=nrep,info=info
  common setup3, flist, nlist, mlist, astep,i, af0,groups
  common files, psfname, inname, maskname, whtname, imaname, dir_out, psf
  common chistuff, isig, numbit, Nof
  common specs, rmin, sigover, fmask, fweight,err_fac,rmax,tklog,maskonly
  ibefore = i
  rmaxbefore = rmax
  rminbefore = rmin

  if keyword_set(adapt) eq 1 and N_ELEMENTS(accrate) eq 0 then message, 'Accrate is missing!'
  
  ;Number of adaption steps 
  if keyword_set(adapt) eq 1 then nadapt = long(adapt)
  if keyword_set(adapt) eq 0 then nadapt = 0
  
  ;Number of parameters 
  index = where(flist[0,*] eq 1)
  npars= n_elements(Pars)
  nfpars = n_elements(Pars[index])
  multii = (2.38D)^2D/(float(nfpars))

  ;Number of groups
  ingroups = [1:max(groups)]
  ngr = n_elements(ingroups)
  fgroups = groups[index]
  

  
  
  ;Matrix to store MC chain 
  oldpars= double(pars)
  pars = MAKE_ARRAY(npars,nstep,/DOUBLE)
  pars[*,0] = oldpars
  parsig = MAKE_ARRAY(npars,nstep,/DOUBLE)
  parsig[index,0] = scale[index]
  sqrdjump = make_array(nfpars,nstep,/DOUBLE)
  avgsqrdjump = make_array(nfpars,((nstep - (0.2*nstep)))/200l)
  pvalue = make_array(nfpars,/DOUBLE)
  
  ;Array to store log densities p(x) or log likelihoods 
  like = MAKE_ARRAY(nstep,/DOUBLE)
  if (err_fac eq 1) then begin
     fX = 0D
     newoldpars = [oldpars,fx]
  endif else newOldPars = oldpars
  like[0] = call_function(like_func,newoldpars)
  alpha = MAKE_ARRAY(nstep,/DOUBLE)
  ;gX = MAKE_ARRAY(nstep,/DOUBLE)
  

  ;Initial S 
  if (n_elements(scale[index]) eq nfpars) then begin
     M = DIAG_MATRIX(scale[index])
  endif else begin
     M = scale[index]
  endelse 
  
  
  ;Check 
  msize = size(M,/DIMENSIONS)
  if (msize[0] ne nfpars) then message, 'Length or dimension of pars and scale do not match!'

  ;*********************************************************************************************
  ;1st Stage: ADAPTIVE METROPOLIS within GIBBS ALGORITHM
  print, '1st Adaption phase activated' 
  
  ;Coeffiecents for AMwG
  bwid_a1 = 100l
  bwid = bwid_a1
  tempAccRate = make_array(nfpars,nstep/bwid,/DOUBLE)
  nacc = make_array(nfpars,nstep,/DOUBLE)
  if (ngr eq nfpars) then begin
     minacc = 0.28D
     maxacc = 0.6D
     accaim = 0.44D
  endif else begin
     minacc = 0.15D
     maxacc = 0.4D
     accaim = 0.23D
  endelse

  endbatch = 2D
  unimult = scale[index] ;make_array(nfpars,/DOUBLE)
  PropSigma = scale[index];make_array(nfpars,/DOUBLE,value=1)
  adpScale = 0.05
  count = 0D ;Counting concecutive pi(x) = 0
  nimprob = 100D ;Max number of concecutive pi(x) = 0
  break_flag = 0
  X = pars[index,0]

  Y = X
  gX = like[0]
  Xveclist = make_array(npars,nstep,/DOUBLE)
  Xveclist[*,0] = oldpars
  PropPars = double(oldpars) 
  fsig = 0.3D
  
  
  for ii=1l, nstep-1l do begin
     
     if (count ge 50) then begin
        print, 'Level 1 failure: zero acceptance'
        stop
        goto, mistake1
     endif
     
     
     
     ;Adjust sigma for each coordinate based on local acceptance rate
     k = ii
     if (k gt 1l) AND ((k MOD bwid) eq 0) then begin
        newindex = ceil(ii/bwid)
        for jj=0, nfpars-1 do begin
           tempAccRate[jj,newindex-1] = total(nacc[jj,(ii-bwid):(ii)])/bwid
        endfor

        ;If acceptance rate between minacc and maxacc then stop adaption
        if ((min(tempAccRate[*,newindex-1])) gt minacc) AND ((max(tempAccrate[*,newindex-1]) lt maxacc)) then begin

           k = bwid + 1l
           bwid = 2D*bwid
           print, 'Hot!', tempAccRate[*,newindex-1]
           if (bwid eq bwid_a1*2D^(endbatch+1)) then begin
              adpstop = ii
              transstop = ii
              break
           endif
        endif else begin

           
           if (max(tempAccRate[*,newindex-1]) lt 0.70) AND (min(tempAccRate[*,newindex-1]) gt 0.15) then print, 'Warm', tempAccRate[*,newindex-1] else begin
              print, 'Cold!', tempAccRate[*,newindex-1]
              count = count+1D              
           endelse


           for j=0, nfpars-1 do begin
              
              if (tempAccRate[j,newindex-1] lt 0.70) AND (tempAccRate[j,newindex-1] gt 0.15) then begin
                 
                 if (tempAccRate[j,newindex-1] gt accaim) then begin
                    unimult[j] = ALOG(propSigma[j]) + adpScale
                 endif else begin
                    unimult[j] = ALOG(propSigma[j]) - adpScale
                 endelse

              endif else begin

                 if (tempAccRate[j,newindex-1] gt accaim) then begin
                    unimult[j] = ALOG(propSigma[j]) + 1D
                 endif else begin
                    unimult[j] = ALOG(propSigma[j]) - 1D
                 endelse

              endelse

              propSigma[j] = exp(unimult[j])
           endfor
           
        endelse

     endif 
     


;Here we can randomise the order using a Fisher and yates algorithm
     inorder = [0:nfpars-1]
     outorder = ingroups
     Fny,outorder

     ;outgroups = ingroups
     ;Fny,outgroups
     for jj=0, ngr-1 do begin
        ;nj = outorder[jj]
        nj = where(fgroups eq outorder[jj])
        
        flagcnt = 0D
        
        REPEAT BEGIN
           Y[nj] = X[nj] + RNORM((n_elements(nj)),SIG=propSigma[nj]) ;PROPOSE A VALUE
           ;Y[nj] = X[nj] + RCAUY((n_elements(nj)),GAM=propSigma[nj])
           
           PropPars[(index[nj])] = Y[nj]
         
           pflag = bagal_flag(PropPars)
           flagcnt = flagcnt + 1
           if (flagcnt eq 200) then message, '200 iters gone by without good parameters'
           
        ENDREP UNTIL (pflag eq 1) 
        ;if ((10D^y[1])/(10D^y[6])) gt 1.678 then stop
        ;if ((10D^proppars[3])/(10D^proppars[10])) gt 1.678 then stop

        ;T = SYSTIME(1)
        gY = call_function(like_func,PropPars)
        ;PRINT, 'MCMC took ', (SYSTIME(1) - T), ' Seconds to complete'
        ;stop

        ;Acceptance probability 
        Pracc = exp(gY - gX)
        
        ;stop
        if (gX eq 0D) then begin
           A = 1D
           count = count+1D
        endif else begin
           A = min([1D,Pracc])
           count = 0D
        endelse
        
        U = RANDOMU(seed)
         ;Accept or regect
        if (U lt A) then begin
           gX = gY
           X[nj] = Y[nj]
           ;if (err_fac eq 1) then fX = fY
           nacc[nj,ii] = nacc[nj,ii]+1l
        endif else begin
           Y[nj] = X[nj]
           ;if (err_fac eq 1) then fY = fX
        endelse
        
     endfor

     ;Save new values
     ;if ((10D^x[1])/(10D^x[6])) gt 1.678 then stop
     PropPars[index] = X
     Xveclist[*,ii] = PropPars ;Save X value
     like[ii] = gX
     ;if ((10D^x[1])/(10D^x[6])) gt 1.678 then stop

     if (ii MOD 50l eq 0) then begin
        ;Iscreen,Xveclist,like,nlist,ii
     endif
     ;if ((10D^x[1])/(10D^x[6])) gt 1.678 then stop
  endfor
  
  
  if (count ge (nimprob*(DOUBLE(nfpars)))) then message, (strcompress(('Past '+string(nimprob)+' X values all have probabilities of 0 under the stationary distribution')))
  if (ii eq nstep) then begin
     print, 'MCMC chain ran for the maximum amount of iterations without convergence'
     goto, mistake1
  endif

  print,''
  print, strcompress('Final acceptance rates: '+ strjoin(string(tempAccRate[*,newindex-1]),', '))
  print, strcompress('Final proposal sigmas: '+ strjoin(string(Propsigma),', '))
  Print, strcompress('1st Adaption phase deactivated at i = '+ string(adpstop))

  ;*********************************************************************************************
  ;Transient phase 
  
  ;Conditions and values
  pvalue = make_array(nfpars,/DOUBLE)
  bwidth = 200l
  trans = make_array(nfpars,(nstep/bwidth),/DOUBLE)
  nreg = 5l
  transstop = adpstop + bwidth*5
  ;Main code
  ;for ii=adpstop, transstop-1 do begin
  for ii=adpstop, nstep-1 do begin
     inorder = [0:nfpars-1]
     ;outorder = inorder
     ;Fny,outorder
     
     outorder = ingroups
     Fny,outorder
     for jj=0, ngr-1 do begin
        ;nj = outorder[jj]
        nj = where(fgroups eq outorder[jj])
        
        REPEAT BEGIN
           Y[nj] = X[nj] + RNORM((n_elements(nj)),SIG=propSigma[nj]) ;PROPOSE A VALUE
           ;Y[nj] = X[nj] + RCAUY((n_elements(nj)),GAM=propSigma[nj])
           PropPars[(index[nj])] = Y[nj]
           pflag = bagal_flag(PropPars)
        ENDREP UNTIL (pflag eq 1)
        gY = call_function(like_func,PropPars)
        ;Acceptance probability 
        Pracc = exp(gY - gX)
        U = (randomu(seed))
        if (gX eq 0D) then begin
           A = 1D
           count = count+1D
        endif else begin
           A = min([1D,Pracc])
           count = 0D
        endelse
        
         ;Accept or regect
        if (U lt A) then begin
           gX = gY
           X[nj] = Y[nj]
           nacc[nj,ii] = nacc[nj,ii]+1l
        endif else begin
           Y[nj] = X[nj]
        endelse
        
     endfor
     ;Save new values
     PropPars[index] = X
     Xveclist[*,ii] = PropPars  ;Save X value
     like[ii] = gX
     if (ii MOD bwidth eq 0) then begin
        newindex = ceil(ii/bwidth)
        
        for jj=0, nfpars-1 do begin
           trans[jj,newindex] = mean(Xveclist[index[jj],(ii-bwidth):ii])
           if (ceil((ii-adpstop)/bwidth) gt (nreg - 1)) then begin
                 responce = trans[jj,(newindex - (nreg - 1)):newindex]
                 exp = indgen(nreg+1)
                 remove,0,exp
                 fit  = LINCORR(responce, exp)
                 pvalue[jj] = fit[1] 
              endif
        endfor
        
        if (min(pvalue) gt 0.1) then begin
           transstop = ii
           break
        endif else if (ii eq transstop) then break
        
     endif
     if (ii MOD 50l eq 0) then begin
        ;Iscreen,Xveclist,like,nlist,ii
     endif
    

  endfor
  
  if (count ge (nimprob*(DOUBLE(nfpars)))) then message, (strcompress(('Past '+string(nimprob)+' X values all have probabilities of 0 under the stationary distribution')))
  if (ii eq nstep) then begin
     print, 'MCMC chain ran for the maximum amount of iterations without convergence'
     goto, mistake1
  endif

  print, ''
  print, 'Transient phase finished at i = ', transstop
  
 

  !p.multi = 0
  ;*********************************************************************************************
  ;2nd ADAPTIVE Stage: ADAPTIVE METROPOLIS ALGORITHM
  print,''
  print, '2nd Adaption phase activated' 
  
  if (n_elements(nlist) le 7) then nplots = 5 else nplots = 9
  !P.Multi = [0, 2,nplots]

  k = 0l
  hh = 0l
  fudg = 1l
  kk = 0l
  nreg = 7l
  stop_flag = 0l
  bwidth = 200l
  pval = make_array(nfpars,/DOUBLE)
  naccept = make_array(nstep-transstop,/DOUBLE)
  minaccept = 0.02D
  maxaccept = 0.6D
  
  for ii=transstop,nstep-1L do begin
     
     if (stop_flag eq 0) then begin
        if total(flist[0,*]) eq 1 then begin
           covar = VARIANCE(Xveclist[index,(transstop - (bwidth)):ii-1])
        endif else begin
           covar = CORRELATE(Xveclist[index,(transstop - (bwidth)):ii-1], /COVARIANCE) 
           ;covar = CORRELATE(Xveclist[index,0:ii-1], /COVARIANCE) 
        endelse
        propSigma = covar * multii 
     endif 
     
     
     pflag = 0
     REPEAT BEGIN 
        if (nfpars eq 1) then begin
           parsi = (RNORM(1,sig=sqrt(propSigma))) 
        endif else begin
           sigsize = size(propSigma,/DIMENSIONS)
           ;if (ii MOD 20 eq 0) then stop
           parsi = rmvtnorm( 1, (make_array(nfpars,/DOUBLE)), propsigma)
           
           ;parsi = (MVTNORM(k=sigsize[0],sig=propSigma))
        endelse 
        parx = FINITE(parsi, /NAN)
        if total(parx) gt 0 then stop  
        Prop =  Xveclist[*,ii-1]
        Propvals = (Xveclist[index,ii-1] + (parsi)); + B*(pars[index,ii-1] + (addFac))
        Prop[index] = Propvals
        pflag = bagal_flag(Prop)
        
     ENDREP UNTIL (pflag eq 1) 
     
     parsig[index,ii] = parsi
     
     ;print, propSigma
 ;Calculate the density at prop
     gY = call_function(like_func,Prop)
  
   
     ;Acceptance probability 
     Pracc = exp(gY - gX)
     alpha[ii] = min([1D,Pracc]) ;For log density/likelihood
     ;if alpha[ii] lt 0.002 then alpha[ii] = 0.002D
     randomnum = (randomu(seed))
     
     ;Accept with P=alpha 
     if (randomnum lt alpha[ii]) then begin  
        Xveclist[*,ii] = Prop    ;Accept 
        gX = gY
        k = k+1l
        naccept[ii-transstop] = naccept[ii-transstop]+1l
     endif else begin
        Xveclist[*,ii] = Xveclist[*,ii-1]    ;Or not 
     endelse
     like[ii] = gX
     sqrdjump[*,ii-transstop] = (Xveclist[index,ii-1] - Xveclist[index,ii])^2D
     
    
     if ((ii-transstop) MOD bwidth eq 0) and (ii - transstop) ne 0 then begin
        newindex = ceil((ii-transstop)/bwidth)
           
        for j=0, nfpars-1 do begin
           avgsqrdjump[j,newindex] = mean(sqrdjump[j,(ii-transstop-bwidth) : (ii-transstop)])
           if (ceil((ii-transstop)/bwidth) gt (nreg - 1)) then begin
              responce = avgsqrdjump[j,(newindex - (nreg - 1)):newindex]
              exp = indgen(nreg+1)
              remove,0,exp
              fit  = LINCORR(responce, exp)
              pval[j] = fit[1]
           endif
        endfor
        if (min(pval) gt 0.1) then begin
           adapt_stop = ii
           print,''
           print, 'Adaption has stopped at i = ', adapt_stop
           break
        endif
     endif 
  
     
     if (ii MOD 100l eq 0) then begin
        ;Iscreen,Xveclist,like,nlist,ii
       
        if (ii-transstop) eq 200l then begin
           accept_rate = TOTAL(naccept[0:ii-transstop])/(DOUBLE(ii-transstop))
        
           if (accept_rate lt minaccept) then begin
              ;multii = multii/(max([2D,nfpars])) 
              multii = multii/(float(nfpars))
              Xveclist[*,ii] = Xveclist[*,transstop]
              like[ii] = like[transstop]
              ii = transstop
              print, ''
              print, strcompress('Current acceptance rate is: '+string(accept_rate))
              print, 'Adaption has failed! Chian restarted at i = ', ii
             
              count = count+1D
           endif
        endif
        
     endif
     
     if (count ge (nimprob*(DOUBLE(nfpars)))) then message, (strcompress(('Past '+string(nimprob)+' X values all have probabilities of 0 under the stationary distribution')))
     ;print, $
     ;   strjoin(strarr(21)+string(byte(8)),''), $
     ;   'MCMC: ',100L*ii/nstep,' % ', $
     ;   format= '($,A,A,I2.2,A)'
  endfor
  
  print,''
  accept_rate = TOTAL(naccept[0:ii-transstop])/(DOUBLE(ii-transstop))
  print, strcompress('Final acceptance rate: '+string(accept_rate))
  print, 'End of 2nd Adaption phase' 
  
  
  !p.multi = 0
  ;*********************************************************************************************
  ;Sampling phase: METROPOLIS ALGORITHM
  !p.multi = [0, 1,nplots]
  count = make_array(nrep,/DOUBLE)
  starttemp = make_array(nfpars,nrep)
  startdist = 0;0.5D
  ;startdist = (startdist-1D)/2D
  a = make_array(2,nfpars,/DOUBLE)
  support = make_array(nfpars,/DOUBLE)
  
  X = make_array(nrep,npars,/DOUBLE)
  gX = make_array(nrep,/DOUBLE)
  
  
  for k=0, nrep-1 do begin
     pflag = 0
     X[k,*] = Xveclist[*,adapt_stop]
     REPEAT BEGIN
        for l=0, nfpars-1 do begin
           
           gew = geweke(Xveclist[index[l],(transstop):adapt_stop], burn = burn, zstat = zstat)
           ;a[0,l] = min(Xveclist[index[l],burn:adapt_stop])
           ;a[1,l] = max(Xveclist[index[l],burn:adapt_stop])
           a[0,l] = (Xveclist[index[l],adapt_stop])
           a[1,l] = (Xveclist[index[l],adapt_stop])
           dist = a[1,l] - a[0,l]
           starttemp[l,*] = Uprior(nrep,(a[0,l]-(startdist*dist)),(a[1,l]+(startdist*dist)))
           ;starttemp[l,*] = [Xveclist[index[l],adapt_stop],Xveclist[index[l],adapt_stop],Xveclist[index[l],adapt_stop]]
        endfor 
        
       
        X[k,index] = starttemp[*,k]
        pflag = bagal_flag(X[kk,*])
     ENDREP UNTIL (pflag eq 1)
     gX[k] = call_function(like_func,X[k,*])
  endfor
  print, ''
  print, 'Starting values for the chian(s): '
  for l=0l, nrep-1 do begin
     print, strcompress('Chain '+string(l+1l)+': ')+strcompress(strjoin(string(transpose(X[l,index])),'  '))
  endfor
  print, ''
  Y = make_array(nrep,npars,/DOUBLE)
  gY = make_array(nrep,/DOUBLE)
  like = make_array(nrep,(nstep-adapt_stop),/DOUBLE)

  Xveclistdim = make_array(nrep,npars,(nstep-adapt_stop),/DOUBLE)
  naccept = make_array(nrep,nstep-adapt_stop,/DOUBLE)
  
  ;Split remaining iterations to have number of iterations per chain 
  ichain = ceil((nstep - adapt_stop)/nrep)
  numper = nstep-adapt_stop
  runmean  = make_array(nrep,npars,(nstep-adapt_stop),/DOUBLE)

  ;if numper lt 0 then numper = nstep-adapt_stop
  avelist = make_array(nrep,nfpars,numper,/DOUBLE)
  s2 = avelist & xbar = avelist
  xbartotal = make_array(nfpars,(numper/bwidth))
  varc = avelist & varmean = xbartotal & vhat = xbartotal & Rvar = xbartotal & varchain = xbartotal & Rhat = xbartotal

  Rmin = 0.9D
  Rmax = 1.1D
  
 
  gew = make_array(nrep,nfpars,/DOUBLE)
  burnin = make_array(nrep,nfpars,/DOUBLE)
 
  acc_rate = make_array(nrep,(nstep-adapt_stop),/DOUBLE)
  chainstart = make_array(nrep,/DOUBLE)
  minacc = 0.1D
  maxacc = 0.6D

  ;chainstop = nstep
  
  ;MAIN CODE
  for ii=LONG(adapt_stop), nstep-1 do begin
     for kk=0, nrep-1 do begin
        
        pflag = 0
        REPEAT BEGIN 
           if (nfpars eq 1) then begin
              parsi = (RNORM(1,sig=sqrt(propSigma))) 
           endif else begin
              ;parsi = (MVTNORM(k=nfpars,sig=propSigma))
              parsi = rmvtnorm(1, (make_array(nfpars,/DOUBLE)), propsigma)
           endelse 
           parx = FINITE(parsi, /NAN)
           if total(parx) gt 0 then stop  
           Y[kk,*] = X[kk,*]
           Propvals = (X[kk,index] + (parsi)) ; + B*(pars[index,ii-1] + (addFac))
           Y[kk,index] = Propvals
           pflag = bagal_flag(Y[kk,*])
           
        ENDREP UNTIL (pflag eq 1)
     
        
     ;print, propSigma
 ;Calculate the density at prop
        gY[kk] = call_function(like_func,Y[kk,*])
  

        ;Acceptance probability 
        Pracc = exp(gY[kk] - gX[kk])
        alpha = min([1D,Pracc]) 
        if alpha lt 0.002 then alpha = 0.002D
        randomnum = (randomu(seed))
     
        ;Accept with P=alpha 
        if (randomnum lt alpha) then begin  
           Xveclistdim[kk,*,ii-adapt_stop] = transpose(Y[kk,*]) ;Accept 
           X[kk,*] = Y[kk,*]
           gX[kk] = gY[kk]
           naccept[kk,ii-adapt_stop] = naccept[kk,ii-adapt_stop]+1l
        endif else begin
           Xveclistdim[kk,*,ii-adapt_stop] = transpose(X[kk,*]) ;Or not 
        endelse
        like[kk,ii-adapt_stop] = gX[kk]
     endfor
     
    
 
     
     if (ii MOD 500l eq 0) and (ii-adapt_stop gt 0) then begin
         ;Iscreen,Xveclistdim[*,*,0:ii-adapt_stop],like[*,0:ii-adapt_stop],nlist,ii-adapt_stop,cmulti=nrep
     endif 

     ;AVERAGE ACCEPTANCE RATE
     if (ii-adapt_stop eq (ceil(ichain/2D))) then begin
           for kk=0, nrep-1 do begin
              newindex = ceil((ii - adapt_stop)/(ceil(ichain/2D)))
              acc_rate[kk,newindex] = TOTAL(naccept[kk,(ceil(ichain/4D)):ii-adapt_stop])/(DOUBLE(ichain/4D))
              if (acc_rate[kk,newindex] lt minacc) or (acc_rate[kk,newindex] gt maxacc) then begin
                 print, ''
                 print, strcompress('Failure in chain: '+string(kk+1))
                 print, strcompress('Acceptance rate: '+string(acc_rate[kk,newindex]))
                 chainstart[kk] = ii-adapt_stop
                 X[kk,*] = Xveclist[*,adapt_stop]
                 X[kk,index] = starttemp[*,kk]
                 gX[kk] = call_function(like_func,X[kk,*])
              endif
        endfor
        
        
     endif

     if (nrep gt 1) then begin
        if (ii-adapt_stop ge ichain) AND (ii MOD bwidth) then begin
           
           newindex = ceil((ii - adapt_stop - ichain)/bwidth)
           
           for jj=0, nfpars-1 do begin
              for kk=0, nrep-1 do begin
                 
                 
                 ;GEWEKE DIAGNOSTIC FOR BURN-In iteration
                 gew[kk,jj] = geweke(Xveclistdim[kk,index[jj],0:ii-adapt_stop], burn = burn, zstat = zstat)
                 burnin[kk,jj] = burn
                 

              endfor
           endfor
           
           
           discard = (max(burnin[*,*]))
           samplesize = ii - adapt_stop - discard
          
           
           if (min(gew) gt 0) OR (ii ge 0.95D*nstep) then begin
              
              for jj=0, nfpars-1 do begin
                 for kk=0, nrep-1 do begin
                    
                 ;Gelman-Rubin DIAGNOSTIc
                    avelist[kk,jj,newindex] = TOTAL(Xveclistdim[kk,index[jj],(ii-adapt_stop-ichain+discard):(ii-adapt_stop)])/samplesize ;xbar_c -> mean of each chain
                    varc[kk,jj,newindex] = TOTAL((Xveclistdim[kk,index[jj],(ii-adapt_stop-ichain+discard):(ii-adapt_stop)] - avelist[kk,jj,newindex])^2D)/(samplesize -1D) ;var of each chain
                 endfor
                 
                 xbartotal[jj,newindex] = mean(Xveclistdim[*,index[jj],(ii-adapt_stop-ichain+discard):(ii-adapt_stop)]) ;xbar -> mean of all chains
                 varchain[jj,newindex] = TOTAL(varc[*,jj,newindex])/nrep                                        ;mean of individual chain vars 
                 varmean[jj,newindex] = TOTAL((avelist[*,jj,newindex] - xbartotal[jj,newindex])^2D)/nrep        ;var of chain means
                 
                 Rhat[jj,newindex] = ((((samplesize - 1D)/samplesize)*varchain[jj,newindex]) + (varmean[jj,newindex]/samplesize)) / varchain[jj,newindex]
              endfor
              
              if (max(Rhat[*,newindex]) lt Rmax) and (min(Rhat[*,newindex]) gt Rmin) then begin
                 chainstop = ii
                 print, ''
                 print, 'Sampling phase finished at i = ',chainstop
                 break
              endif
           endif 
        endif  
     endif  
     ;print, $
     ;   strjoin(strarr(21)+string(byte(8)),''), $
     ;   'MCMC: ',100L*ii/nstep,' % ', $
     ;   format= '($,A,A,I2.2,A)'
  endfor
  
  !p.multi = 0
  print, strjoin(strarr(21)+string(byte(8)),'')+'MCMC: done      '
  print, ''
  ;*********************************************************************************************
   
  ;Save Stuff for plotting later 
  SAVE,Xveclistdim,filename=dir_out+'TraceHistory.sav'
  SAVE,Xveclist,filename=dir_out+'Pre_Trace.sav'
  ;SAVE,ACF,filename=dir_out+'ACF.sav'
  ;SAVE,acfsig,filename=dir_out+'ACFsigma.sav'
  ;SAVE,ESS,filename=dir_out+'ESS.sav'
  sumchain = [adpstop,transstop,adapt_stop,discard,chainstop]
  SAVE,sumchain,filename=dir_out+'ChainInfo.sav'
  

  sampsize = (chainstop - adapt_stop)
  repits = (sampsize - discard)
  print, strcompress('R-hat for free parameters: '+strjoin(string(Rhat[*,newindex]),', '))
  print, 'Acceptance rate(s): '
  
  acc_rate = make_array(nrep,/DOUBLE)
  for i=0l, nrep-1 do begin
     acc_rate[i] = TOTAL(naccept[i,(discard):sampsize])/(DOUBLE(repits))
     print, strcompress('Chain '+string(i+1l)+': '+string(acc_rate[i]))
  endfor
  print, ''


 ;CALCULATE THE ACF and ESS
  nlags = 50l
  ACF = make_array(nrep,nfpars,nlags+1,/DOUBLE)
  acfsig = make_array(nrep,nfpars,/DOUBLE)
  ESS = make_array(nrep,nfpars,/DOUBLE)
  for kk=0,nrep-1 do begin
     for jj=0, nfpars-1 do begin
        ACF[kk,jj,*] =  autocorr(Xveclistdim[kk,index[jj],burnin[kk,jj]:sampsize],nlags=nlags)
        acfsig[kk,jj] = 2D*sqrt((1D/(sampsize-burnin[kk,jj]))*(1D + 2D*(TOTAL(ACF[kk,jj,*]^2D))))
        essindex = where(ACF[kk,jj,*] ge (2D*acfsig[kk,jj]))
        ESS[kk,jj] = DOUBLE(sampsize-burnin[kk,jj]) /(1D + (2D *TOTAL((acf[kk,jj,*])[essindex])))
     endfor
  endfor



  ;Rejecting chains with poor acceptance rates   
  chnacc = make_array(nrep,/double)
  accindex = where(acc_rate ge minacc)
  if (min(accindex) ne -1) then chnacc[accindex] = 1
  chburn = make_array(nrep,/double)
  
  if (min(accindex) ne -1) then begin
     newpars = Xveclistdim[accindex,*,*]
     newlike = like[accindex,*]
     newnrep = n_elements(accindex)
     Print, strcompress('Chain(s): '+strjoin(string(accindex+1),' & ')+' have been accepted for the final sample')
     pars = make_array(npars,1,/DOUBLE)
     like = make_array(1,/DOUBLE)
     burn = make_array(newnrep,/DOUBLE)
     for i=0l, newnrep-1 do begin
        true_burn = make_array(nfpars,/DOUBLE)
        for j=0, nfpars-1 do begin
           gew = geweke(newpars[i,index[j],chainstart[i]:sampsize-1],burn=burninn)
           true_burn[j] = burninn
        endfor
        burnmax = max(true_burn)
        burn[i] = burnmax
        pars = [[pars],[reform(newpars[i,*,(chainstart[i]+burnmax):sampsize-1])]]
        like = [like,(reform(newlike[i,(chainstart[i]+burnmax):sampsize-1]))]
        print, strcompress('Burn-in for chain '+string(i+1)+': '+string(chainstart[i]+burnmax))
     endfor
     chburn[accindex] = burn
  endif else begin
     Print, 'All the chains have failed'
     newnrep = nrep
     newpars = Xveclistdim[*,*,*]
     newlike = like[*,*]
     pars = make_array(npars,1,/DOUBLE)
     like = make_array(1,/DOUBLE)
     for ii=0l, newnrep-1 do begin
        pars = [[pars],[reform(newpars[ii,*,(chainstart[ii]):sampsize-1])]]
        like = [like,(reform(newlike[ii,(chainstart[ii]):sampsize-1]))]
     endfor
  endelse
  pars = pars[*,1:*]
  like = like[1:*]
  
 
  print, strcompress('Sample total: '+string(n_elements(like)))

  nsamp = (n_elements(like))

  mistake1:

  Info = [nrep,nsamp,adpstop,transstop,adapt_stop,chainstop]
  for ii=0,nrep-1 do begin
     Info = [info,chnacc[ii],chburn[ii]]
  endfor

  i = ibefore 
  rmax = rmaxbefore
  rmin = rminbefore
END



function psd, time, data, freq, plot=plot, log=log

  dt = time[1]-time[0]
  n = n_elements(time)
  
  freq = findgen(round(n/2))/dt/n
  fx = fft(data)
  psd = abs(fx[0:round(n/2)-1])^2
  
  if keyword_set(log) then psd = alog10(psd)
  if keyword_set(log) then ytitlep = ' (logarithmic)' else ytitlep = ''
  
  
  if keyword_set(plot) then plot, freq, psd, xtitle='Frequency (in Hz)', ytitle = 'Power Spectral Density'+ytitlep
  
  return, psd
end


Pro Iscreen, pars,like,nlist,ii,Cmulti=Cmulti
  common chistuff, isig, numbit, Nof, Pri
  common setup2, scale, zcalib, tpsf, betaa, fFFT, seeing, sky, sigsky,fixcen
  
  sorder = where(nlist eq 'serX0', sernn)
  eorder = where(nlist eq 'expX0', expnn)
  cenorder = where(nlist eq 'X0', cennn)
  skyorder = where(nlist eq 'sky', skynn)

  if (n_elements(nlist) le 7) then nplots = 5 else nplots = 9
  !P.Multi = [0,1,nplots]

  charsize = 1.7
  if (keyword_set(Cmulti) eq 0) then begin
     if (expnn eq 0) then begin
           miny = min((10D^pars[2,0:ii]))-(0.01*min((10D^pars[2,0:ii])))
           maxy = max((10D^pars[2,0:ii]))+(0.01*max((10D^pars[2,0:ii])))
           cgplot, (10D^pars[2,0:ii]),yr=[miny,maxy], yt='Ie',xt='Iteration',position=[0.1,0.6,0.3,0.8],charsize=charsize
           
           miny = min(10D^pars[3,0:ii])-(0.01*min(10D^pars[3,0:ii]))
           maxy = max(10D^pars[3,0:ii])+(0.01*max(10D^pars[3,0:ii]))
           cgplot, 10D^(pars[3,0:ii]), yt='Re',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.35,0.6,0.55,0.8],charsize=charsize
           
           miny = min(10D^pars[4,0:ii])-(0.01*min(10D^pars[4,0:ii]))
           maxy = max(10D^pars[4,0:ii])+(0.01*max(10D^pars[4,0:ii]))
           cgplot, (10D^pars[4,0:ii]), yt='n',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.6,0.6,0.8,0.8],charsize=charsize

           miny = min(pars[5,0:ii])-(0.01*min(pars[5,0:ii]))
           maxy = max(pars[5,0:ii])+(0.01*max(pars[5,0:ii]))
           cgplot, pars[5,0:ii], yt='b/a',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.1,0.35,0.3,0.55],charsize=charsize
           
           pa = pars[6,0:ii];acos(pars[6,0:ii])*180D/!DPI
           ;pa = acos(pars[6,0:ii])*180D/!DPI
           miny = min(pa)-(0.01*min(pa))
           maxy = max(pa)+(0.01*max(pa))
           cgplot, pa, yt='PA',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.35,0.35,0.55,0.55],charsize=charsize
           ;cgplot, pars[6,0:ii],/overplot
        endif else begin
           miny = min(10^(pars[2,0:ii]))-(0.01*min(10^(pars[2,0:ii])))
           maxy = max(10^(pars[2,0:ii]))+(0.01*max(10^(pars[2,0:ii])))
           cgplot, 10^(pars[2,0:ii]),yr=[miny,maxy], yt='Ie',xt='Iteration',position=[0.04,0.6,0.24,0.8],charsize=charsize
           
           miny = min(10D^pars[3,0:ii])-(0.01*min(10D^pars[3,0:ii]))
           maxy = max(10D^pars[3,0:ii])+(0.01*max(10D^pars[3,0:ii]))
           cgplot, 10D^(pars[3,0:ii]), yt='Re',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.29,0.6,0.49,0.8],charsize=charsize
           
           miny = min(10D^pars[4,0:ii])-(0.01*min(10D^pars[4,0:ii]))
           maxy = max(10D^pars[4,0:ii])+(0.01*max(10D^pars[4,0:ii]))
           cgplot, (10D^pars[4,0:ii]), yt='n',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.6,0.74,0.8],charsize=charsize

           miny = min(pars[5,0:ii])-(0.01*min(pars[5,0:ii]))
           maxy = max(pars[5,0:ii])+(0.01*max(pars[5,0:ii]))
           cgplot, pars[5,0:ii], yt='b/a',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.79,0.6,0.99,0.8],charsize=charsize
           
           pa = pars[6,0:ii];acos(pars[6,0:ii])*180D/!DPI
           ;pa = acos(pars[6,0:ii])*180D/!DPI
           miny = min(pa)-(0.01*min(pa))
           maxy = max(pa)+(0.01*max(pa))
           cgplot, pa, yt='PA',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.04,0.35,0.24,0.55],charsize=charsize
           
           miny = min(10D^pars[9,0:ii])-(0.01*min(10D^pars[9,0:ii]))
           maxy = max(10D^pars[9,0:ii])+(0.01*max(10D^pars[9,0:ii]))
           cgplot, 10D^pars[9,0:ii], yt='I0',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.29,0.35,0.49,0.55],charsize=charsize

           miny = min(10D^pars[10,0:ii])-(0.01*min(10D^pars[10,0:ii]))
           maxy = max(10D^pars[10,0:ii])+(0.01*max(10D^pars[10,0:ii]))
           cgplot, (10D^pars[10,0:ii]), yt='h',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.35,0.74,0.55],charsize=charsize
           
           miny = min(pars[11,0:ii])-(0.01*min(pars[11,0:ii]))
           maxy = max(pars[11,0:ii])+(0.01*max(pars[11,0:ii]))
           cgplot, pars[11,0:ii],yt='b/a',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.79,0.35,0.99,0.55],charsize=charsize
           
           ;pa = acos(pars[12,0:ii])*180D/!DPI
           pa = pars[12,0:ii];acos(pars[12,0:ii])*180D/!DPI
           miny = min(pa)-(0.01*min(pa))
           maxy = max(pa)+(0.01*max(pa))
           cgplot, pa,yt='PA',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.04,0.1,0.24,0.3],charsize=charsize
           ;cgplot, pars[12,0:ii],/overplot
        endelse

        ;!P.Multi = [0,2,1]
  
        ;Centre coordinates 
        if fixcen eq 0 then begin
           
           minx = min(pars[0,ii])-(0.01*min(pars[0,ii]))
           maxx = max(pars[0,ii])+(0.01*max(pars[0,ii]))
           miny = min(pars[1,ii])-(0.01*min(pars[1,ii]))
           maxy = max(pars[1,ii])+(0.01*max(pars[1,ii]))
           
           cgplot,pars[0,0:ii],yt='X',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.29,0.1,0.49,0.3],charsize=charsize
           cgplot, pars[1,0:ii],yt='Y',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.1,0.74,0.3],charsize=charsize

           
        endif else begin
           cenorder = where(nlist eq 'X0', cennn)
          
         
           minx = min(pars[cenorder,0:ii])-(0.001*min(pars[cenorder,0:ii]))
           maxx = max(pars[cenorder,0:ii])+(0.001*max(pars[cenorder,0:ii]))
           miny = min(pars[cenorder+1,0:ii])-(0.001*min(pars[cenorder+1,0:ii]))
           maxy = max(pars[cenorder+1,0:ii])+(0.001*max(pars[cenorder+1,0:ii]))
           
           cgplot,pars[cenorder,0:ii],yt='X',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.29,0.1,0.49,0.3],charsize=charsize
           cgplot, pars[cenorder+1,0:ii],yt='Y',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.1,0.74,0.3],charsize=charsize

        endelse

        if skynn gt 0 then begin
           
           minx = min(pars[skyorder,0:ii])-(2D*min(pars[skyorder,0:ii]))
           maxx = max(pars[skyorder,0:ii])+(2D*max(pars[skyorder,0:ii]))
           
           cgplot,pars[skyorder,0:ii],yt='sky',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.79,0.1,0.99,0.3],charsize=charsize
        endif

        ;Badness of fit and red-chisq
        
        badness = -2D*((like[0:ii] - pri[0])/10e4)
        redchi = -2D*(like[0:ii] - pri[0] + isig + numbit) / Nof        
        badness = -2D*((like[0:ii])/10e4)
        redchi = -2D*(like[0:ii] + isig + numbit) / Nof        

        cgplot, badness,yt='-2 ln(L) [x10$\up4$]',xt='Iteration',charsize=2.,/noerase,position=[0.05,0.85,0.27,0.95],symsize=0.2
        cgplot, redchi, yt='$\chi$$\downred$$\up2$ ',xt='Iteration',charsize=2.,/noerase,position=[0.54,0.85,0.79,0.95],symsize=0.2,xr=[0,n_elements(redchi)],yr=[min(redchi) - (min(redchi/10)), max(redchi) + (max(redchi/10))]

        cgplot, badness[(ii-49):ii],yt='',xt='Iteration',charsize=1.,/noerase,position=[0.95,0.95,0.99,0.99],/nodata,ytickformat='(A1)'
        cgplot, badness[(ii-49):ii],/overplot
  endif else begin

     nrep = cmulti
     xx = range(0,n_elements(pars[0,0,*]),n_elements(pars[0,0,*]))
     colors = Round(cgScaleVector(Findgen(nrep), 0, 255))
     cgLoadCT, 34
     color = plotcolor(indgen(nrep))
     Xveclistdim = pars
     
     if (expnn eq 0) then begin
        
        miny = min(10D^(Xveclistdim[*,2,*]))-(0.01*min(Xveclistdim[*,2,0:*]))
        maxy = max(10D^(Xveclistdim[*,2,*]))+(0.01*max(Xveclistdim[*,2,0:*]))
        cgplot, xx, 10D^(Xveclistdim[*,2,*]),yr=[miny,maxy], yt='Ie',xt='Iteration',/nodata,position=[0.1,0.6,0.3,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, 10D^(Xveclistdim[kk,2,*]),/overplot,color=color[kk]
        endfor
        
        miny = min(10D^Xveclistdim[*,3,0:*])-(0.01*min(10D^Xveclistdim[*,3,*]))
        maxy = max(10D^Xveclistdim[*,3,0:*])+(0.01*max(10D^Xveclistdim[*,3,*]))
        cgplot, xx, 10D^(Xveclistdim[0,3,*]), yt='Re',xt='Iteration',yr=[miny,maxy],/nodata,/noerase,position=[0.35,0.6,0.55,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, 10D^(Xveclistdim[kk,3,0:*]),/overplot,color=color[kk]
        endfor
        
        
        miny = min(10D^Xveclistdim[*,4,0:*])-(0.01*min(10D^Xveclistdim[*,4,*]))
        maxy = max(10D^Xveclistdim[*,4,0:*])+(0.01*max(10D^Xveclistdim[*,4,*]))
        cgplot, xx, (10D^Xveclistdim[0,4,*]),yt='n',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.6,0.6,0.8,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, (10D^Xveclistdim[kk,4,*]),/overplot,color=color[kk]
        endfor
        
        miny = min(Xveclistdim[*,5,*])-(0.01*min(Xveclistdim[*,5,*]))
        maxy = max(Xveclistdim[*,5,*])+(0.01*max(Xveclistdim[*,5,*]))
        cgplot, xx, Xveclistdim[0,5,*],yr=[miny,maxy], yt='b/a',xt='Iteration',/nodata,/noerase,position=[0.1,0.35,0.3,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, Xveclistdim[kk,5,*],/overplot,color=color[kk]
        endfor
        pa = Xveclistdim[*,6,*];acos(Xveclistdim[*,6,*])*180D/!DPI
        ;pa = acos(Xveclistdim[*,6,*])*180D/!DPI
        miny = min(pa)-(0.01*min(pa))
        maxy = max(pa)+(0.01*max(pa))
        cgplot, xx, pa,yt='ser PA',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.35,0.35,0.55,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           pa = Xveclistdim[kk,6,*];acos(Xveclistdim[kk,6,*])*180D/!DPI
           ;pa = acos(Xveclistdim[kk,6,*])*180D/!DPI
           cgplot, xx, pa,/overplot,color=color[kk]
        endfor      

        !P.Multi = [0,2,1]
        
        ;Centre coordinates 
        
        if fixcen eq 0 then begin
           xc0 = pars[0,0,ii]-5
           xc1 = pars[0,0,ii]+5
           yc0 = pars[0,1,ii]-5
           yc1 = pars[0,1,ii]+5

           cgplot, pars[0,0,ii],pars[0,1,ii],psym=2.,/noerase,symsize=1.5,position=[0.1,0.85,0.2,0.95],/nodata,xr=[xc0,xc1],yr=[yc0,yc1],charsize=0.7,color='dark grey'
           for kk=0, nrep-1 do begin
              cgplot, pars[kk,0,ii],pars[kk,1,ii],psym=2.,/noerase,symsize=1.5,color=color[kk],/overplot
              if (n_elements(nlist) gt 7) then cgplot, pars[kk,7,ii],pars[kk,8,ii],psym=1,/overplot,symsize=2.,color=color[kk]
           endfor
        endif else begin
           cenorder = where(nlist eq 'X0', cennn)
          
           
           minx = min(pars[*,cenorder,0:ii])-(0.001*min(pars[*,cenorder,0:ii]))
           maxx = max(pars[*,cenorder,0:ii])+(0.001*max(pars[*,cenorder,0:ii]))
           miny = min(pars[*,cenorder+1,0:ii])-(0.001*min(pars[*,cenorder+1,0:ii]))
           maxy = max(pars[*,cenorder+1,0:ii])+(0.001*max(pars[*,cenorder+1,0:ii]))
           
           cgplot,xx,pars[0,cenorder,0:ii],yt='X',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.29,0.1,0.49,0.3],charsize=charsize,/nodata
           for kk=0, nrep-1 do begin
              cgplot, xx,pars[kk,cenorder,0:ii],/overplot,color=color[kk]
           endfor

           cgplot,xx,pars[0,cenorder+1,0:ii],yt='Y',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.1,0.74,0.3],charsize=charsize,/nodata
           for kk=0, nrep-1 do begin
              cgplot, xx,pars[kk,cenorder+1,0:ii],/overplot,color=color[kk]
           endfor
        endelse
        ;Badness of fit
        badness = -2D*(like[*,*]/10e4)
        
        cgplot, badness[0,*],yt='-2 ln(L) [x10$\up4$]',xt='Iteration',charsize=1.,/noerase,position=[0.3,0.85,0.55,0.95],psym=-16,symsize=0.4,/nodata
        for kk=0, nrep-1 do begin
           cgplot, badness[kk,*],/overplot,color=color[kk]
        endfor
        
        redchi = -2D*(like[*,*] + isig + numbit) / NoF
       
        cgplot, redchi, yt='$\chi$$\downred$$\up2$ ',xt='Iteration',charsize=1.,/noerase,position=[0.85,0.6,0.99,0.75],psym=16,symsize=0.2,/nodata,yr=[0,max(redchi[*,ceil(n_elements(redchi[0,*])/10D):n_elements(redchi[0,*])-1])],xr=[0,n_elements(redchi[0,*])]
        for kk=0, nrep-1 do begin
           cgplot, redchi[kk,*],/overplot,color=color[kk]
        endfor

        newbadness = badness[*,(ii-49):ii];/(max(ABS(badness[*,(ii-49):ii])))
        cgplot, newbadness[0,*],yt='',xt='Iteration',charsize=1.,/noerase,position=[0.6,0.85,0.8,0.95],ytickformat='(A1)',/nodata
        for kk=0, nrep-1 do begin
           cgplot, badness[kk,*],/overplot,color=color[kk]
        endfor
     
     endif else begin
        
        miny = min(10^Xveclistdim[*,2,*])-(0.01*min(10^Xveclistdim[*,2,*]))
        maxy = max(10^Xveclistdim[*,2,*])+(0.01*max(10^Xveclistdim[*,2,*]))
        cgplot, xx, 10^Xveclistdim[0,2,*],yr=[miny,maxy], yt='Ie',xt='Iteration',/nodata,position=[0.04,0.6,0.24,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, 10^Xveclistdim[kk,2,*],/overplot,color=color[kk]
        endfor
        
        miny = min(10D^Xveclistdim[*,3,*])-(0.01*min(10D^Xveclistdim[*,3,*]))
        maxy = max(10D^Xveclistdim[*,3,*])+(0.01*max(10D^Xveclistdim[*,3,*]))
        cgplot, xx, 10D^(Xveclistdim[0,3,*]), yt='Re',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.29,0.6,0.49,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, 10D^(Xveclistdim[kk,3,*]),/overplot,color=color[kk]
        endfor
        
        miny = min(10D^Xveclistdim[*,4,*])-(0.01*min(10D^Xveclistdim[*,4,*]))
        maxy = max(10D^Xveclistdim[*,4,*])+(0.01*max(10D^Xveclistdim[*,4,*]))
        cgplot, xx, (10D^Xveclistdim[0,4,*]),yt='n',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.54,0.6,0.74,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, (10D^Xveclistdim[kk,4,0:*]),/overplot,color=color[kk]
        endfor
        
        miny = min(Xveclistdim[*,5,*])-(0.01*min(Xveclistdim[*,5,*]))
        maxy = max(Xveclistdim[*,5,*])+(0.01*max(Xveclistdim[*,5,*]))
        cgplot, xx, Xveclistdim[0,5,*],yr=[miny,maxy], yt='b/a',xt='Iteration',/nodata,/noerase,position=[0.79,0.6,0.99,0.8],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, Xveclistdim[kk,5,*],/overplot,color=color[kk]
        endfor
           
        pa = Xveclistdim[*,6,*];acos(Xveclistdim[*,6,*])*180D/!DPI
        ;pa = acos(Xveclistdim[*,6,*])*180D/!DPI
        miny = min(pa)-(0.01*min(pa))
        maxy = max(pa)+(0.01*max(pa))
        cgplot, xx, pa,yt='ser PA',xt='Iteration',yr=[miny,maxy],/nodata,position=[0.04,0.35,0.24,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           pa = Xveclistdim[kk,6,*];acos(Xveclistdim[kk,6,*])*180D/!DPI
           ;pa = acos(Xveclistdim[kk,6,*])*180D/!DPI
           cgplot, xx, pa,/overplot,color=color[kk]
        endfor
        
        miny = min(10D^Xveclistdim[*,9,*])-(0.01*min(10D^Xveclistdim[*,9,*]))
        maxy = max(10D^Xveclistdim[*,9,*])+(0.01*max(10D^Xveclistdim[*,9,*]))
        cgplot, xx, 10D^Xveclistdim[0,9,*],yr=[miny,maxy], yt='I0',xt='Iteration',/nodata,/noerase,position=[0.29,0.35,0.49,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, 10D^Xveclistdim[kk,9,*],/overplot,color=color[kk]
        endfor
        
        
        miny = min(10D^Xveclistdim[*,10,*])-(0.01*min(10D^Xveclistdim[*,10,*]))
        maxy = max(10D^Xveclistdim[*,10,*])+(0.01*max(10D^Xveclistdim[*,10,*]))
        cgplot, xx, (10D^Xveclistdim[0,10,*]),yt='h',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.54,0.35,0.74,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, (10D^Xveclistdim[kk,10,*]),/overplot,color=color[kk]
        endfor         
        
        miny = min(Xveclistdim[*,11,*])-(0.01*min(Xveclistdim[*,11,*]))
        maxy = max(Xveclistdim[*,11,*])+(0.01*max(Xveclistdim[*,11,*]))
        cgplot, xx, Xveclistdim[0,11,*],yr=[miny,maxy],yt='exp b/a',xt='Iteration',/nodata,/noerase,position=[0.79,0.35,0.99,0.55],charsize=charsize
        for kk=0, nrep-1 do begin
           cgplot, xx, Xveclistdim[kk,11,*],/overplot,color=color[kk]
        endfor
        pa = Xveclistdim[*,12,*];acos(Xveclistdim[*,12,*])*180D/!DPI
        ;pa = acos(Xveclistdim[*,12,*])*180D/!DPI
        miny = min(pa)-(0.01*min(pa))
        maxy = max(pa)+(0.01*max(pa))
        cgplot, xx, pa,yt='exp PA',xt='Iteration',/nodata,yr=[miny,maxy],/noerase,position=[0.04,0.1,0.24,0.3],charsize=charsize
        for kk=0, nrep-1 do begin
           pa = Xveclistdim[kk,12,*];acos(Xveclistdim[kk,12,*])*180D/!DPI
           ;pa = acos(Xveclistdim[kk,12,*])*180D/!DPI
           cgplot, xx, pa,/overplot,color=color[kk]
        endfor
        
        
        
        ;Centre coordinates 
        if fixcen eq 0 then begin
          
           
        endif else begin
           cenorder = where(nlist eq 'X0', cennn)
           
           
           minx = min(pars[*,cenorder,0:ii])-(0.001*min(pars[*,cenorder,0:ii]))
           maxx = max(pars[*,cenorder,0:ii])+(0.001*max(pars[*,cenorder,0:ii]))
           miny = min(pars[*,cenorder+1,0:ii])-(0.001*min(pars[*,cenorder+1,0:ii]))
           maxy = max(pars[*,cenorder+1,0:ii])+(0.001*max(pars[*,cenorder+1,0:ii]))
           
           cgplot,xx,pars[0,cenorder,0:ii],yt='X',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.29,0.1,0.49,0.3],charsize=charsize,/nodata
           for kk=0, nrep-1 do begin
              cgplot, xx,pars[kk,cenorder,0:ii],/overplot,color=color[kk]
           endfor

           cgplot,xx,pars[0,cenorder+1,0:ii],yt='Y',xt='Iteration',yr=[miny,maxy],/noerase,position=[0.54,0.1,0.74,0.3],charsize=charsize,/nodata
           for kk=0, nrep-1 do begin
              cgplot, xx,pars[kk,cenorder+1,0:ii],/overplot,color=color[kk]
           endfor

       
        endelse 
        
        ;Centre coordinates 
       
        if skynn gt 0 then begin
           
           minx = min(pars[*,skyorder,0:ii])-(0.1*min(pars[*,skyorder,0:ii]))
           maxx = max(pars[*,skyorder,0:ii])+(0.1*max(pars[*,skyorder,0:ii]))
           
           cgplot,pars[0,skyorder,0:ii],yt='sky',xt='Iteration',yr=[minx,maxx],/noerase,position=[0.79,0.1,0.99,0.3],charsize=charsize,/nodata
           for kk=0, nrep-1 do begin
              cgplot, xx,pars[kk,skyorder,0:ii],/overplot,color=color[kk]
           endfor

        endif

        ;Badness of fit and red-chisq
        
       
      
        ;Badness of fit
        badness = -2D*(like[*,*]/10e4)
        
        cgplot, badness[0,*],yt='-2 ln(L) [x10$\up4$]',xt='Iteration',charsize=2.,/noerase,position=[0.05,0.85,0.27,0.95],symsize=0.4,/nodata
        for kk=0, nrep-1 do begin
           cgplot, badness[kk,*],/overplot,color=color[kk]
        endfor
        
        redchi = -2D*(like[*,*] - pri[0] + isig + numbit) / NoF
       
        cgplot, redchi, yt='$\chi$$\downred$$\up2$ ',xt='Iteration',charsize=2.,/noerase,position=[0.54,0.85,0.79,0.95],symsize=0.2,/nodata,xr=[0,n_elements(redchi[0,*])],yr=[0,max(redchi)]
        for kk=0, nrep-1 do begin
           cgplot, redchi[kk,*],/overplot,color=color[kk]
        endfor

        newbadness = badness[*,(ii-49):ii];/(max(ABS(badness[*,(ii-49):ii])))
        cgplot, newbadness[0,*],yt='',xt='Iteration',charsize=0.1,/noerase,position=[0.95,0.95,0.99,0.99],ytickformat='(A1)',/nodata
        for kk=0, nrep-1 do begin
           cgplot, badness[kk,*],/overplot,color=color[kk]
        endfor
        
     endelse
     
  endelse
  
  
  
END
