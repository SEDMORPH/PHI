pro BGmkfile,file

  readcol,file,a,FORMAT='A',SKIPLINE=1,/QUICK,/SILENT
  nfile = n_elements(a)
  OPENR, lun, file, /GET_LUN
  line = ''
  READF, lun, line
  FREE_LUN, lun
  order = STRSPLIT(line, /EXTRACT)
  nstr = N_elements(order)
  strlist = ['NAME','seeing','seeing2','betaa',$
             'serX0', 'serY0', 'Ie', 'Re', 'n', 'serE', 'serPA',$
             'expX0', 'expY0', 'I0', 'h', 'expE', 'expPA']
  Forlist = ['A',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F']
  
  nlist = n_elements(strlist)

  NAME = MAKE_ARRAY(nfile, /STRING)
  seeing = fltarr(nfile)
  seeing2 = fltarr(nfile)
  betaa = fltarr(nfile)
  serX0 = fltarr(nfile)
  serY0 = fltarr(nfile)
  Ie = fltarr(nfile)
  Re = fltarr(nfile)
  n = fltarr(nfile)
  serE = fltarr(nfile)
  serPA = fltarr(nfile)
  expX0 = fltarr(nfile)
  expY0 = fltarr(nfile)
  I0 = fltarr(nfile)
  h = fltarr(nfile)
  expE = fltarr(nfile)
  expPA = fltarr(nfile)
 
  Result = MAKE_ARRAY(nstr, nlist, /INTEGER)
  for i=0l,nstr-1 do begin
     for ii=0l,nlist-1 do begin
        Result[i,ii] = STRCMP(order[i], strlist[ii])
     endfor
  endfor
  
  Order = MAKE_ARRAY(nstr, /STRING)
  Format = MAKE_ARRAY(nstr, /STRING)
  
  
  for i=0l,nstr-1 do begin
     index = where(Result[i,*] eq 1)
     Order[i] = strlist[index]
     Format[i] = Forlist[index]
  endfor
  
  Order = STRJOIN(Order,',')
  Format = STRJOIN(Format,',')
  
  
  redat = 'readcol,file,'+Order+',FORMAT='+string(39B)+Format+string(39B)+',SKIPLINE=1,/SILENT'
  ;CALL_PROCEDURE,'readcol',file,NAME,seein,serX0,serY0,Ie,Re,n,serE,serPA,expX0,expY0,I0,h,expE,expPA,FORMAT=FORMAT,SKIPLINE=1,SILENT=1 
  read = execute(redat)
  
  
  eno = MAKE_ARRAY(nlist, /INTEGER)
  index = MAKE_ARRAY(nstr, /INTEGER)

  for i=0l,nstr-1 do begin
     index[i] = where(Result[i,*] eq 1)
  endfor
  
  eno[index] = 1

  HEADER = STRJOIN(strlist,' ')
  openw,2,'./inputs/bagal_inputs.txt',width=500
  printf,2,HEADER


  
  for j=0l, nfile-1 do begin 
     for i=0, nlist-1 do begin
        case i of 
           0:begin
              if eno[i] eq 0 then begin
                 message, 'Need to supply file names'
              endif  
           end
           1:begin
              if eno[i] eq 0 then begin
                 seeing[j] = -1
              endif 
           end
           2:begin
              if eno[i] eq 0 then begin
                 seeing2[j] = -1
              endif 
           end
           3:begin
              if eno[i] eq 0 then begin
                 betaa[j] = -1
              endif 
           end
           4:begin
              if eno[i] eq 0 then begin
                 serX0[j] = -1
              endif 
           end
           5:begin
              if eno[i] eq 0 then begin
                 serY0[j] = -1
              endif 
           end
           6:begin
              if eno[i] eq 0 then begin
                 Ie[j] = -1
              endif 
           end
           7:begin
              if eno[i] eq 0 then begin
                 Re[j] = -1
              endif  
           end
           8:begin
              if eno[i] eq 0 then begin
                 n[j] = -1
              endif 
           end
           9:begin
              if eno[i] eq 0 then begin
                 serE[j] = -1
              endif
           end
           10:begin
              if eno[i] eq 0 then begin
                 serPA[j] = -1
              endif  
           end
           11:begin
              if eno[i] eq 0 then begin
                 expX0[j] = -1
              endif 
           end
           12:begin
              if eno[i] eq 0 then begin
                 expY0[j] = -1
              endif  
           end
           13:begin
              if eno[i] eq 0 then begin
                 I0[j] = -1
              endif  
           end
           14:begin
              if eno[i] eq 0 then begin
                 h[j] = -1
              endif  
           end
           15:begin
              if eno[i] eq 0 then begin
                 expE[j] = -1
              endif 
           end
           16:begin
              if eno[i] eq 0 then begin
                 expPA[j] = -1
              endif 
           end
        endcase 
     endfor
     printf,2,NAME[j],seeing[j],seeing2[j],betaa[j],$
             serX0[j], serY0[j], Ie[j], Re[j], n[j], serE[j], serPA[j],$
             expX0[j], expY0[j], I0[j], h[j], expE[j], expPA[j]
  endfor
  close,2

  
end
