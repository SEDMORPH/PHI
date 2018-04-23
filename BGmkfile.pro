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
             'expX0', 'expY0', 'I0', 'h', 'expE', 'expPA',$
             'ser2X0', 'ser2Y0', 'Ie2', 'Re2', 'n2', 'ser2E', 'ser2PA',$
             'ser3X0', 'ser3Y0', 'Ie3', 'Re3', 'n3', 'ser3E', 'ser3PA',$
             'ser4X0', 'ser4Y0', 'Ie4', 'Re4', 'n4', 'ser4E', 'ser4PA',$
             'bexpX0','bexpY0','bI0','bh1','bh2','bexpa','rb','bexpE','bexpPA',$
             'ferX0','ferY0','ferI0','l','fern','ferE','ferPA','ferc',$
             'fer2X0','fer2Y0','fer2I0','l2','fer2n','fer2E','fer2PA','fer2c'$
            ]
  
  Forlist = ['A',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F',$
             'F','F','F','F','F','F','F','F'$
            ]
  
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
  
  ser2X0 = fltarr(nfile)
  ser2Y0 = fltarr(nfile)
  Ie2 = fltarr(nfile)
  Re2 = fltarr(nfile)
  n2 = fltarr(nfile)
  ser2E = fltarr(nfile)
  ser2PA = fltarr(nfile)

  ser3X0 = fltarr(nfile)
  ser3Y0 = fltarr(nfile)
  Ie3 = fltarr(nfile)
  Re3 = fltarr(nfile)
  n3 = fltarr(nfile)
  ser3E = fltarr(nfile)
  ser3PA = fltarr(nfile)

  ser4X0 = fltarr(nfile)
  ser4Y0 = fltarr(nfile)
  Ie4 = fltarr(nfile)
  Re4 = fltarr(nfile)
  n4 = fltarr(nfile)
  ser4E = fltarr(nfile)
  ser4PA = fltarr(nfile)

  bexpX0 = fltarr(nfile)
  bexpY0 = fltarr(nfile)
  bI0 = fltarr(nfile)
  bh1 = fltarr(nfile)
  bh2 = fltarr(nfile)
  bexpa = fltarr(nfile)
  rb = fltarr(nfile)
  bexpE = fltarr(nfile)
  bexpPA = fltarr(nfile)

  ferX0 = fltarr(nfile)
  ferY0 = fltarr(nfile)
  ferI0 = fltarr(nfile)
  l = fltarr(nfile)
  fern = fltarr(nfile)
  ferE = fltarr(nfile)
  ferPA = fltarr(nfile)
  ferc = fltarr(nfile)

  fer2X0 = fltarr(nfile)
  fer2Y0 = fltarr(nfile)
  fer2I0 = fltarr(nfile)
  l2 = fltarr(nfile)
  fer2n = fltarr(nfile)
  fer2E = fltarr(nfile)
  fer2PA = fltarr(nfile)
  fer2c = fltarr(nfile)
  
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
           17:begin
              if eno[i] eq 0 then begin
                 ser2X0[j] = -1
              endif 
           end
           18:begin
              if eno[i] eq 0 then begin
                 ser2Y0[j] = -1
              endif 
           end
           19:begin
              if eno[i] eq 0 then begin
                 Ie2[j] = -1
              endif 
           end
           20:begin
              if eno[i] eq 0 then begin
                 Re2[j] = -1
              endif  
           end
           21:begin
              if eno[i] eq 0 then begin
                 n2[j] = -1
              endif 
           end
           22:begin
              if eno[i] eq 0 then begin
                 ser2E[j] = -1
              endif
           end
           23:begin
              if eno[i] eq 0 then begin
                 ser2PA[j] = -1
              endif  
           end
           24:begin
              if eno[i] eq 0 then begin
                 ser3X0[j] = -1
              endif 
           end
           25:begin
              if eno[i] eq 0 then begin
                 ser3Y0[j] = -1
              endif 
           end
           26:begin
              if eno[i] eq 0 then begin
                 Ie3[j] = -1
              endif 
           end
           27:begin
              if eno[i] eq 0 then begin
                 Re3[j] = -1
              endif  
           end
           28:begin
              if eno[i] eq 0 then begin
                 n3[j] = -1
              endif 
           end
           29:begin
              if eno[i] eq 0 then begin
                 ser3E[j] = -1
              endif
           end
           30:begin
              if eno[i] eq 0 then begin
                 ser3PA[j] = -1
              endif  
           end
           31:begin
              if eno[i] eq 0 then begin
                 ser4X0[j] = -1
              endif 
           end
           32:begin
              if eno[i] eq 0 then begin
                 ser4Y0[j] = -1
              endif 
           end
           33:begin
              if eno[i] eq 0 then begin
                 Ie4[j] = -1
              endif 
           end
           34:begin
              if eno[i] eq 0 then begin
                 Re4[j] = -1
              endif  
           end
           35:begin
              if eno[i] eq 0 then begin
                 n4[j] = -1
              endif 
           end
           36:begin
              if eno[i] eq 0 then begin
                 ser4E[j] = -1
              endif
           end
           37:begin
              if eno[i] eq 0 then begin
                 ser4PA[j] = -1
              endif  
           end
           38:begin
              if eno[i] eq 0 then begin
                 bexpX0[j] = -1
              endif 
           end
           39:begin
              if eno[i] eq 0 then begin
                 bexpY0[j] = -1
              endif  
           end
           40:begin
              if eno[i] eq 0 then begin
                 bI0[j] = -1
              endif  
           end
           41:begin
              if eno[i] eq 0 then begin
                 bh1[j] = -1
              endif  
           end
           42:begin
              if eno[i] eq 0 then begin
                 bh2[j] = -1
              endif  
           end
           43:begin
              if eno[i] eq 0 then begin
                 bexpa[j] = -1
              endif 
           end
           44:begin
              if eno[i] eq 0 then begin
                 rb[j] = -1
              endif 
           end
           45:begin
              if eno[i] eq 0 then begin
                 bexpE[j] = -1
              endif 
           end
           46:begin
              if eno[i] eq 0 then begin
                 bexpPA[j] = -1
              endif 
           end
           47:begin
              if eno[i] eq 0 then begin
                 ferX0[j] = -1
              endif 
           end
           48:begin
              if eno[i] eq 0 then begin
                 ferY0[j] = -1
              endif 
           end
           49:begin
              if eno[i] eq 0 then begin
                 ferI0[j] = -1
              endif 
           end
           50:begin
              if eno[i] eq 0 then begin
                 l[j] = -1
              endif  
           end
           51:begin
              if eno[i] eq 0 then begin
                 fern[j] = -1
              endif 
           end
           52:begin
              if eno[i] eq 0 then begin
                 ferE[j] = -1
              endif
           end
           53:begin
              if eno[i] eq 0 then begin
                 ferPA[j] = -1
              endif  
           end
           54:begin
              if eno[i] eq 0 then begin
                 ferc[j] = -1
              endif
           end
           55:begin
              if eno[i] eq 0 then begin
                 fer2X0[j] = -1
              endif 
           end
           56:begin
              if eno[i] eq 0 then begin
                 fer2Y0[j] = -1
              endif 
           end
           57:begin
              if eno[i] eq 0 then begin
                 fer2I0[j] = -1
              endif 
           end
           58:begin
              if eno[i] eq 0 then begin
                 l2[j] = -1
              endif  
           end
           59:begin
              if eno[i] eq 0 then begin
                 fer2n[j] = -1
              endif 
           end
           60:begin
              if eno[i] eq 0 then begin
                 fer2E[j] = -1
              endif
           end
           61:begin
              if eno[i] eq 0 then begin
                 fer2PA[j] = -1
              endif  
           end
           62:begin
              if eno[i] eq 0 then begin
                 fer2c[j] = -1
              endif
           end
        endcase
     endfor
     printf,2,NAME[j],seeing[j],seeing2[j],betaa[j],$
            serX0[j], serY0[j], Ie[j], Re[j], n[j], serE[j], serPA[j],$
            expX0[j], expY0[j], I0[j], h[j], expE[j], expPA[j],$
            ser2X0[j], ser2Y0[j], Ie2[j], Re2[j], n2[j], ser2E[j], ser2PA[j],$
            ser3X0[j], ser3Y0[j], Ie3[j], Re3[j], n3[j], ser3E[j], ser3PA[j],$
            ser4X0[j], ser4Y0[j], Ie4[j], Re4[j], n4[j], ser4E[j], ser4PA[j],$
            bexpX0[j], bexpY0[j], bI0[j], bh1[j],bh2[j], bexpa[j], rb[j], bexpE[j], bexpPA[j],$
            ferX0[j], ferY0[j], ferI0[j], l[j], fern[j], ferE[j], ferPA[j], ferc[j],$
            fer2X0[j], fer2Y0[j], fer2I0[j], l2[j], fer2n[j], fer2E[j], fer2PA[j], fer2c[j]
            
  endfor
  
  close,2

  
end
