;+
; NAME:
;    outer_sub
;
; PURPOSE:
;    The outer subtraction of the arrays X and Y is the array A with
;    dimension c(dim(X), dim(Y)) where element A[c(arrayindex.x, 
;    arrayindex.y)] = FUN(X[arrayindex.x], Y[arrayindex.y], ...).
;    
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    d = outer(X, Y, FUN = "*")
;
; INPUTS:
;    x, y: First and second vectors for function FUN.
;
;    FUN: a function to use on the outer products, found via match.fun (except for the special case "*")
;
;
; OUTPUTS:
;    d: probability density estimated at each value specified by t
;


FUNCTION outer_sub, X, Y
  
  xsize = size(X, /DIMENSIONS)
  ysize = size(Y, /DIMENSIONS)
  
  x = double(x)
  y = double(y) 
  
  stop
  z = make_array(xsize[0],ysize[0],/DOUBLE)
  
  
  for ix = 0l, xsize[0]-1 do begin
     for iy = 0l, ysize[0]-1 do z[ix,iy] = X[ix] - Y[iy]
  endfor
  
  
  return, z
END
