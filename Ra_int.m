function [I3] = Ra_int(Im, x0, y0)
[a1 b1]=size(x0);
re_x0=reshape(x0,a1*b1,1);
re_y0=reshape(y0,a1*b1,1);
eIm1=Im(sub2ind(size(Im),floor(re_x0'+1),floor(re_y0'+1)));
Im1=reshape(eIm1',size(Im));
eIm2=Im(sub2ind(size(Im),floor(re_x0'+1),floor(re_y0')));
Im2=reshape(eIm2',size(Im));
eIm3=Im(sub2ind(size(Im),floor(re_x0'),floor(re_y0'+1)));
Im3=reshape(eIm3',size(Im));
eIm4=Im(sub2ind(size(Im),floor(re_x0'),floor(re_y0')));
Im4=reshape(eIm4',size(Im));
I3 = (x0 - floor(x0) ).* ((y0 - floor(y0) ).* Im1 + (floor(y0+1)-y0).* Im2) + (floor(x0+1) - x0 ).* ((y0 - floor(y0) ).* Im3 + (floor(y0)+1-y0).* Im4);