function [dd, maxNCC,RangeNCC] = disp_est(I1, I2, num_win, len_win,lat_len_win, Range)
% Copyright Hassan Rivaz 2006-2016 all rights reserved
[mm,nn] = size(I1); %mm is much bigger than nn

dis_win = round(mm/num_win);
minNCC = 0.3;
dd = zeros (num_win,nn);
maxNCC = zeros (num_win,nn);
RangeNCC = zeros (2*Range + 1,1);
lat_len=(lat_len_win-1)/2;

% for k = 25:nn-lat_len-1-Rangel-6

for k = 1:nn
    
    prevDispa = 0;
    
  
    
    if lat_len+1<=k  & k<= nn-lat_len
       RF1 = I1(:,k-lat_len:k+lat_len);
      RF2 = I2(:,k-lat_len:k+lat_len);
    elseif k<lat_len+1
       RF1 = I1(:,1:k+lat_len);
    RF2 = I2(:,1:k+lat_len);
    else
       RF1 = I1(:,k-lat_len:nn);
   RF2 = I2(:,k-lat_len:nn);
    end
    
     for i = 1:num_win
                  
             if (k == 1) || (maxNCC(i,k-1)<minNCC)
                u_max  = prevDispa - Range; % left max search range
                d_max = prevDispa + Range; % right max search range
            else
                u_max  = min(prevDispa,round(dd(i,k-1))) - Range;
                d_max = max(prevDispa,round(dd(i,k-1))) + Range;
            end
         
             if ((len_win+(i-1)*dis_win+ Range+d_max+1)>mm) || (len_win+(i-1)*dis_win+ Range +1)>mm
                continue;    
             end
            
             
            for j = u_max:d_max

                RangeNCC(j-u_max+1,1) = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+j,len_win);
                                
                if RangeNCC(j-u_max+1,1) >= maxNCC(i,k)
     				maxNCC(i,k) = RangeNCC(j-u_max+1,1);
                    dd(i,k) = j;
                  
                end
            end
   
            %calculating subpixel displ 
              posa = dd(i,k);
            
%  
            
%%            
            if (maxNCC(i,k) >= minNCC) && (dd(i,k) > u_max) && (dd(i,k) < d_max) % max happens internally

                 if (posa-u_max) >= 1
                    a1 = RangeNCC(posa-u_max,1);
                 else
                    a1 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa-1),len_win);
                 end
                
                b = RangeNCC(posa-u_max+1,1);
                
              
                if (d_max-posa) >= 1  
                    c1 = RangeNCC(posa-u_max+2,1);
               else
                    c1 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa+1),len_win);
               end           
                delta_d = sbupixel(a1,b,c1,1);
                dd(i,k) = (posa+delta_d);
                
                prevDispa = round(posa+delta_d); 

%                     dd(i,k) = sbupixel2(u_max,d_max,RangeNCC);
%                    prevDispa = round(dd(i,k)); 

            else
                dd(i,k) = prevDispa;
            end
            
            
            
%%             
%              if (maxNCC(i,k) >= minNCC) && (dl(i,k) > l_max) && (dl(i,k) < r_max) % max happens internally
% 
% %                  if (posl-l_max) >= 1
%                      a = RangeNCC(posa-u_max+1,posl-l_max);
% %                  else
% %                    RF2 = I2(:,k+posl-lat_len-1:k+posl+lat_len-1); 
% %                      a = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
% %                   end;
%                 
%                 b = RangeNCC(posa-u_max+1,posl-l_max+1);
%                 
%                
% %                  if (r_max-posl) >= 1  
%                      c = RangeNCC(posa-u_max+1,posl-l_max+2);
% %                  else
% %                      RF2 = I2(:,k+posl-lat_len+1:k+posl+lat_len+1); 
% %                     c = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
% %                  end;
%             
%                 delta_dl = sbupixel(a,b,c,1);
%                 dl(i,k) = (posl+delta_dl);
%                 prevDispl = round(posl+delta_dl); 
%              else
%                  dl(i,k) = prevDispl;
%              end
%% spline interp.
           
%               if (maxNCC(i,k) >= minNCC) && (dl(i,k) > l_max) && (dl(i,k) < r_max) % max happens internally
%                     
%                   RF2 = I2(:,k+posl-lat_len-2:k+posl+lat_len-2); 
%                   a2 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
%                 
%                   RF2 = I2(:,k+posl-lat_len-1:k+posl+lat_len-1); 
%                   a1 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
%                             
%                   b = RangeNCC(posa-u_max+1,posl-l_max+1);
%                 
%                   RF2 = I2(:,k+posl-lat_len+1:k+posl+lat_len+1); 
%                   c1 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
%                
%                   RF2 = I2(:,k+posl-lat_len+2:k+posl+lat_len+2); 
%                   c2 = NCC(RF1,RF2,(i-1)*dis_win+ Range+1,(i-1)*dis_win+Range+1+(posa),len_win);
% %                 delta_dl = sbupixel(a1,b,c1,1);
%                  delta_dl = sbupixell(a2,a1,b,c1,c2);
%                 dl(i,k) = (posl+delta_dl);
%                 prevDispl = round(posl+delta_dl); 
%              else
%                  dl(i,k) = prevDispl;
%              end;
            
%%            
     end
end  

   