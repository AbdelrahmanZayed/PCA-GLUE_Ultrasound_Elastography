function normal_dot_prod = NCC (I1, I2, start1, start2, length)
% Copyright Hassan Rivaz 2006-2016 all rights reserved
% normal_dot_prod is the normalized dot product, or NCC
if (start1<=0 | start2<=0 | (round(start1) ~= start1) | (round(start2) ~= start2))
    normal_dot_prod = 0;
    return
end

[s1,s2]=size(I1(start1:start1+length,:));
I1=reshape(I1(start1:start1+length,:),s1*s2,1);
I2=reshape(I2(start2:start2+length,:),s1*s2,1);

if I1==0 & I2==0
    normal_dot_prod = 1; 
else

dot_prod = I1'*I2;
var1 = I1'*I1;
var2 = I2'*I2;

denom = sqrt(var1*var2);

if (denom == 0)
    normal_dot_prod = 0; 
else
    normal_dot_prod = (dot_prod)/denom;
end

end