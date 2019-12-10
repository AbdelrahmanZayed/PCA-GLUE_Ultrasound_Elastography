% Copyright Hassan Rivaz 2006-2016 all rights reserved
function delta_d = sbupixel(a,b,c,n)

if n==1 % parabolic interpolation
    denum = -a+2*b-c;
    if denum
        delta_d = 0.5*(c-a)/denum;
    else
        delta_d = 0;
    end;
else % cosine interpolation
    omega = acos((a+c)/(2*b));
    teta = atan((a-c)/(2*b*sin(omega)));
    delta_d = -teta/omega;
end;
