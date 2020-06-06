function bb= sbupixel2(u_max,d_max,RangeNCC)
ab=RangeNCC';
x=u_max:d_max;
gam=1;
xq=u_max:gam:d_max;
b=spline(x,ab,xq);
d=find(b==max(b));
bb=u_max+(d-1)*gam;

d=2;
