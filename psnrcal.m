clear
clc
f_sample = 40e6; %  Sampling frequency [Hz]
c_sound = 1540e3;
load('Im0_with_tstart')
I1 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im1_with_tstart')
I2 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im2_with_tstart')
I3 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im3_with_tstart')
I4 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im4_with_tstart')
I5 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im5_with_tstart')
I6 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im6_with_tstart')
I7 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );


load('Im0p_with_tstart')
I8 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im1p_with_tstart')
I9 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im2p_with_tstart')
I10 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im3p_with_tstart')
I11 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im4p_with_tstart')
I12 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im5p_with_tstart')
I13 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );
load('Im6p_with_tstart')
I14 = Im0(5*2*round(f_sample/c_sound):25*2*round(f_sample/c_sound) , : );

maxI = max(I1(:));
[a,b]=size(I1);

I1 = (I1/maxI);
I2 = (I2/maxI);
I3 = (I3/maxI);
I4 = (I4/maxI);
I5 = (I5/maxI);
I6 = (I6/maxI);
I7 = (I7/maxI);

I8 = (I8/maxI);
I9 = (I9/maxI);
I10 = (I10/maxI);
I11 = (I11/maxI);
I12 = (I12/maxI);
I13 = (I13/maxI);
I14 = (I14/maxI);

noise1=0.7*rand(a,b);
noise2=0.3*rand(a,b);
noise3=0.3*rand(a,b);
noise4=0.3*rand(a,b);
noise5=0.3*rand(a,b);
noise6=0.3*rand(a,b);
noise7=0.3*rand(a,b);
noise8=0.3*rand(a,b);
noise9=0.3*rand(a,b);
noise10=0.3*rand(a,b);
noise11=0.3*rand(a,b);
noise12=0.3*rand(a,b);
noise13=0.3*rand(a,b);
noise14=0.3*rand(a,b);

I15 = I1+noise1;
I16 = I2;
I17 = I3;
I18 = I4;
I19 = I5;
I20 = I6;
I21 = I7+noise7;

I22 = I8;
I23 = I9;
I24 = I10;
I25 = I11;
I26 = I12;
I27 = I13;
I28 = I14;
peaksnr = psnr(I15,I1)