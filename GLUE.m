% GLUE        Calculates the displacement between two RF frames
%       Please cite the following paper if you are using this code:
%       Hashemi, H., Rivaz, H., Global Time-Delay Estimation in Ultrasound
%       Elastography, IEEE Trans. UFFC, 2017
%
%       Regularization coefficients:
%       alfa1 = axial continuity in a A-line
%       alfa2 = axial continuity at the neighboring samples
%       beta1 = lateral continuity in a A-line
%       beta2 = lateral continuity at the neighboring samples
%       The best set is 5, 1, 5, 1 or if you need smoother displacements,
%       put them 20, 1, 20, 1 respectively.
%
% Example:
%       GLUE (initialAxialDisp, initialLateralDisp, Im1, Im2, 20, 1, 20, 1) 
%
%
% Copyright 2015-2016 - Hoda S. Hashemi (c) Aug 2016



function [ax_d, lat_d] = GLUE (ax0, lat0, Im1, Im2, alfa1, alfa2, beta1, beta2)

[m, n]=size(Im1);

a = ax0;
l = lat0;

for i=1:n
    a(:, i) = step2linear( a(:, i), m ); % interpolation of DP
    l(:, i) = step2linear( l(:, i), m );
end
[Gx,Gy] = imgradientxy(Im2,'centraldifference'); % X axis: lateral, Y axis: axial
V = Gy; % axial derivative
W = Gx; % lateral derivative

dp_axial = reshape(a',m*n,1); % putting all the axial disp. in one vertical vector
dp_lateral = reshape(l',m*n,1); % putting all the lateral disp. in one vertical vector

% merging two vectors in order to obtain that SPECIFIED vector for disp:
y = [dp_axial'; dp_lateral'];
dp = reshape(y, 2*m*n, 1);

% we need derivative of Im2 w.r.t. a(i,j) at the point (i+a(i,j), j+l(i,j)):
Der2_a = zeros(m, n);
Der2_l = zeros(m, n);
DiffI = zeros(m, n);
a_calc=zeros(size(a));
a_calc(1:m-1,1:n-1)=a(1:m-1,1:n-1);
l_calc=zeros(size(l));
l_calc(1:m-1,1:n-1)=l(1:m-1,1:n-1);


[row1 col1]=find(a_calc==a_calc);
re_row1=reshape(row1,m,n);
re_col1=reshape(col1,m,n);
X1 = re_row1+a_calc;
Y1 = re_col1+l_calc;
X1(X1 >= m-1)=m-1;
Y1(Y1>= n-1)=n-1;
X1(X1 <= 1)=1;
Y1(Y1 <= 1)=1;
Der2_a = Ra_int(V, X1, Y1) ; % axial der. @ point X1, Y1
Der2_l = Ra_int(W, X1, Y1) ; % lateral der. @ point X1, Y1
DiffI = Im1- Ra_int(Im2, X1, Y1);


Der2_a(:,end)=0;
Der2_a(end,:)=0;
Der2_l(:,end)=0;
Der2_l(end,:)=0;
DiffI(end,:)=0;
DiffI(:,end)=0;


DiffI_Vec = reshape(DiffI', m*n, 1);
t = [DiffI_Vec'; DiffI_Vec'];
DiffI_Vec2 = reshape(t, 2*m*n, 1);

der2_axial= reshape(Der2_a', m*n, 1); % putting all the axial der. of Im2 in one vertical vector
der2_lateral= reshape(Der2_l', m*n, 1); % putting all the lateral der. of Im2 in one vertical vector
% merging two vectors in order to obtain that SPECIFIED vector for der:
y = [der2_axial'; der2_lateral'];
Der2 = reshape(y, 2*m*n, 1);
aa = reshape ([(der2_axial.^2)'; (der2_lateral.^2)'], 2*m*n,1);
usum2 = sparse(1:2*m*n,1:2*m*n,aa',2*m*n,2*m*n);
mm = der2_axial.*der2_lateral ;
mm = reshape ([mm'; zeros(1,m*n)], 2*m*n,1);
mm(2*m*n) = [];
mmUP1diag = sparse(2:2*m*n,1:2*m*n-1,mm',2*m*n,2*m*n);
usum2 = usum2 + mmUP1diag + mmUP1diag'; % 2nd and 3rd terms are related to upper diagonal elements
MainDiag=zeros(1,2*n*m);
% First 2n elements on the diagonal:
MainDiag(1) = alfa1 + alfa2;
MainDiag(2) = beta1 + beta2;
MainDiag(2*n-1) = alfa1 + alfa2;
MainDiag(2*n) = beta1 + beta2;


MainDiag(3:2:2*n-3) = alfa1 + 2*alfa2;
MainDiag(4:2:2*n-2) = beta1 + 2*beta2;


% last 2n elements on the diagonal:
MainDiag(2*m*n-2*n+1:2*m*n) = MainDiag(1:2*n);
% blue elements:
MainDiag(2*n + 1) = 2*alfa1+alfa2; % fisrt one
MainDiag(2*n +2) = 2*beta1+beta2; % second one
MainDiag(2*m*n-2*n-1) = 2*alfa1+alfa2; % before last one
MainDiag(2*m*n-2*n) = 2*beta1+beta2; % last one


MainDiag(2*n*(2:m-2) - 1) = 2*alfa1+alfa2;
MainDiag(2*n*(2:m-2)) = 2*beta1+beta2;
MainDiag(2*n*(2:m-2) + 1) = 2*alfa1+alfa2;
MainDiag(2*n*(2:m-2) + 2) = 2*beta1+beta2;


% black elements:
for i = 1 : m-2
    for j = 2*n*i+3 : 2 : 2*n*(i+1)-2-1
        MainDiag(j) = 2*alfa1 + 2*alfa2;
        MainDiag(j+1) = 2*beta1 + 2*beta2;
    end
end

% usum1 = diag(MainDiag);
usum1 = sparse(1:2*m*n,1:2*m*n,MainDiag,2*m*n,2*m*n);
D = usum1 + usum2;


% U: Upper triangular matrix
ue1=zeros(1, 2*n*(m-1));

ue1(1:2:2*n*(m-1))=-alfa1;
ue1(2:2:2*n*(m-1)+1)=-beta1;


% u1=diag(ue1, 2*n);
u1 = sparse(2*n+1:2*m*n,1:2*m*n-2*n,ue1,2*m*n,2*m*n);

ue2=zeros(1, 2*(n-1));
ue2(1:2:2*n-2)=-alfa2;
ue2(2:2:2*n-1)=-beta2;
u2 = [ue2, 0, 0];
usum3 = u2;
for k=1:m-1
    usum3 = [ usum3 , u2 ];
end
usum3(2*m*n) = [];
usum3(2*m*n - 1) = [];
% usum4 = diag( usum3 , 2 );
usum4 = sparse(3:2*m*n,1:2*m*n-2,usum3,2*m*n,2*m*n);
U = u1 + usum4; % upper tri. matrix
% L: lower triangular matrix
L=U';
% Ax=b
A= U+ D+ L;
% b:
Anew= usum1+ U + L;
b= -Anew*dp + Der2.*DiffI_Vec2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adding bias:
Ax_eps = 1/(m-1) *( a(m, :) - a(1, :) ) ;
Lat_eps = 1/(n-1) *( l(:, n) - l(:, 1) ) ;

if max(Ax_eps) >= 0
    eps1 = max(Ax_eps) ;
else
    eps1 = min(Ax_eps) ;
end

if max(Lat_eps) >= 0
    eps2 = max(Lat_eps) ;
else
    eps2 = min(Lat_eps) ;
end

bias = sparse(2*m*n , 1);


bias(2*n*(1:m-1)) = beta2*eps2 ;
bias(2*n*(1:m-1) - 1) = alfa2*eps1 ;
bias(2*m*n-2*n+1 : 2 : 2*m*n-3) = alfa1*eps1 ;
bias(2*m*n-2*n+2 : 2 : 2*m*n-2) = beta1*eps2 ;

bias(2*m*n-1) = alfa1*eps1 + alfa2*eps1 ;
bias(2*m*n) = beta1*eps2 + beta2*eps2 ;


for i= 1 : 2*n-3
    bias(i) = -alfa1*eps1 ;
    bias(i+1) = -beta1*eps2 ;
end


bias(2*n-1) = -alfa1*eps1 + alfa2*eps1 ;
bias(2*n) = -beta1*eps2 + beta2*eps2 ;

bias(2*n*(0:m-1)+1) = -alfa2*eps1 ;
bias(2*n*(0:m-1)+2) = -beta2*eps2 ;


bias(1) = (-alfa1-alfa2)*eps1 ;
bias(2) = (-beta1-beta2)*eps2 ;
bias(2*n*(m-1)+1) = (alfa1-alfa2)*eps1 ;
bias(2*n*(m-1)+2) = (beta1-beta2)*eps2 ;
b = b + bias ;
x = A\b ;


% converting x to the initial matrix shape:
axi=x(1:2:end-1);
late=x(2:2:end);
deltaD_a=(reshape(axi,n,m))';
deltaD_l=(reshape(late,n,m))';


ax_d =  a+ deltaD_a; %axial: d + delta d
lat_d = l+ deltaD_l; %lateral: d + delta d