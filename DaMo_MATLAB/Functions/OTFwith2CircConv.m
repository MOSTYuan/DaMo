function [Out] = OTFwith2CircConv(NA,lambda,xy_num,dxy,Choose)
dp = 1/xy_num/dxy;
Kdx = lambda/2/NA;
Kdp = 2*NA/lambda;
KApodp = 2*Kdp;
[X,Y] = meshgrid(-(xy_num/2-1):xy_num/2,-(xy_num/2-1):xy_num/2);
PRadius = sqrt(X.^2+Y.^2)*dp;
OTF = (2*acos(PRadius/Kdp)-sin(2*acos(PRadius/Kdp)))/pi;
OTF(PRadius>=Kdp) = 0;
OTFApo = (2*acos(PRadius/KApodp)-sin(2*acos(PRadius/KApodp)))/pi;
OTFApo(PRadius>=KApodp) = 0;
switch(Choose)
    case('OTF')
        Out = OTF;
    case('OTFApo')
        Out = OTFApo;
end
end

% figure()
% imshow(abs(OTF));