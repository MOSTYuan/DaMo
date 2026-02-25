%% ---------------------------Single-Spectrum reconstruction (SSR) algorithm---------------------------%
%%********************************************Introduction*********************************************%
% SSR is a rapid super-resolution reconstruction algorithm designed for Digital array 
% Modulation microscopy (DaMo), allowing reconstruction from illumination modulated images using 
% only a single high-order spectrum. 
%%*****************************************************************************************************%
% It is a part of publication:
% Sijie Li et al. Three-dimensional super-resolution imaging with suppressed background via digital 
% array modulation microscopy, Nature Photonics.
%%****************************************System requirements******************************************%
% Our software is developed on Windows 10 / Windows 11 with MATLAB (MathWorks®). MATLAB of R2021b or a
% later version is preferred. As a lightweight technology, SSR necessitates no specific computational or
% storage resources from the underlying computer hardware.
%%*****************************************************************************************************%
% If you have any questions, please contact yuanj@hust.edu.cn for help.
%%*****************************************************************************************************%
% Copyright 2025 Sijie Li et al.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%*****************************************************************************************************%

close all
clear all
SSR_mode = 'HS';  %HC(High Contrast)/HF(High Fedility)/HS(High Speed);
Blockm = 2048;      %size of raw data
Blockn = 2048;  
PixelSize = 58;     %unit:nm
NA = 1.5;           % Objective numerical aperture
lambda = 0.52;      % Emission wavelength (μm)
mainFolder = fileparts(mfilename('fullpath'));  
addpath(fullfile(mainFolder, 'Functions'));
path=fullfile(mainFolder,'\Demo\Raw_data\');
pathout=fullfile(mainFolder,'\Demo\SSR_result\');
filename='tubulin';
switch(SSR_mode)
    case('HF')
        IsFmask = 1;            % Half-spectrum filtering enabled (1) or disabled (0)
        IsApodize = 1;          % Apodization method: 0=None, 1=Sine(0-π/2), 2=Cosine(0-π). Both suppress high-frequency noise (Cosine stronger)
        IsUseDeconvolution = 1; % Deconvolution enabled (1) or disabled (0)
        IsPadImage = 1;         % 2× upsampling enabled (1) or disabled (0)
        WienerType.Type = 2;    % Wiener filter method (Experimental: Use 2 only). 1=WF uses fitted signal power, 2=Signal:1, Noise:0.1, 3=Actual signal power, 4=High orders use fitted power. 
        WienerType.Co = 0.2; 
    case('HC')
        IsFmask = 1;   
        IsApodize = 1; 
        IsUseDeconvolution = 0; 
        IsPadImage = 0; 
        WienerType.Type = 2; 
        WienerType.Co = 0.2; 
    case('HS')
        IsFmask = 0;   
        IsApodize = 0; 
        IsUseDeconvolution = 0; 
        IsPadImage = 0; 
        WienerType.Type = 2; 
        WienerType.Co = 0.2; 
end
WienerType.MF = 1;      % 0: m=ModulationFactor; 1: m=ModulationFactor_RealFq
% Notch filtering
IsAttOTF.AMP = 0;       % Notch filter intensity in Wiener deconvolution (0=disable)
IsAttOTF.PixelSize = PixelSize;   
IsAttOTF.Vari = 1.5;    % Gaussian filter variance (/μm, typically 1.5)
IsAttOTF.AttLowFreq = 0.0;  % Low-frequency suppression amplitude (typically 0)
PhaseMatch = 1;
levelOffset = -0.1;
lineCount = 1;          % Stripe index for raw stripe size extraction
Subarray = 8;
T = 4;            % Slit width, matches number of input modulation pictures
kA = [-Blockm/T,0];
OTF = OTFwith2CircConv(NA,lambda,Blockm,PixelSize,'OTF');
OTF = OTF/max(OTF(:));  % Normalize OTF to [0,1] range
% Generating filters
[ Filter ] = MakeFilterExperiment(Blockm,IsPadImage,OTF,kA,IsUseDeconvolution,WienerType,IsFmask,IsApodize,IsAttOTF);

SignLX1 = -ones(T,2);
for i = 1:T
    if mod(i-T/2,T) < T/2
        SignLX1(i,1) = 1;
    end
    if mod(i-T/2+1,T) < T/2
        SignLX1(i,2) = 1;
    end
end
SignLX1 = 0.5*SignLX1;

if(~mkdir(pathout))
    mkdir(pathout);
end

a = imread(strcat(path,filename,'_',num2str(lineCount),'.tif'));   
[m,n] = size(a);
LineShifted = zeros(m+T,n,T,'double');
LineOriginal = zeros(m,n,T,'double');
IllIntens = zeros(1,T);
for i = 1:T
    LineOriginal(:,:,i) = double(imread(strcat(path,filename,'_',num2str(i),'.tif')))-80+10*(i==T/2||i==T/2+1);
    IllIntens(i) = sum(sum(LineOriginal(:,:,i).*(LineOriginal(:,:,1)>300)))/sum(sum((LineOriginal(:,:,1)>300)));       
    LineShifted(i:i+m-1,:,i) = LineOriginal(:,:,i);   
end
LineOriginal(LineOriginal < 0) = 0;
LineShifted(LineShifted < 0) = 0;
LineShifted = LineShifted(1:m,:,:);
m = floor(m/Blockm)*Blockm;
SubBlockNum = m/Blockm;
IRcorrectDeconvFinal = zeros(m*(1+IsPadImage),n*(1+IsPadImage),'double');
IRcorrectDeconvParFor = zeros(Blockm*(1+IsPadImage),Blockn*(1+IsPadImage),SubBlockNum,'double');

for SubBlockCount = 1:SubBlockNum   
    [BlockLine,TwolineConfocal,InitialPhase1,Xishu,a1,b1,muban,ave] = LiOfSeAdjIndex(LineOriginal,T,T,T,kA,Blockm,SubBlockCount,levelOffset);
    Index = zeros(1,T);
    if kA(1) < 0
        for i = 1:T
            Index(i) = (a1+1i*b1)*(SignLX1(i,1)+1i*(SignLX1(i,2)))*1;
        end
    else
        for i = 1:T
            Index(i) = (a1-1i*b1)*(SignLX1(i,1)-1i*(SignLX1(i,2)))*1;
        end
    end
    I1Shifted = zeros(Blockm,Blockn);
    for i = 1:T
        I1Shifted  = I1Shifted + Index(i).*BlockLine(:,:,i);
    end
    IF1Shifted = circshift( ifftshift( fft2( fftshift( circshift( I1Shifted , [1 1])))), [-1 -1] );
    if IsPadImage == 1
        IF1ShiftedPad = zeros(Blockm*2,Blockm*2);
        if kA(1) < 0
            IF1ShiftedPad(Blockm/2+1+Blockm/T:Blockm/2+Blockm,Blockm/2+1:Blockm/2+Blockm) = IF1Shifted(Blockm/T+1:Blockm,:);
            IF1ShiftedPad(Blockm*1.5+1:Blockm*1.5+Blockm/T,Blockm/2+1:Blockm/2+Blockm) = IF1Shifted(1:Blockm/T,:);
        else
            IF1ShiftedPad(Blockm/2+1:Blockm/2+Blockm-Blockm/T,Blockm/2+1:Blockm/2+Blockm) = IF1Shifted(1:Blockm-Blockm/T,:);
            IF1ShiftedPad(Blockm/2-Blockm/T+1:Blockm/2,Blockm/2+1:Blockm/2+Blockm) = IF1Shifted(Blockm-Blockm/T+1:Blockm,:);
        end
        OTF_Pad = padarray(OTF,[Blockm/2,Blockm/2]);
    else
        IF1ShiftedPad = IF1Shifted;
        OTF_Pad = OTF;
    end
    IFcorrectDeconv = IF1ShiftedPad.*Filter;
    if PhaseMatch == 1
        I0 = ave;
        I0F = circshift( ifftshift( fft2( fftshift( circshift( I0 , [1 1])))), [-1 -1] );
        if IsPadImage == 1
            I0F = padarray(I0F,[Blockm/2,Blockm/2]);
        end
        % phase matching
        [X,Y] = meshgrid(0:size(IFcorrectDeconv,1)-1,0:size(IFcorrectDeconv,1)-1);
        Cv = (X-size(IFcorrectDeconv,1)/2) + 1i*(Y-size(IFcorrectDeconv,1)/2);
        Ro = abs(Cv);
        kv = kA(2) + 1i*kA(1);
        RNega = abs(Cv+kv);
        k2 = sqrt(kA*kA');
        % frequency range over which corrective phase is determined
        Zmask = (Ro < 0.8*k2).*(RNega < 0.8*k2);
        % corrective phase
        Angle0 = angle( sum(sum( I0F.*conj(IFcorrectDeconv).*Zmask )) );
        % phase correction
        IFcorrectDeconv = exp(+1i*Angle0).*IFcorrectDeconv;
    end
    IRcorrectDeconv = circshift( fftshift( ifft2( ifftshift( circshift( IFcorrectDeconv , [1 1])))), [-1 -1] );
    IRcorrectDeconv = 2*real(IRcorrectDeconv);
    IRcorrectDeconv(IRcorrectDeconv<0) = 0;
    IRcorrectDeconv(1:Blockm-1,:) = IRcorrectDeconv(1+1:Blockm,:);
    IRcorrectDeconvParFor(:,:,SubBlockCount) = IRcorrectDeconv;
    for SubBlockCount = 1:SubBlockNum
        IRcorrectDeconvFinal((SubBlockCount-1)*Blockm*(1+IsPadImage)+1:SubBlockCount*Blockm*(1+IsPadImage),:) = IRcorrectDeconvParFor(:,:,SubBlockCount);
    end
    imwrite(uint16((IRcorrectDeconvFinal)),strcat(pathout,filename,'_',SSR_mode,'_SR.tif'));

end