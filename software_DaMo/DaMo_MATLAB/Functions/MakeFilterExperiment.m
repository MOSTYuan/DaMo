function [ Filter ] = MakeFilterExperiment(Blockm,IsPadImage,OTF,kA,FastSpReIsFilter,WienerType,IsFmask,IsApodize,IsAttOTF)
if IsPadImage == 1
    OTF_Pad = padarray(OTF,[Blockm/2,Blockm/2]);
    BlockmPad = Blockm*2;
else
    OTF_Pad = OTF;
    BlockmPad = Blockm;
end
wo = BlockmPad/2;
if FastSpReIsFilter == 1
    SFo = 1;    %modulation depth
    SignalPower = 1;
    NoisePower = WienerType.Co;
    if IsAttOTF.AMP > 0
        [X,Y] = meshgrid(1:BlockmPad,1:BlockmPad);
        Cv = (X-wo) + 1i*(Y-wo);
        R = abs(Cv);
        cyclesPerMicron = 1/(wo*(IsAttOTF.PixelSize));
        att_fun= 1-(IsAttOTF.AMP-0.00001)*exp(-R.^2/(2*(IsAttOTF.Vari/cyclesPerMicron)^2));
        att_OTF = att_fun.*OTF_Pad;
    else
        att_OTF = OTF_Pad;
    end
    MoveMax = ceil(max(abs(kA)));
    OTFoPad = padarray(OTF_Pad,[MoveMax,MoveMax]);
    OTFoPad = circshift(OTFoPad,round(-kA));
    OTFoShifted = OTFoPad(MoveMax+1:MoveMax+BlockmPad,MoveMax+1:MoveMax+BlockmPad);
    att_OTFPad = padarray(att_OTF,[MoveMax,MoveMax]);
    att_OTFPad = circshift(att_OTFPad,round(-kA));
    att_OTFShifted = att_OTFPad(MoveMax+1:MoveMax+BlockmPad,MoveMax+1:MoveMax+BlockmPad);
    DeFilter = SFo.*conj(att_OTFShifted)./((SFo.^2).*OTFoShifted.*conj(att_OTFShifted) + NoisePower./SignalPower);
else
    DeFilter = ones(BlockmPad,BlockmPad);
end

if IsFmask == 1
    Fmask = zeros(BlockmPad,BlockmPad);
    if kA(1)<0
        for i = ceil(BlockmPad/2):BlockmPad
            Fmask(i,:) = 1;
        end
        Fmask(ceil(BlockmPad/2),ceil(BlockmPad/2)) = 1/2;
        Fmask(ceil(BlockmPad/2),1:ceil(BlockmPad/2)-1) = 0;
    else
        for i = 1:floor(BlockmPad/2)
            Fmask(i,:) = 1;
        end
        Fmask(floor(BlockmPad/2),floor(BlockmPad/2)) = 1/2;
        Fmask(floor(BlockmPad/2),floor(BlockmPad/2)+1:BlockmPad) = 0;
    end
else
    Fmask = ones(BlockmPad,BlockmPad);
end

if IsApodize ~= 0
    nfac = 2;
    [X] = meshgrid(1:BlockmPad);
    x_im = X/(BlockmPad);
    if IsApodize == 1
        FFTApomask = sin(pi*nfac/2*x_im);  
    elseif IsApodize == 2
        FFTApomask = 0.5-0.5*cos(pi*nfac*x_im);      
    end
    FFTApomask(x_im>1/nfac) = 1;
    FFTApomask = FFTApomask.*rot90(FFTApomask);
    FFTApomask = FFTApomask.*rot90(FFTApomask,2);
else
    FFTApomask = ones(BlockmPad,BlockmPad);
end

% Low frequency suppression
if IsAttOTF.AttLowFreq ~= 0
    wo = BlockmPad/2;
    [X,Y] = meshgrid(1:BlockmPad,1:BlockmPad);
    Cv = (X-wo) + 1i*(Y-wo);
    R = abs(Cv);
    cyclesPerMicron = 1/(wo*(IsAttOTF.PixelSize));
    att_LowFreq = 1-(IsAttOTF.AttLowFreq)*exp(-R.^2/(2*(IsAttOTF.Vari/cyclesPerMicron)^2));
else
    att_LowFreq = ones(BlockmPad,BlockmPad);
end

Filter = DeFilter.*Fmask(:,:).*FFTApomask.*att_LowFreq;

