function [BlockLine,TwolineConfocal,InitialPhase1,xishu,a1,b1,muban,ave] = LiOfSeAdjIndex(LineOriginal,LineZongShu,XiaFeng,TPixel,kA,Blockm,SubBlockCount,levelOffset)
BlockLine = double(LineOriginal((SubBlockCount-1)*Blockm+1:SubBlockCount*Blockm,:,(1+(LineZongShu-XiaFeng)/2):(XiaFeng+(LineZongShu-XiaFeng)/2)));
%The number of strips in BlockLine is XiaFeng
TwolineConfocal = zeros(Blockm,Blockm,TPixel);
for VSMPhaseCounter = 1:TPixel
    for i = 1:XiaFeng
        if mod(i-XiaFeng/2+VSMPhaseCounter-1,TPixel) < TPixel/2
            TwolineConfocal(:,:,VSMPhaseCounter) = TwolineConfocal(:,:,VSMPhaseCounter)+BlockLine(:,:,i);
        end
    end
end
ave = sum(TwolineConfocal(:,:,:),3)/XiaFeng;
level = (max(max(ave))+min(min(ave)))/2/max(max(ave))+levelOffset;
muban  = imbinarize(ave/max(max(ave)),level);
%figure()
%imshow(muban);
xishu(1) = sum(sum(TwolineConfocal(:,:,TPixel/4*0+1).*muban))-sum(sum(ave.*muban));
xishu(2) = sum(sum(TwolineConfocal(:,:,TPixel/4*1+1).*muban))-sum(sum(ave.*muban));
xishu(3) = sum(sum(TwolineConfocal(:,:,TPixel/4*2+1).*muban))-sum(sum(ave.*muban));
xishu(4) = sum(sum(TwolineConfocal(:,:,TPixel/4*3+1).*muban))-sum(sum(ave.*muban));
xishu  = xishu/max(xishu);
InitialPhase1 = atan2(-xishu(2),xishu(1));
a1 = cos(InitialPhase1);
b1 = sin(InitialPhase1);
