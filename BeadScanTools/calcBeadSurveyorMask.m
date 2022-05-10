function mask0 = calcBeadSurveyorMask(xyr,dxy,xxL,yyL,szimg)
% [xxL,yyL] = meshgrid(-10:10,-10:10);  % local grid; center +/- 10 pixels

% calculate whole and fractional center
xyWhole = floor(dxy+xyr(1:2));
fCtr = dxy + xyr(1:2) - xyWhole;

% calculate a local mask, based on the fractional center
masklocal = hypot(yyL-fCtr(2),xxL-fCtr(1))<xyr(3);

% now create a full mask and copy the local mask into the correct subarray
mask0=zeros(fliplr(szimg));

% assign subarray, avoiding edge errors
lowlim = xyWhole - [10 10];
hilim  = xyWhole + [10 10];
lowlim = max(lowlim,[1 1]);
hilim  = min(hilim, szimg);
ix1 = lowlim(1):hilim(1);
ix2 = lowlim(2):hilim(2);

mask0(ix2,ix1) = masklocal(ix2-xyWhole(2)+11,ix1-xyWhole(1)+11);
