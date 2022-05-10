function pas = beadPastiche(roiImages,sel,scl)
if nargin>=2
    selImgs = roiImages(sel);
else
    selImgs = roiImages;
end
if nargin<3, scl=0; end
% calculate grid size (grd x grd) and separation
n   = numel(selImgs);
grd = ceil(sqrt(n));
sep = max(cellfun(@(c) max(size(c{1})),selImgs));
pas = zeros(grd*sep);
for ix=1:n
    ix1 = 1+sep*rem((ix-1),grd);
    ix2 = 1+sep*floor((ix-1)/grd);
    img = selImgs{ix}{1};
    if scl, img = img*128/median(img(:))/scl; end
    sz  = size(img)-1;
    pas(ix1+(0:sz(1)),ix2+(0:sz(2))) = img;
end
    
