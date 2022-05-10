function biOutput = beadImageAnalysis(roiImages,arg2)
% invoke as 
%   bia = beadImageAnalysis(roiImages)  % returns all output
%   bia = beadImageAnalysis(roiImages,sel);
%   bia = beadImageAnalysis(roiImages,trainedModel);  % calculates .cs
global beaddata
if nargin>1 && isnumeric(arg2)
    selImg = roiImages(arg2);
else
    selImg = roiImages;
end
dispImg = false;
for k=1:numel(selImg)
    dat = selImg{k};
    img = double(dat{1});
    sz  = size(img);
    ctr = dat{2};
    radius = dat{3};
    [xxL,yyL] = meshgrid(1:sz(2),1:sz(1));
    r = hypot(xxL-ctr(1),yyL-ctr(2));
    msk = r < radius;
    if dispImg
        figure(44); clf;
        subplot(2,1,1); image(img*64/median(img(:)))
        viscircles(ctr,radius);
        subplot(2,1,2); imagesc(double(msk));
    end
    % calculate the contained intensity values and their radial position
    vals = img(msk);
    rads = r(msk);
    [ordrads,idx] = sort(rads);
    ordvals = vals(idx);
    nhalf   = round(numel(vals)/2);
    n90     = find(ordrads>radius*.9,1,'first');
    med  = median(vals);
    biOutput(k,:) = [median(vals) iqr(vals) skewness(vals) kurtosis(vals) ...
        sum(vals>(med*3))/numel(vals) mean(ordvals(1:nhalf))/mean(ordvals(nhalf+1:end))];
    if dispImg
        figure(45); clf; plot(ordrads(1:n90),ordvals(1:n90)); % hist(vals,40);
        annotation('textbox',[.6 .6 .3 .3],'String',sprintf('%0.1f ',biOutput(k,:)),'FontSize',8);
    end
end
% if a model is included, calculate the classifier
if nargin>1 && ~isnumeric(arg2)
    beaddata.cs = arg2.predictFcn(biOutput);
end

