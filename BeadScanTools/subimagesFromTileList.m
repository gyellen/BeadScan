function newsubs = subimagesFromTileList(tileIJList,xmlFilename)
% identify the tiles of interest, by tileIJ [col row]
tinfo = readThorimageExperimentFile(xmlFilename);
% is there a single master Tile array?
if numel(tinfo.tiles)~=1
    error('Did not find single master tile array in experiment');
end
% create the new subimages
XY0 = tinfo.tile0xy/1000 - tinfo.homeOffset; % for transOffsetXMM,YMM 
dXY = tinfo.tileDxy/1000;
ovr = tinfo.tileOvr;
newsubs = {};
for k=1:size(tileIJList,1)
    IJ = tileIJList(k,:);
    XY = XY0 + dXY .* (IJ - [1 1]);
    str = ['      <SubImages ' ...
        'name="c' num2str(IJ(1)) 'r' num2str(IJ(2)) '" ' ...
        'isEnabled="True" subRows="1" subColumns="1" ' ...
        'transOffsetXMM="' num2str(XY(1),'%.4f') '" ' ...
        'transOffsetYMM="' num2str(XY(2),'%.4f') '" ' ...
        'transOffsetZMM="0.0000" ' ...
        'overlapX="' num2str(ovr(1)) '" overlapY="' num2str(ovr(2)) ...
        '" />'];
    newsubs{k}=str;
end