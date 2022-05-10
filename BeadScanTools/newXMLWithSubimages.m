function xmlOut = newXMLWithSubimages(oldXMLName,tileIJList)
% outputs a cell array with the new XML file
% sample invocation:
%  xmlOut = newXMLWithSubimages('experiment.xml', ...
%                 [1 1; 3 2; 11 15]);
% each row of the tile IJ list has [col row]
newSubs = subimagesFromTileList(tileIJList,oldXMLName);
ff = fopen(oldXMLName);
xmlOut = {};
lin = fgetl(ff);
% copy until SubImages
while ~contains(lin,'SubImages')
    xmlOut{end+1} = lin;
    lin = fgetl(ff);
end
% skip any unused
while ~contains(lin,'isEnabled="True"')
    lin = fgetl(ff);
end
pos = strfind(lin,'isEnabled="True"');
lin = [lin(1:pos+10) 'False' lin(pos+15:end)];
xmlOut{end+1} = lin;
% now add the new subimages
xmlOut = [xmlOut newSubs(:)'];
% and copy until the end (except unused SubImages)
ok = 1;
while ok
    lin = fgetl(ff);
    if lin==-1
        ok = 0;
    elseif ~contains(lin,'SubImages')
        xmlOut{end+1} = lin;
    end
end
xmlOut = xmlOut(:);  % easier to look at when a column
fclose(ff);

    
