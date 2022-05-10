function acqorder = tileAcqorderFromIJ(ij,tilesXY)
% given that tiles are acquired in a serpentine pattern from top left, with 
%   tilesXY(1) columns and tilesXY(2) rows, convert serial number to column
%   and row.  Numbering of tiles is NOT serpentine, only their order of acq
% note that I is column number and J is row number
j = ij(2) - 1;  % row number (0-based)
i = ij(1);  % column number if not serpentine (0-based)
if floor(j/2)~=j/2 % odd numbered rows are scanned backwards
    i = 1+tilesXY(1)-i;
end
acqorder = i+j*tilesXY(1);

    