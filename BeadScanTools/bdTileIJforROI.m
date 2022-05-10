function ij = bdTileIJforROI(roi)
global beaddata
bdti  = beaddata.tileInfo;
hiROI = bdti(:,1)+bdti(:,2)-1;
row = find(roi>=bdti(:,1) & roi<=hiROI,1);
ij = bdti(row,3:4);