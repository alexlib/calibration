function [a] = write2_gdf(data, filename)

MAGIC_int = int32(082991);

sz = size(data);
idl_sz = int32([length(sz), sz(2), sz(1), 5, sz(1)]);

fid = fopen(filename, 'w');
if fid == -1, error('Could not open file'), end;

count = fwrite(fid, MAGIC_int, 'int32');
if count ~= 1,  error('failed to write Magic number'), end;
    
count = fwrite(fid, idl_sz , 'int32'); 
if count < 4,  error('failed to write array info'), end;
    
count = fwrite(fid, data', 'float64');
% data_int=uint32(data(:,3:6));
% for i=1:size(data,1)
%             fwrite(fid, data(i,1:2)' , 'float32'); 
%             fwrite(fid, data_int(i,:)' , 'uint32'); 
% end
if count < length(data),  error('failed to write data'), end;
 
fclose(fid);
a=1;

end
    
    
    

