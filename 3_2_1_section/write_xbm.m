function write_xbm(IMAGE, filename)

[height, width] = size(IMAGE);
IMAGE = uint8(IMAGE');
tmp = zeros(ceil(width/8)*8,height,'uint8');
tmp(1:width,1:height) = IMAGE;
tmp=reshape(tmp>0,8,height*ceil(width/8))';

data = zeros(height*ceil(width/8),1,'uint8');
for i = 1:8
  data = bitor(data, bitshift(uint8(tmp(:,i)), i-1));
end

%% Write file
fid = fopen(filename, 'wt');
fprintf(fid, '#define im_width %d\n', width);
fprintf(fid, '#define im_height %d\n', height);
fprintf(fid, 'static char im_bits[] = {\n');
fprintf(fid, [repmat('0x%02x,', [1 15]) '\n'], data(1:end-1));
fprintf(fid, '%02x\n};\n',data(end));
fclose(fid);
