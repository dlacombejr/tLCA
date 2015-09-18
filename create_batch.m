function I=create_batch(Images,patch_size,batch_size)

[~, imsize, num_Images] = size(Images);

border=10;
patch_side = sqrt(patch_size);

I = zeros(patch_size,batch_size);

im_num= ceil(num_Images * rand());

for i=1:batch_size
    
    row = border + ceil((imsize-patch_side-2*border) * rand());
    col = border + ceil((imsize-patch_side-2*border) * rand());
    
    I(:,i) = reshape(Images(row:row+patch_side-1, col:col+patch_side-1, im_num),[patch_size, 1]);
    
    
end


end