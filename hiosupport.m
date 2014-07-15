function Sup = hiosupport(image_size, support_size)

% image dimension
img_h = image_size(1);
if numel(image_size) > 1
    img_w = image_size(2);
else
    img_w = image_size;
end

% support dimenstion
sp_h = support_size(1);
if numel(support_size) > 1
    sp_w = support_size(2);
else
    sp_w = support_size;
end



Sup = false(img_h, img_w);
Sup(1:sp_h, 1:sp_w) = true;

Sup = circshift(Sup, [round((img_h-sp_h)/2) round((img_w-sp_w)/2)]);
