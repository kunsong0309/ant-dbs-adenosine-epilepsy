function [L, A] = imMaskDAPI(im, xfold)
if nargin < 2
    xfold = 40;
end

im = double(im);

se = strel('disk', round(10 / 40 * xfold));
M = imopen(imgaussfilt(im, 5 / 40 * xfold), se);
% figure; imshow(M/1000, [0, 1]);

fudgeFactor = 1;
method = 'canny';
[~,threshold] = edge(M, method);
BWs = edge(M, method, threshold * fudgeFactor);
% figure; imshow(BWs);

se360 = strel('disk', round(3 / 40 * xfold));
BWsdil = imdilate(BWs, se360);
% figure; imshow(BWsdil);

BWdfill = imfill(BWsdil, 'holes');
% figure; imshow(BWdfill);

BWfinal = BWdfill & ~BWsdil;
seD = strel('diamond', 1);
BWfinal = imerode(BWfinal, seD);
BWfinal = imerode(BWfinal, seD);
% figure; imshow(labeloverlay(I/1000, BWfinal, 'Transparency', 0.5));

[B,L] = bwboundaries(BWfinal,'noholes');

amin = 100 / ((40 / xfold) .^ 2);
Lmsk = L;
A = nan(size(B));
isv = true(size(B));
for kk = 1:length(B)
    A(kk) = bwarea(L == kk);
    if A(kk) < amin
        Lmsk(L == kk) = 0;
        isv(kk) = false;
    end
end
% figure; imshow(labeloverlay(im, Lmsk, 'Transparency', 0.5));

[~, L] = bwboundaries(Lmsk,'noholes');
A = A(isv);

end