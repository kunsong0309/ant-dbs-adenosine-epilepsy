function B = chanMaskCell(im, xfold)
if nargin < 2
    xfold = 1;
end

M = imopen(imgaussfilt(single(im), 5 * xfold), ...
    strel('disk', round(10 * xfold)));

[~, threshold] = edge(M, 'canny');
B = edge(M, 'canny', threshold);

B = imdilate(B, strel('disk', round(3 * xfold)));

B = imfill(B, 'holes') & ~B;

B = imerode(B, strel('diamond', 1));
B = imerode(B, strel('diamond', 1));
end