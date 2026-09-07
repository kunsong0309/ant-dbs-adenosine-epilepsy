function [N, L] = imSegmentDots(im, thr)
[B, L] = bwboundaries(single(im) > thr,'noholes');
N = numel(B);
end