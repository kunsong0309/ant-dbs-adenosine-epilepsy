function [L, A, P] = imMaskVglut2(im)
im = double(im);
thr = mad(im, 1, 'all') * 1.4826 * 3 + median(im, 'all');
% figure; imshow(im);

se = strel('disk', 20);
M = imopen(imgaussfilt((im > thr) * 1, 10), se);
% figure; histogram(M(:)); xline(0.05, ':r');
BWfinal = M > 0.05;
% figure; imshow(M > 0.1125, [0, 1]);

% fudgeFactor = 0.5;
% method = 'canny';
% [~,threshold] = edge(M, method);
% BWs = edge(M, method, threshold * fudgeFactor);
% % figure; imshow(BWs);
% 
% se360 = strel('disk', 1);
% BWsdil = BWs;
% for ii = 1:3
%     BWsdil = imdilate(BWsdil, se360);
% end
% % figure; imshow(BWsdil);
% 
% BWdfill = imfill(BWsdil, 'holes');
% % figure; imshow(BWdfill);
% 
% BWfinal = BWdfill & ~BWsdil;
% % seD = strel('diamond', 1);
% % BWfinal = imerode(BWfinal, seD);
% se360 = strel('disk', 1);
% for ii = 1:3
%     BWfinal = imdilate(BWfinal, se360);
% end
% % figure; imshow(labeloverlay(im, BWfinal, 'Transparency', 0.5));

[B, L] = bwboundaries(BWfinal,'noholes');
% figure;
% imshow(label2rgb(L, @jet, [.5 .5 .5]));
% hold on
% for k = 1:length(B)
%     boundary = B{k};
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.1);
% end

Lmsk = L;
M = im > thr;
mm = mean(M, 'all');
isv = true(size(B));
A = nan(size(B));
P = nan(size(B));
for kk = 1:length(B)
    A(kk) = bwarea(L == kk);
    if A(kk) < 1000
        Lmsk(L == kk) = 0;
        isv(kk) = false;
    end
    P(kk) = 1 - cdf('Poisson', sum(M(L == kk), 'all'), A(kk) * mm); 
    % 1 - cdf('Binomial', sum(M(L == kk), 'all'), A(kk), mm);
    if P(kk) >= 0.05
        Lmsk(L == kk) = 0;
        isv(kk) = false;
    end
end
% figure; imshow(labeloverlay(im, Lmsk, 'Transparency', 0.5));

[~, L] = bwboundaries(Lmsk,'noholes');
A = A(isv);
P = P(isv);
end