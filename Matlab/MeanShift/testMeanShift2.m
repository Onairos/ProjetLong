clear all;
close all;

OImg=imread('Oimg2.jpg');

cform = makecform('srgb2lab');
lab_he = applycform(OImg,cform);

ab=double(lab_he(:,:,2:3));
subplot 211
imshow(OImg);
subplot 212
imshow(lab_he);

nrows=size(ab,1);
ncols=size(ab,2);

ab=reshape(ab,nrows*ncols,2);
nColors = 3;
%repeat the clustering 3 times to avoid
[clusterCent, point2cluster, clustMembsCell]  = MeanShiftCluster(ab',30);

pixel_labels = reshape(point2cluster, nrows, ncols);

numClust=length(clustMembsCell);

figure
imshow(pixel_labels, []);

figure

segmented_images = cell(1,numClust);
rgb_label = repmat(pixel_labels, [1 1 3]);
for k=1:nColors
    color=OImg;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end
imshow(segmented_images{1}), title('objects in cluster 1')
figure
imshow(segmented_images{2}), title('objects in cluster 2')
figure
imshow(segmented_images{3}), title('objects in cluster 3')
