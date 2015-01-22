
% Génération des datas pour SC parallèle


clc
clear all
close all



I=imread('cameraman.tif'); % exemple d'image à charger

fid=fopen('Cameraman.txt','w'); % fichier data à créer;


n=size(I,1);m=size(I,2);
fprintf(fid, '%i %i \n', 2,1); % matrice 2D = en niveau de gris
fprintf(fid, '%i %i \n', n,m); % dimensions de la matrice

for i=1:n
    for j=1:m
        fprintf(fid, '%12.8f\n',I(i,j));
    end
end
fclose(fid);
disp('fini');