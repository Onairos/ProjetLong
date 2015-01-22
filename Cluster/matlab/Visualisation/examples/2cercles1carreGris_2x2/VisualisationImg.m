%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des images
close all
clear all

% Chargement des données
%fichierImag=importdata('nomfichier');
%I1=imread(fichier);


load 2Cercles1CarreGris.mat

I=double(I1);

disp('nbre de points')
[n1,m1,p]=size(I)

% Nombre de classe à charger
disp('nombre de classes')
nbcluster=importdata('nbclusters');



if p==1
    disp('image niveau de gris')
   % Récupération des clusters et définition du niveau de gris moyen par
   % classe
    IDX=zeros(n1*m1,1); 
    I1=I';
    Imag=I1(:);
    for i=1:nbcluster
        A=importdata(strcat('cluster.final.',num2str(i)));
        LevelColor=mean(Imag(A(2:end,1),:));
        IDX(A(2:end,1))=repmat(LevelColor,size(A,1)-1,1);
    end

  
     % Visualisation de la comparaison entre image brute et résultat de
     % classification
    figure
    hold on
    subplot(1,2,1)
    I=uint8(I);
    imshow(I)
    title(['Original data']);
    
    subplot(1,2,2)
    I4=reshape(IDX,m1,n1);
    I4=ipermute(I4,[2,1,3]);
    I5=uint8(I4);
    imshow(I5)
    rotate3d on
    title(['Clusters number = ' int2str(nbcluster)]);
     
   
else
    
    disp('image couleur')
   % Récupération des clusters et définition du niveau de gris moyen par
   % classe
    IDX=zeros(n1*m1,3);
    Imag=zeros(n1*m1,3);
    G1=I(:,:,1)';
    G2=I(:,:,2)';
    G3=I(:,:,3)';
    Imag(:,1)=G1(:);
    Imag(:,2)=G2(:);
    Imag(:,3)=G3(:);
    for i=1:nbcluster
        A=importdata(strcat('cluster.final.',num2str(i)));
        LevelColor=mean(Imag(A(2:end,1),:));
        IDX(A(2:end,1),:)=repmat(floor(LevelColor),size(A,1)-1,1);
    end
    
    % Visualisation de la comparaison entre image brute et résultat de
    % classification
    figure
    hold on
    subplot(1,2,1)
    I=uint8(I);
    imshow(I)
    title(['Original data ' ]);
    
    subplot(1,2,2)
    I4=reshape(IDX,m1,n1,3);
    I4=ipermute(I4,[2,1,3]);
    I5=uint8(I4);
    imshow(I5)
    rotate3d on
    title(['Clusters number = ' int2str(nbcluster)]);
     


end



return

