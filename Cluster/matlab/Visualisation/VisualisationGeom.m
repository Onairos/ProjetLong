%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des figures g�om�triques
close all
clear all

% Chargement des donn�es
Data0=importdata('3blocsSimple.txt');
Data=Data0(2:end, :);

% Taille des donn�es
disp('nbre de points')
[n1,p]=size(Data)

% Nombre de classe � charger
disp('nombre de classes')
nbcluster=importdata('nbclusters');



if p==2
   disp('Figure 2D')
   cm=colormap(hsv(10*(nbcluster+1)));
   % Visualisation des donn�es brutes
    subplot(1,2,1)  
    hold on
    plot(Data(:,1),Data(:,2),'k.','MarkerSize',20);
     title(['Original data ' ]);
   % R�cup�ration et Visualisation des clusters 
    subplot(1,2,2)
    hold on
    for i=1:nbcluster
        ii=importdata(strcat('cluster.final.',num2str(i)));
        plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',20);
    end
     title(['Clusters number = ' int2str(nbcluster)]);
    
else
    disp('Figure 3D')
    cm=colormap(hsv(10*(nbcluster+1)));
    % Visualisation des donn�es brutes
   
    
    subplot(1,2,1)    
    hold on 
    plot3(Data(:,1),Data(:,2),Data(:,3),'g.','MarkerSize',15);
     title(['Original data ' ]);
   % R�cup�ration et Visualisation des clusters 
    subplot(1,2,2)
        hold on 
    for i=1:nbcluster
        ii=importdata(strcat('cluster.final.',num2str(i)));
        plot3(Data(ii(2:end),1),Data(ii(2:end),2),Data(ii(2:end),3),'.','color',[cm(i*10,:)],'MarkerSize',15);
    end
     title(['Clusters number = ' int2str(nbcluster)]);
    
end

return
