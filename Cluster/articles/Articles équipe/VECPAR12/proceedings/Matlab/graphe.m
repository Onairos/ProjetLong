clear all
close all

%benjamin
%arthus1
arthus2
data=[sparse1(:),sparse3(:),sparse5(:),sparse8(:),sparse_10(:),sparse_12(:),sparse_15(:)];
%data=[sparse3(:),sparse5(:),sparse8(:),sparse_10(:),sparse_12(:),sparse_15(:)];
datamoy=mean(data);
datamin=min(data);
datamax=max(data);

figure
hold on
seuil=[1,3,5,8,10,12,15];
%seuil=[3,5,8,10,12,15];
plot(seuil,datamoy,'k-+','LineWidth',2,'MarkerSize',14)
plot(seuil,datamin,'k.-','LineWidth',2,'MarkerSize',14)
plot(seuil,datamax,'k--','LineWidth',2,'MarkerSize',14)
xlabel('Factor','FontSize',14)
ylabel('nnz(Affinity matrix)','FontSize',14)
%legend('Mean',4)
legend('Mean','Min','Max',4)

ti = get(gca,'TightInset')
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

saveas(1, 'toto.pdf', 'pdf');
