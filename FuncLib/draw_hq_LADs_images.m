% create figure for LADs
% 
%myDir = 'LocsLib/Stiffness'; % gets directory
%myFiles = dir(fullfile(myDir,'*.txt')); %gets all wav files in struct
cell_index = 8;
% lads_new = data{cell_index,1}.lads_new;
% lads_flat = data{cell_index,1}.lads_flat;
name = extractBefore(data{cell_index,1}.name,'.');


%% create figure for LADs
lads = data{cell_index,1}.lads;
monitor2 = figure();
set(gcf,'name','Storm Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.5 0.7]);
subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);    
subplot(1,1,1)
scatter(lads(:,1),lads(:,2),3,'MarkerEdgeColor','y','MarkerFaceColor','y'); hold on

% pixels_um_x = 3*1000;
% pixels_um_y = 0.3*1000;
% rectangle('Position',[min(x)-900 min(y)-900 pixels_um_x pixels_um_y],'facecolor','w'); hold on
axis equal
x = lads(:,1);
y = lads(:,2);
% xlim([min(x)-1000 max(x)+1000])
% ylim([min(y)-1000 max(y)+1000])
L=20000;
xlim([min(x)-1000 min(x)+L+1000])
ylim([min(y)-1000 min(y)+L+1000])
set(gcf, 'InvertHardCopy', 'off');
set(gca,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
print(gcf,['LADsLib/',['LADs_',name],'.png'],'-dpng','-r800');
% pause(2)
close(monitor2)


%%
% figure('Position',  [300, 300, 1400, 100])
% set(gcf,'NumberTitle','off','color','k','units','normalized');
% x = lads_flat(:,1);
% y = lads_flat(:,2);
% yy = smoothdata(y,'gaussian',15);
% bar(x,yy,1,'FaceColor',[0.8 0.5 0.5]); hold on
% scatter(lads_new(:,1),lads_new(:,2),1,'.y'); hold on
% 
% axis equal
% xlabel('Circumference / nm')
% ylabel('Height / nm')
% xlim([0,max(x)])
% ylim([0,1200])
% legend('off')
% set(gcf, 'InvertHardCopy', 'off');
% set(gca,'color','k','box','on','BoxStyle','full','XColor','w','YColor','w');
%%
%lmt = data{cell_index,1}.lmt;
%affinity = data{cell_index,1}.affinity;
%figure('Position',  [300, 300, 1400, 100])
%set(gcf,'NumberTitle','off','color','k','units','normalized');
%for i = 1:length(affinity)
%    area(lmt(i,:),[affinity(i),affinity(i)],'FaceColor',[0.8 0.5 0.5],'EdgeColor',[0.8 0.5 0.5]); hold on
%end

%xlabel('Circumference / nm')
%ylabel('Affinity')

%xlim([0,max(x)])
%ylim([0,1.2])
%set(gcf, 'InvertHardCopy', 'off');
%set(gca,'color','k','box','on','BoxStyle','full','XColor','w','YColor','w');

%print(gcf,['LADsLib/',['Affinity_',name],'.png'],'-dpng','-r800');