% Display Error
close all;clear all;
cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp');
load('Phase129810');

ColorVector=hsv(10);
Impacts=[28.5638528937160,26.0770511040048,23.5211360663909,20.9028818439732,18.2292277203549,15.5072598081303,12.7441922682240,9.94734818985714,7.12414018181586,4.28205072646230,1.42861234855574,-1.42861234855574,-4.28205072646230,-7.12414018181586,-12.7441922682241,-15.5072598081303,-18.2292277203549,-20.9028818439732,-23.5211360663909];
Impacts=int64(Impacts)';
S = get(0, 'ScreenSize');
fntsz=13;
h1 = figure('Visible','on','Name','Data Comparison','Position',[.25*S(3) 35 .6*S(3) .5*S(3)],...
        'Color', [1 1 1]);
ax(1) = axes('Parent', h1, 'Position', [.1 .1 .8 .35], 'FontSize', fntsz);
ax(2) = axes('Parent', h1, 'Position', [.1 .56 .8 .35], 'FontSize', fntsz);
xlabel(ax(1),'Impact Parameter');
%xlabel(ax(2),'Impact Parameter');
ylabel(ax(1),'Error, [km/s]');
ylabel(ax(2),'Error, [km/s]');

title(ax(2),'Average Total Uncertainty, 14.5 kHz');
title(ax(1),'Average Total Uncertainty, 53.5 and 68.5 kHz');

hold(ax(1),'on');
hold(ax(2),'on');
cd('T:\IDS\Data Repository\TEMP');
quad53Pos=load('quad53Pos');
quad53Pos=quad53Pos.TEMP;
quad53Neg=load('quad58Neg');
quad53Neg=quad53Neg.TEMP;
quad68=load('quad68');
quad68=quad68.temp;
quad14=load('quad14');
quad14=quad14.quad;



L1=[18.510,199.88,54.385,234.16,91.124,269.11,127.42,305.29,163.73,339.05];

% legend_str=[];
% for i=1:10
%     
%     plot(ax(1),linspace(-23,28,19),temp(end:-1:1),'color',ColorVector(i,:));
%     legend_str = [legend_str; {num2str(i)}];
% end
% columnlegend(3, legend_str, 'Location', 'NorthWest');
% set(ax(1),'XTick',Impacts(end:-1:1));
% set(ax(2),'XTick',Impacts(end:-1:1));
%set(ax(1),'XLim',[-23 28]);


Phase=(PhaseVelocity.Phase*180);
Phase1=zeros(5,2);
Phase1(:,1)=Phase(1:5);
Phase1(:,2)=Phase(6:10);
%legend(ax(1),L1);
legend_str = []; 
     for i=1:10,  
          temp=quad14(:,i);
          plot(ax(2),linspace(-23,28,19),temp(end:-1:1),'color',ColorVector(i,:),'LineWidth',2); hold on; 
          legend_str = [legend_str; {num2str(L1(i))}];
     end
     %columnlegend(5, legend_str, 'Location', 'Best','FontSize',3,'box','on');

     
plot(ax(1),linspace(-23,28,19),mean(quad53Pos(end:-1:1,:),2),'color','blue','LineWidth',2);
plot(ax(1),linspace(-23,28,19),mean(quad53Neg(end:-1:1,:),2),'color','red','LineWidth',2);
plot(ax(1),linspace(-23,28,19),mean(quad68(end:-1:1),2),'color','green','LineWidth',2);

legend(ax(1),'53 kHz, Positive','53 kHz, Negative','68 kHz, Negative','Location','Best');
set(ax(2),'XLim',[-25 28]);
set(ax(1),'XLim',[-25 28]);
set(ax(2),'YLim',[0 7]);
grid(ax(2),'on');
grid(ax(1),'on');
box(ax(1),'on');
box(ax(2),'on');