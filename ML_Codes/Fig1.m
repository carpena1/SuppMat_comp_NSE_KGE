close all; clear
% Figure name, resolution and extension, marker size, font size
figName='Fig1'; figRes=300; figExt='png'; 

fig1=figure; set(fig1,'Position',[400 400 1000 350], 'color','w');
    posX1=0.05;
    posY1=0.13;
    width=0.28;
    heigth=.85;
    gapX=0.05;
    gapY=0.075;

   for j=1:3  
        for i=1:1
          panel=2*(i-1)+j;
          h=axes(fig1,'Position',[posX1+(width+gapX)*(j-1) posY1-(heigth+gapY)*(i-1) width heigth]);
          eval(['axes' num2str(panel) '=h;']);clear h
        end
   end

clear fnm

 fnm{1}='Fig1a_Example.txt';
 fnm{2}='Fig1b_Example.txt';

yLbl{1}='Target variable';
yLbl{2}=yLbl{1};
yLbl{3}='Standardized target variable';



for pp=1:2
   eval(['axes(axes' num2str(pp) '); hold all']);
   eval(['drawGraph(fnm{' num2str(pp) '},pp,[0 8],yLbl{' num2str(pp) '});']);
end
axes(axes3);
AUXDATA=load(fnm{2});
Yo=AUXDATA(:,1); SDo=nanstd(Yo); AVo=nanmean(Yo);
Yp=AUXDATA(:,2); SDp=nanstd(Yp); AVp=nanmean(Yp);
Yo=(Yo - AVo)/SDo;
Yp=(Yp - AVp)/SDp;

writeArray([Yo Yp],20,'Example_b_std.in');
drawGraph('Example_b_std.in',3,[-3 4],yLbl{3});
delete Example_b_std.in
print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
eval(print_str); 
   
function drawGraph(filename,panelNr,yL,yLblTxt)

AUXDATA=load(filename);
Yo=AUXDATA(:,1);
Yp=AUXDATA(:,2);

Yerror=(Yo-Yp);
nanind=~isnan(Yerror); Yerror=Yerror(nanind);Yp=Yp(nanind);Yo=Yo(nanind);

SSQ=sum(Yerror.^2);
n=size(Yo,1);
RMSE=(SSQ/n)^.5;
SD=std(Yo);
NRMSE=RMSE/SD*100;
NSE= 1 - SSQ/var(Yo)/(n-1);

r=corrcoef(Yp,Yo); r=r(2,1);
alpha=std(Yp)/std(Yo);
beta=mean(Yp)/mean(Yo); 
KGE=1-sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2);


x=(1:1:n)';
lnwdth=1.1;
fntsz=14;
fntwght='bold';
serCol=[1 1 1]*.8;
plot([0 50],[0 0],'-','color',[1 1 1],'linewidth',2);

 y1=(Yp+RMSE)'; y2=(Yp-RMSE)'; 
 h=area(x',[y1;y2-y1]');  hold all
   set(h,'LineStyle','none','HandleVisibility','off');
   set(h(1),'FaceColor','none'); 
   set(h(2),'FaceColor',serCol); 
   uistack(h,'bottom');
   set(gca,'Layer','top');

p(1)=plot(x,Yo,'bo','markerfacecolor','none','markersize',6,'linewidth',1); hold all
p(2)=plot(x,Yp,'-','markersize',4,'linewidth',1.5,'Color',[0.494 0.184 0.556]);
plot(x,mean(Yo)*x./x,'k--','markerfacecolor','b','markersize',4,'linewidth',1);

x0=[x(1) x(end)];
p(3)=plot(x0,[1 1]*(mean(Yo)+SD),'k--','linewidth',lnwdth-.3);
plot(x0,[1 1]*(mean(Yo)-SD),'k--','linewidth',lnwdth-.3);

xlabel('Position within the series','Fontsize',fntsz,'fontweight',fntwght);
ylabel(yLblTxt,'Fontsize',fntsz,'fontweight',fntwght);
ht=gca; ht.YRuler.TickLabelFormat='%3.0f'; 

LegLoc='NE'; LegOri='vertical'; LegBox='on';

if panelNr==1
legTxt={'\it o_i \rm','\it p_i \rm',['\it o \rm' char(177) '\its_o \rm'],['\it o \rm' char(177) 'RMSE_{} ']};
leg1=legend([p(1) p(2) p(3) h(2)],legTxt,...
    'Location',LegLoc,'box',LegBox,'Fontsize',fntsz-3,'Orientation',LegOri,'color','none');
  legPos=get(leg1,'Position'); legPos(1)=legPos(1)+0.01; legPos(2)=legPos(2)+0.01;
  set(leg1,'Position',legPos,'box','off');
  leg1.ItemTokenSize = [22,12];
  LinePosX=40.1;LinePosY=6.72; 
  plot([LinePosX LinePosX+1],LinePosY*[1 1],'k-','linewidth',0.8); % Add line to denote mean O
  plot([LinePosX LinePosX+1],(LinePosY-.5)*[1 1],'k-','linewidth',0.8); % Add line to denote mean O
end
  ylim(yL);
    RMSE=rounddec(RMSE,2);
  NSE=rounddec(NSE,2);
  KGE=rounddec(KGE,2);

frmt='%1.2f'; j=1;
txt{j}=['\bf{RMSE}\rm= ' sprintf(frmt,RMSE)]; j=j+1;
txt{j}=['\bf{NSE}\rm= ' sprintf(frmt,NSE)]; j=j+1;
txt{j}=['\bf{KGE}\rm= ' sprintf(frmt,KGE)]; j=j+1;
clear j

AddText(txt,fntsz-2,0.05,0.9,'left',fntwght);
AddText([char(96+panelNr) ')'],fntsz+2,-0.15,1,'left',fntwght);
set(gca,'Fontsize',fntsz);
box on

end %----

function hTxt=AddText(txt_label,fntsz,xv,yv,align,fntwght)
  lim_x=xlim; lim_y=ylim; 
  xpos=(lim_x(2)-lim_x(1))*xv + lim_x(1);
  ypos=(lim_y(2)-lim_y(1))*yv + lim_y(1);
  hTxt=text(xpos,ypos,txt_label,'Fontsize',fntsz,'HorizontalAlignment',align,'Fontweight',fntwght);
end %----
