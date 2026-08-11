close all; clear
figRes=600; figExt='png'; figName='FigA3';

fig1=figure('Position',[400 400 1000 350], 'color','w');
    posX1=0.05; posY1=0.13;
    width=0.28; heigth=.85;
    gapX=0.05; gapY=0.075;

   for j=1:3  
    for i=1:1
      panel=2*(i-1)+j;
      h=axes(fig1,'Position',[posX1+(width+gapX)*(j-1) posY1-(heigth+gapY)*(i-1) width heigth]);
      eval(['axes' num2str(panel) '=h;']);clear h
    end
   end

fnm{1}='FigA3a_05495500_48.in';
fnm{2}='FigA3b_02349900_48.in';
fnm{3}='FigA3c_04059500_48.in';
 
for pp=1:3
   eval(['axes(axes' num2str(pp) '); hold all']);
   eval(['drawGraph(fnm{' num2str(pp) '},pp);']);
end
print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
  eval(print_str);

fig2=figure('Position',[10 558 440 350], 'color','w');
    posX=0.133; 

for pp=1:3
   axes(fig2,'Position',[0.133 0.13 (0.97-posX) 0.85]); hold all
   eval(['drawGraph(fnm{' num2str(pp) '},pp);']);
   figName1=[figName char(96+pp)];
   print_str=[ 'print(gcf,' char(39) figName1 char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
   eval(print_str); 
   delete(gca);
end  
close(fig2);  
  
function drawGraph(filename,panelNr)
 AUXDATA=load(filename);
 Yo=AUXDATA(:,1);
 Yp=AUXDATA(:,2);

 Yerror=(Yo-Yp);
 nanind=~isnan(Yerror); Yerror=Yerror(nanind);Yp=Yp(nanind);Yo=Yo(nanind);

 SSQ=sum(Yerror.^2);
 n=size(Yo,1);
 NSE= 1 - SSQ/var(Yo)/(n-1);

 r=corrcoef(Yp,Yo); r=r(2,1);
 alpha=std(Yp)/std(Yo);
 beta=(1+(mean(Yp)-mean(Yo))/std(Yo)); 
 KGEm=1-sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2);

 x=(1:1:n)';
 fntsz=14;
 fntwght='bold';
 set(gca,'Fontsize',fntsz,'Layer','top');

 plot(x,Yp,'-','markersize',4,'linewidth',1.5,'Color',0.4*ones(1,3));
 plot(x,Yo,'b.','markerfacecolor','none','markersize',8,'linewidth',1);
 box on

 xlabel('Position within the series','Fontsize',fntsz,'fontweight',fntwght);
 ylabel('Runoff (mm d^{-1})','Fontsize',fntsz,'fontweight',fntwght);
 ht=gca; ht.YRuler.TickLabelFormat='%3.0f'; 

 xlim([0 size(AUXDATA,1)+5]);
 yT=get(gca,'YTick');
 set(gca,'YTick',yT,'Ylim',[0 yT(end)+(yT(end)-yT(end-1))/2]);

 frmt='%1.3f'; j=1;
 txt{j}=['\bf{NSE}\rm= ' sprintf(frmt,NSE)]; j=j+1;
 txt{j}=['\bf{KGE*}\rm= ' sprintf(frmt,KGEm)]; j=j+1;
 txt{j}=['\bf{\rho/\alpha}\rm= ' sprintf(frmt,r/alpha)]; j=j+1;
 txt{j}=['\bf{\rho}\rm= ' sprintf(frmt,r)]; j=j+1;
 txt{j}=['\bf{\alpha}\rm= ' sprintf(frmt,alpha)]; j=j+1;
 clear j

 AddText(txt,fntsz-2,0.05,0.85,'left',fntwght);
 AddText([char(96+panelNr) ')'],fntsz+2,-0.15,1,'left',fntwght);
end %----

function hTxt=AddText(txt_label,fntsz,xv,yv,align,fntwght)
  lim_x=xlim; lim_y=ylim; 
  xpos=(lim_x(2)-lim_x(1))*xv + lim_x(1);
  ypos=(lim_y(2)-lim_y(1))*yv + lim_y(1);
  hTxt=text(xpos,ypos,txt_label,'Fontsize',fntsz,'HorizontalAlignment',align,'Fontweight',fntwght);
end %----
