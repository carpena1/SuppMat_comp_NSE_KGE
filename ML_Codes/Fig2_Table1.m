close all; clear 
global  Mksz fntSz figExt figRes figName

clear DATA AUXDATA CAMELS

% Figure resolution and extension, marker size, font size
figRes=300; figExt='png'; Mksz=6; fntSz=16; 
   
colNSE=4; colKGE=5; colN=12; colRBias=13; colSetGroupNr=15;

% Datasets 1 - 9 
   fnm='NSEvsKGE_datasets_0-9.txt'; fprintf('Reading %s ... ',fnm);
   A=importdata(fnm); DATA=A.data; clear A fnm
   fprintf(1,'%i rows and %i columns data available\n',size(DATA)); 

 nseth=0;
 auxind=DATA(:,colNSE)>nseth;
 disp(['Nr. of NSE values <= 0 discarded: ' num2str(length(auxind)-sum(auxind))]);
 DATA=DATA(auxind,:); clear auxind
 fprintf(1,'Number of total cases after discarding NSE<=%1.1f: %i\n',nseth,sum(~isnan(DATA(:,colNSE))));

 AUXDATA=DATA; clear DATA
 fprintf(1,'Total number of pair values: %i\n',sum(AUXDATA(:,colN)));

biasth=5; 
ksth=0.05;
nseth=0.2; 

 DiscardedId=find(abs(AUXDATA(:,colRBias))>biasth); 
 disp(['Cases with bias > ' num2str(biasth) '%: ' num2str(length(DiscardedId))]);
 disp(' ');
 
%%%% CAMELS DATA SET %%%%%%% 
 fprintf('Reading CAMELS dataset ... ');
 fnm='NSEvsKGE_datasets_camels_48.txt'; 
   A=importdata(fnm); CDATA=A.data; clear A fnm
   fprintf(1,'%i rows and %i columns data available\n',size(CDATA)); 
 
    auxind=CDATA(:,colNSE)>nseth;
    disp(['Nr. of NSE values <= ' num2str(nseth) ' discarded: ' num2str(length(auxind)-sum(auxind))]);
    CAMELS=CDATA(auxind,:); clear auxind CDATA
    fprintf(1,'Number of total cases after discarding NSE<%1.1f: %i\n',nseth,sum(~isnan(CAMELS(:,colNSE))));
 
    fprintf(1,'Total number of pair values: %i\n',sum(CAMELS(:,colN)));

 CAMELS_BiasedId=find(abs(CAMELS(:,colRBias))>biasth); 
 disp(['Cases with bias > ' num2str(biasth) '%: ' num2str(length(CAMELS_BiasedId))]);
 disp(' ');
 
%%% -----------------------------------------
   IntX=[.2 0.95]; IntY=IntX;
   axPos=[0.137 0.135 0.835 0.85];
   figPos=[10 558 440 390];
   MaxGroupNr=9;
   legTxt=cell(MaxGroupNr,1);
   outsymb='x';
%%% -----------------------------------------

 %%%% Fig 2 - Boxplots %%%%
 VarCol=[colNSE colKGE]; VarYlim=[0 1.05; -1.1 1.1]; VarName={'NSE','KGE'};
 figure('Position',[10 558 570 450],'Color','w');
 Hght=.449; xp1=.105; xp2=.85; yp1=.545; yp2=.1; Wd1=.75; Wd2=.13;
 pos(1,:)=[xp1 yp1 Wd1 Hght];
 pos(2,:)=[xp1 yp2 Wd1 Hght];
 pos(3,:)=[xp2 yp1 Wd2 Hght];
 pos(4,:)=[xp2 yp2 Wd2 Hght]; clear xp* yp* Wd* Hght

   axes1=axes('Position',pos(1,:));
   axes2=axes('Position',pos(2,:));
   axes3=axes('Position',pos(3,:));
   axes4=axes('Position',pos(4,:));
   
 for i=1:2
  for ss=1:MaxGroupNr 
    legTxt{ss}=sprintf('Set %i',ss);
  end; clear ss
  
  SetIndex=AUXDATA(:,colSetGroupNr); % Data set group flag

  eval(['axes(axes' num2str(i) ' )']);
   boxplot(AUXDATA(:,VarCol(:,i)),SetIndex,'Symbol',outsymb); hold all
    h=findobj(gca,'LineStyle','--'); set(h, 'LineStyle','-');
    for ss=1:MaxGroupNr
     auxid=find(SetIndex==ss);
     fprintf(1,'Set %i: %s [%1.3f - %1.3f]\n',ss,VarName{i},min(AUXDATA(auxid,VarCol(:,i))),max(AUXDATA(auxid,VarCol(:,i))));
      plot(ss,nanmean(AUXDATA(auxid,VarCol(:,i))),'Color','b','Marker','o','Markersize',4);
      n(ss)=length(auxid);
    end; clear auxid
    disp(' ');
    set(gca,'xticklabel',legTxt);

    formatFig(gcf,pos(i,:),IntX,IntY,'Data set',VarName{i}); 
    xlim([.5 MaxGroupNr+.5]); ylim(VarYlim(i,:));
    set(gca,'Xtick',1:1:MaxGroupNr);
    set(gca,'Ytick',floor(VarYlim(i,1)):i/10:floor(VarYlim(i,2)));
    yTks=get(gca,'Ytick');
    
   ss=10; 
   fprintf(1,'Set %i: %s [%1.3f - %1.3f]\n\n',ss,VarName{i},min(CAMELS(:,VarCol(:,i))),max(CAMELS(:,VarCol(:,i))));

    
  eval(['axes(axes' num2str(i+2) ' )']);
   nL=length(SetIndex); nC=size(CAMELS,1);
   boxplot([AUXDATA(:,VarCol(:,i)); CAMELS(:,VarCol(:,i))],[ones(nL,1);2*ones(nC,1)],'Symbol',outsymb); hold all
     h=findobj(gca,'LineStyle','--'); set(h, 'LineStyle','-');
   plot([1 2],[nanmean(AUXDATA(:,VarCol(:,i))) nanmean(CAMELS(:,VarCol(:,i)))],...
       'bo','Markersize',4);

    ylim(VarYlim(i,:)); xp=xlim; xlim([xp(1)-.5 xp(2)]);
    set(gca,'Ytick',floor(VarYlim(i,1)):i/10:floor(VarYlim(i,2)));
    set(gca,'Fontsize',fntSz-2,'Ytick',yTks,'xticklabel',{'1 - 9','10'},'Yticklabels','');
    
   eval(['axes(axes' num2str(i) ' )']);
   AddText(gcf,[char(96+i) ')'],fntSz+2,-.13,.98,'left','bold');
 end; clear i
set(axes1,'xticklabel',''); 
yTcksLbs=get(axes1,'Yticklabel');
yTcksLbs{1}=strrep(yTcksLbs{1},'0.0',''); set(axes1,'Yticklabel',yTcksLbs);

 figName='Fig2_NSE_KGE_bxplts_sets';
 print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str); figName='';
    disp(' ');

    
function hTxt=AddText(figh,txt_label,fntsz,xv,yv,align,varargin)
cl='k'; fntwgt='normal';

switch nargin
    case 7, if length(varargin{1})==1, cl=varargin{1}; else fntwgt=varargin{1}; end
    case 8, fntwgt=varargin{1}; cl=varargin{2};
end

try figure(figh);catch, subplot(figh); end
  lim_x=xlim; lim_y=ylim; 
  xpos=(lim_x(2)-lim_x(1))*xv + lim_x(1);
  ypos=(lim_y(2)-lim_y(1))*yv + lim_y(1);
  hTxt=text(xpos,ypos,txt_label,'Fontsize',fntsz,'Fontweight',fntwgt,'HorizontalAlignment',align,'Color',cl);
end %----


function formatFig(hFig,axPos,rngX,rngY,xLbTxt,Eq6LbTxt,varargin)
 legendTxt=''; letter='';
 switch nargin
    case 7, letter=varargin{1}; 
    case 8, letter=varargin{1}; legendTxt=varargin{2}; LOC='NW';
    case 9, letter=varargin{1}; legendTxt=varargin{2}; LOC=varargin{3}; 
 end

 global fntSz figExt figName figRes
 figure(hFig);
  set(gca,'Fontsize',fntSz-1,'Position',axPos);
  %axis equal
  xlim(rngX);  ylim(rngY); 
  set(gca,'Xtick',rngX(1):0.1:rngX(2),'Ytick',rngY(1):0.1:rngY(2));
  hA=gca; hA.YRuler.TickLabelFormat='%5.1f'; % Y-axis in %
  hA.XRuler.TickLabelFormat='%1.1f'; % X-axis in %
  if ~isempty(legendTxt); legend(legendTxt,'Location',LOC,'box','off'); end
  ylabel(Eq6LbTxt,'Fontsize',fntSz,'fontweight','bold');
  xlabel(xLbTxt,'Fontsize',fntSz,'fontweight','bold');
  AddText(gcf,letter,fntSz+2,-.15,1,'left','bold');
  if ~isempty(figName)
   print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str); figName='';
 end
end%--

