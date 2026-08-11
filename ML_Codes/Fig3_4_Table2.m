close all; clear 
global  Mksz fntSz ser_color ser_symbol figExt figRes figName
ser_color=['b' 'm' 'r' 'k' 'g' 'r' 'b' 'k' 'm' 'g' 'b' 'm']; % Define color for the series to plot in each panel
ser_symbol=['o' '^' 'd' 'v' '<' '>' 'p' 'h' 's' '+' '*' 'x']; % Define symbol for the series to plot in each panel

% Figure resolution and extension, marker size, font size
   figRes=300; figExt='png';  Mksz=6; fntSz=16; 

   fnm='NSEvsKGE_datasets_0-9.txt'; fprintf('Reading %s ... ',fnm);
   A=importdata(fnm); DATA=A.data; clear A fnm
   fprintf(1,'%i rows and %i columns data available\n',size(DATA)); 

   
j=1;
colSet=j; j=j+1;
colSerNr=j; j=j+1;
colNormzd=j; j=j+1;
colNSE=j; j=j+1; 
colKGE=j; j=j+1;
colRMSE=j; j=j+1; 
colmYo=j; j=j+1;colmYp=j; j=j+1;
colsdYo=j; j=j+1; colsdYp=j; j=j+1;
col_r=j; j=j+1;
colN=j; j=j+1;
colRBias=j; j=j+1;
colKGEm=j; j=j+1;
colSetGroupNr=j;

% Datasets 1 - 9 
 MaxGroupNr=9;
 DATA=DATA(DATA(:,colSetGroupNr)<=MaxGroupNr,:); 
 fprintf(1,'Number of total cases considered for the %i sets selected: %i\n',...
     MaxGroupNr,sum(~isnan(DATA(:,colNSE))));

 
 auxind=DATA(:,colNSE)>0;
 disp(['Nr. of NSE values <= 0 discarded: ' num2str(length(auxind)-sum(auxind))]);
 DATA=DATA(auxind,:); clear auxind
 fprintf(1,'Number of total cases after discarding NSE<0: %i\n',sum(~isnan(DATA(:,colNSE))));

 AUXDATA=DATA; clear DATA
 
 AUXDATA(:,colNSE)=rounddec(AUXDATA(:,colNSE),5); 
 AUXDATA(:,colKGEm)=rounddec(AUXDATA(:,colKGEm),3); 

r=rounddec(AUXDATA(:,col_r),5); % Pearson correlation coefficient
alpha=AUXDATA(:,colsdYp)./AUXDATA(:,colsdYo); % Variability component
beta=1+((AUXDATA(:,colmYp)-AUXDATA(:,colmYo))./AUXDATA(:,colsdYo));

NSE_prd=rounddec((r-1).^2+2.*alpha.*(r-1)+1 - (1-AUXDATA(:,colKGEm)).^2,3);
KGEm_prd=rounddec(1-sqrt(max(0,(r-1).^2+2*alpha.*(r-1)+1-AUXDATA(:,colNSE))),3);
 AUXDATA(:,end+1)=KGEm_prd; colKGEm_prd=size(AUXDATA,2); % KGE* from NSE
 AUXDATA(:,end+1)=NSE_prd; colNSE_prd=size(AUXDATA,2); % NSE from KGE*  

Rho2Alpha=r./alpha; 
%%% ---
 biasth=5;
 BiasedId=find(abs(AUXDATA(:,colRBias))>biasth); 
 disp(['Cases with bias > ' num2str(biasth) '%: ' num2str(length(BiasedId))]);
 DiscardedId=find(beta<0.95); 
 disp(['Cases with beta* < 0.95: ' num2str(length(DiscardedId))]);
  
%%% -----------------------------------------
   IntX=[0 1]; IntY=[-1 1];
   axPos=[0.137 0.135 0.835 0.85];
   figPos=[10 558 440 390];
   legTxt=cell(MaxGroupNr,1);

   Eq2LbTxt='KGE (Eq. 2)';
   Eq6LbTxt='KGE* (Eq. 6)';
   Eq4LbTxt='KGE* (Eq. 4)';
   Eq7LbTxt='KGE* (Eq. 7)';
   func_L='1-sqrt(1.911*(x)-3.909*sqrt(x)+2)';
   func_U='1-sqrt(2.114*(x)-4.111*sqrt(x)+2)';

  % Discard cases if outside uncertaintybound 
   clear DiscardedId
   KGEm=AUXDATA(:,colKGEm);
   x=AUXDATA(:,colNSE); eval(['yLB= ' func_L '; yUB= ' func_U ';']);
   DiscardedId=and(KGEm<0.95,or(KGEm>yUB,KGEm<yLB));
   clear x yLB yUB KGEm
   disp(['Cases (KGE*<0.95) lying outside de Uncertainty Band: ' num2str(sum(DiscardedId))]);
   disp(' ');
   
   UBmaxVal=1;
   x=(0.15:0.01:UBmaxVal); eval(['yLB= ' func_L '; yUB= ' func_U ';']); 
   
%%% -----------------------------------------
 SetIndex=AUXDATA(:,colSetGroupNr); % Data set group flag


 %%%%  Fig 3a  %%%%
   disp('KGE vs Eq.6');
   figure('Position',figPos,'Color','w');
      for ss=1:MaxGroupNr
       selected=find(SetIndex==ss);
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
       KGE=rounddec(AUXDATA(selected,colKGE),3); KGE(KGE<IntY(1))=IntY(1); 
       plot(AUXDATA(selected,colKGEm_prd),KGE,'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
     legTxt{ss}=sprintf('Set %i',ss);
     end; clear ss y KGE
     plot(IntY,IntY,'k--');  
     figName='Fig3a_KGE_vs_Eq6';
     formatFig(gcf,axPos,IntX,IntY,Eq6LbTxt,Eq2LbTxt,'a)',legTxt,'SW'); 
     
 
 %%%%  Fig 3b  %%%%
  IntX=[0.2 1]; IntY=IntX;
  disp('KGE* vs Eq.6');
   figPos(2)=50; 
   figure('Position',figPos,'Color','w');
     for ss=1:MaxGroupNr
       selected=find(SetIndex==ss);     
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
      plot(AUXDATA(selected,colKGEm_prd),AUXDATA(selected,colKGEm),'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      legTxt{ss}=sprintf('Set %i',ss);
     end; clear ss y
     plot(IntX,IntY,'k--'); hold all
     figName='Fig3b_KGEm_vs_Eq6';
     formatFig(gcf,axPos,IntX,IntY,Eq6LbTxt,Eq4LbTxt,'b)',legTxt);
  
   
   AUXDATA(:,colNSE)=rounddec(AUXDATA(:,colNSE),3); 
     
   %%%%  Fig 4a  %%%%   
  disp('KGE vs NSE');
   figPos(1)=500; figPos(2)=558;  
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);
    for ss=1:MaxGroupNr
      selected=find(SetIndex==ss);
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
      KGE=rounddec(AUXDATA(selected,colKGE),3); KGE(KGE<IntY(1))=IntY(1); 
      plot(AUXDATA(selected,colNSE),KGE,'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      legTxt{ss}=sprintf('Set %i',ss);
    end; clear ss KGE

     hBand=area(x,[yUB;yLB-yUB]'); 
       set(hBand,'LineStyle','none','HandleVisibility','off');
       set(hBand(1),'FaceColor','none'); 
       set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
       uistack(hBand,'bottom');
    

    figName='Fig4a_KGE_vs_NSE';
    formatFig(gcf,axPos,IntX,IntY,'NSE^{ }',Eq2LbTxt,'a)',legTxt);
 
    feq='1 - sqrt(2).*(1-sqrt(x))'; fitTxt{1}=['1 - ' char(1140) '2(1 - ' char(1140) 'NSE)'];
    feq2='0.5*(sqrt(2)- (1-(x)) ).^2';   

  %%%%  Fig 4b  %%%%   
  disp('KGE* vs NSE');
   figPos(1)=500; figPos(2)=50;  
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);
    for ss=1:MaxGroupNr
      selected=find(SetIndex==ss);
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
      plot(AUXDATA(selected,colNSE),AUXDATA(selected,colKGEm),'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      legTxt{ss}=sprintf('Set %i',ss);
      legTxt{ss+1}='Eq. 7';

        NrDiscardedSet(ss)=sum(DiscardedId(selected)==1);
      fprintf(1,' Cases lying outside UB in Set %i: %i\n',ss,NrDiscardedSet(ss));


    end; clear ss
    eval(['fplot(@(x) ' feq ',[0.2 1],' char(39) 'k-' char(39) ');']);
 
   SepDiscarded=[AUXDATA(DiscardedId,colNSE),AUXDATA(DiscardedId,colKGEm) Rho2Alpha(DiscardedId)];
   Above105=SepDiscarded(SepDiscarded(:,end)>1.05,1:3);
   Below95=SepDiscarded(SepDiscarded(:,end)<.95,1:3); 
   plot(Above105(:,1),Above105(:,2),'ro','Markersize',10,'linewidth',1.2);
   plot(Below95(:,1),Below95(:,2),'rd','Markersize',10,'linewidth',1.2);
 
   hBand=area(x,[yUB;yLB-yUB]'); 
    set(hBand,'LineStyle','none','HandleVisibility','off');
    set(hBand(1),'FaceColor','none'); 
    set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
    uistack(hBand,'bottom');
     
    figName='Fig4b_KGEm_vs_NSE';
    formatFig(gcf,axPos,IntX,IntY,'NSE^{ }',Eq4LbTxt,'b)',legTxt);
    
 %%  ----- 
    
NSE=AUXDATA(:,colNSE); NSE(DiscardedId)=nan;
fprintf(1,'\n Number of total cases considered after discarding cases lying outside UB: %i of %i\n\n',...
     sum(~isnan(NSE)),sum(~isnan(AUXDATA(:,colNSE)))); 
 
 fprintf(1,'Beta* of %i selected cases: %1.2f - %1.2f\n',...
     sum(~isnan(NSE)),min(beta(~isnan(NSE))),max(beta(~isnan(NSE))));
 fprintf(1,'Beta* of %i discarded cases: %1.2f - %1.2f\n',...
     sum(isnan(NSE)),min(beta(isnan(NSE))),max(beta(isnan(NSE))));
 
  fprintf(1,'Beta_n of %i selected cases: %1.2f - %1.2f\n',...
     sum(~isnan(NSE)),min(beta(~isnan(NSE)))-1,max(beta(~isnan(NSE)))-1);
 fprintf(1,'Beta_n of %i discarded cases: %1.2f - %1.2f\n',...
     sum(isnan(NSE)),min(beta(isnan(NSE)))-1,max(beta(isnan(NSE)))-1);
 disp(' ');
 
%%%%  Fig 4c  %%%%
   disp('KGE* vs. Eq.7'); IntX=[.2 1]; IntY=IntX; YoYp=[];
   figPos(1)=950; figPos(2)=558; 
   fig3=figure('Position',figPos,'Color','w');
     for ss=1:MaxGroupNr
       selected=find(SetIndex==ss);
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
       x=NSE(selected); eval(['y=' feq ';']); YoYp=[YoYp; AUXDATA(selected,colKGEm) y];
       plot(y,AUXDATA(selected,colKGEm),'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      legTxt{ss}=sprintf('Set %i',ss);
     end; clear ss y
     plot(IntX,IntY,'k--');
     legTxt{10}='1:1';
     xLbl=Eq7LbTxt;

     x=(0.15:0.01:UBmaxVal); eval(['y=' feq ';']); clear x
     hBand=area(y,[yUB;yLB-yUB]'); 
      set(hBand,'LineStyle','none','HandleVisibility','off');
      set(hBand(1),'FaceColor','none'); 
      set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
      uistack(hBand,'bottom');

     figName='Fig4c_KGEm_vs_Eq7';
     
     for aa=1:length(NrDiscardedSet)
         if NrDiscardedSet(aa)>0, legTxt{aa}=[legTxt{aa} '*']; end
     end; clear aa
     formatFig(gcf,axPos,IntX,IntY,xLbl,Eq4LbTxt,'c)',legTxt);
     
   %%%%  Fig 4d  %%%%  
   disp('KGE vs. Eq.7'); legTxt=cell(4,1); YoYp=[];
   figPos(1)=950; figPos(2)=50; 
   fig3=figure('Position',figPos,'Color','w');
     for ss=6:MaxGroupNr
       selected=find(SetIndex==ss);
       if sum(AUXDATA(selected,colNormzd))==length(selected)
          serColor='k'; MrkrCol='m';
        else
          serColor='b'; MrkrCol='none';
       end
       x=NSE(selected); eval(['y=' feq ';']); YoYp=[YoYp; AUXDATA(selected,colKGE) y];
       plot(y,AUXDATA(selected,colKGE),'Color',serColor,...
        'Marker',ser_symbol(ss),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      legTxt{ss-5}=sprintf('Set %i',ss);
      legTxt{ss-4}='1:1';

     end; clear ss y
     plot(IntX,IntY,'k--');
     x=(0.15:0.01:UBmaxVal); eval(['y=' feq ';']); clear x
     hBand=area(y,[yUB;yLB-yUB]'); 
      set(hBand,'LineStyle','none','HandleVisibility','off');
      set(hBand(1),'FaceColor','none'); 
      set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
      uistack(hBand,'bottom');
  
     xLbl=Eq7LbTxt;
 
     figName='Fig4d_KGE_vs_Eq7';
     legTxt{1}=[legTxt{1} '*'];
     legTxt{2}=[legTxt{2} '*'];
     formatFig(gcf,axPos,IntX,IntY,xLbl,Eq2LbTxt,'d)',legTxt); 
    
     x=NSE; eval(['y=' feq ';']);
     [nse,kge_mod,rmse]=GoF_indices([AUXDATA(:,colKGEm) y]);
     fprintf(1,' NSE= %1.3f\n KGE*= %1.3f\n RMSE= %1.3f\n\n',nse,kge_mod,rmse);
    
    disp('Model acceptance thresholds: ');  
    x=[.5 .65 .8 .9 1]; eval(['KGEmth=rounddec(' feq ',2);']);
    x=KGEmth; eval(['NSEth=rounddec(' feq2 ',2);']);
     fprintf(' NSE:\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n',NSEth);
     fprintf(' KGE*:\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n',KGEmth);

 
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

  


function AddText(figh,txt_label,fntsz,xv,yv,align,varargin)
cl='k'; fntwgt='normal';

switch nargin
    case 7, if length(varargin{1})==1, cl=varargin{1}; else fntwgt=varargin{1}; end
    case 8, fntwgt=varargin{1}; cl=varargin{2};
end

try figure(figh);catch, subplot(figh); end
  lim_x=xlim; lim_y=ylim; 
  xpos=(lim_x(2)-lim_x(1))*xv + lim_x(1);
  ypos=(lim_y(2)-lim_y(1))*yv + lim_y(1);
  text(xpos,ypos,txt_label,'Fontsize',fntsz,'Fontweight',fntwgt,'HorizontalAlignment',align,'Color',cl);
end %----
    
     
  function [nse,kge_mod,rmse]=GoF_indices(YoYp)
   noNaNs=~isnan(sum(YoYp,2)); N=sum(noNaNs);
   Yprd=YoYp(noNaNs,2);
   Yobs=YoYp(noNaNs,1);
   r=corrcoef(Yprd,Yobs); r=r(2,1);
   
   alpha=std(Yprd)/std(Yobs);
%   beta=mean(Yprd)/mean(Yobs);
   beta_n=(mean(Yprd)-mean(Yobs))/std(Yobs); % see Clark et al. 2021 eq. 11
   
   %kge=rounddec(1-sqrt((1-r).^2+(1-alpha).^2 + (1-beta).^2),3);  
   kge_mod=rounddec(1-sqrt((1-r).^2+(1-alpha).^2 + (beta_n).^2),3);
   
   SSQ=sum((Yobs-Yprd).^2);
   SS=sum((Yobs-mean(Yobs)).^2);
   nse=1-SSQ/SS; 
   rmse=sqrt(SSQ/length(Yobs));
end%---