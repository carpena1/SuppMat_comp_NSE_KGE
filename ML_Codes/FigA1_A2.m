close all; clear 
global  Mksz fntSz ser_color ser_symbol figExt figRes figName
ser_color=['b' 'm' 'r' 'k' 'g' 'r' 'b' 'k' 'm' 'g' 'b' 'm']; % Define color for the series to plot in each panel
ser_symbol=['o' '^' 'd' 'v' '<' '>' 'p' 'h' 's' '+' '*' 'x']; % Define symbol for the series to plot in each panel

clear DATA AUXDATA CAMELS

% Figure resolution and extension, marker size, font size
figRes=300; figExt='png'; Mksz=5; fntSz=16; 
 colNSE=4;  colKGE=5; 
 colmYo=7; 
 colmYp=8; 
 colsdYo=9; 
 colsdYp=10;
 col_r=11;
 colRBias=13;
 colKGEm=14;
 colId=16;
 colElev=17;
 colSlope=18;
 colArea=19;
 colRO_ratio=20;
 colBF_index=21;
 colStreamElas=22;

 feq='1 - sqrt(2)*(1-sqrt(x))';
    func_L='1-sqrt(1.911*(x)-3.909*sqrt(x)+2)';
    func_U='1-sqrt(2.114*(x)-4.111*sqrt(x)+2)';

%%%% CAMELS DATA SET %%%%%%% 
fprintf('Reading CAMELS dataset ... ');
 fnm='NSEvsKGE_datasets_camels_48.txt'; 
   A=importdata(fnm); CDATA=A.data; clear A fnm
   fprintf(1,'%i rows and %i columns data available\n',size(CDATA)); 
 
    nseth=0.2; auxind=CDATA(:,colNSE)>nseth;
    disp(['Nr. of NSE values <= ' num2str(nseth) ' discarded: ' num2str(length(auxind)-sum(auxind))]);
    CAMELS=CDATA(auxind,:); clear auxind CDATA
    fprintf(1,'Number of total cases after discarding NSE<%1.1f: %i\n',nseth,sum(~isnan(CAMELS(:,colNSE))));
 
   CAMELS(:,colNSE)=rounddec(CAMELS(:,colNSE),5); 
   CAMELS(:,colKGEm)=rounddec(CAMELS(:,colKGEm),3); 

   CAMELS_rho=rounddec(CAMELS(:,col_r),5); % Pearson correlation coefficient
   CAMELS_alpha=CAMELS(:,colsdYp)./CAMELS(:,colsdYo); % Variability component
   CAMELS_beta=1+((CAMELS(:,colmYp)-CAMELS(:,colmYo))./CAMELS(:,colsdYo));
   
   CAMELS_Rho2Alpha=CAMELS_rho./CAMELS_alpha;
   CAMELS_RhoMinusAlpha=CAMELS_rho-CAMELS_alpha;

   KGEm_prd=rounddec(1-sqrt(max(0,(CAMELS_rho-1).^2+2*CAMELS_alpha.*(CAMELS_rho-1)+1-CAMELS(:,colNSE))),3);
   CAMELS(:,end+1)=KGEm_prd; colKGEm_prd=size(CAMELS,2); % KGE* from NSE


 % Discard cases if outside uncertaintybound 
   clear CAMELS_OutsideId CAMELS_InsideId
   KGEm=CAMELS(:,colKGEm);
   x=CAMELS(:,colNSE); eval(['yLB= ' func_L '; yUB= ' func_U ';']);
   eval(['y= ' feq ';']);
   
   CAMELS_InsideId=(KGEm<0.95) & (KGEm>=yLB & KGEm<=yUB);
   CAMELS_OutsideId_95= (KGEm<0.95) & (KGEm<yLB | KGEm>yUB) & CAMELS_Rho2Alpha<0.95;
   CAMELS_OutsideId_105=(KGEm<0.95) & (KGEm<yLB | KGEm>yUB) & CAMELS_Rho2Alpha>1.05;
  
   clear x yLB yUB KGEm
   
ExNm=[5495500 2349900 4059500];
for i=1:length(ExNm)
 SeExId(i)=find(CAMELS(:,colId)==ExNm(i));
end; clear i ExNm

%%% -----------------------------------------
   IntX=[.2 0.95]; IntY=IntX;
   axPos=[0.137 0.135 0.835 0.85];
   figPos=[10 558 440 390];
   Eq6LbTxt='KGE* (Eq.6)';
   Eq4LbTxt='KGE* (Eq. 4)';
   Eq7LbTxt='KGE* (Eq. 7)';
   func_L='1-sqrt(1.9114*(x)-3.9092*sqrt(x)+2)';
   func_U='1-sqrt(2.1139*(x)-4.111*sqrt(x)+2)';
   
   UBmaxVal=1;
   x=(0.15:0.01:UBmaxVal); eval(['yLB= ' func_L '; yUB= ' func_U ';']); 
   eval(['y1=' feq ';']); clear x

%%% -----------------------------------------

NrTotalCases=sum(~isnan(CAMELS(:,colNSE)));

vN={'NSE','KGE*','beta_n','rho','alpha','rho/alpha','rho-alpha',...
    'Catchment area','Mean Elevation','Mean Slope',...
    'Runoff ratio','Baseflow index','Streamflow elasticity',...
    'KGEm-UB_upp','UB_low-KGEm'};
TEMP=[CAMELS(:,colNSE) CAMELS(:,colKGEm) CAMELS_beta-1 ...
      CAMELS_rho CAMELS_alpha CAMELS_Rho2Alpha ...
      CAMELS(:,colArea) CAMELS(:,colElev) CAMELS(:,colSlope) ...
      CAMELS(:,colRO_ratio) CAMELS(:,colBF_index) CAMELS(:,colStreamElas) ...
    ];
% rho/alpha
% rho - alpha
% difference of KGE* and upper uncertainty bound (UB_upp)
% difference of lower uncertainty bound (UB_low) and KGE*
%catchment area
%catchment mean elevation
%catchment mean slope
% runoff ratio (ratio of mean daily discharge to mean daily precipitation)
% baseflow index (ratio of mean daily baseflow to mean daily discharge)
% streamflow precipitation elasticity (sensitivity of streamflow to changes in precipitation at the annual time scale)

INSIDE=TEMP(CAMELS_InsideId,:);
OUTSIDE=TEMP(~CAMELS_InsideId,:);
OUTSIDE_95=TEMP(CAMELS_OutsideId_95,:);
OUTSIDE_105=TEMP(CAMELS_OutsideId_105,:);

fprintf(1,'Number of total cases lying INSIDE the UB: %i of %i (%1.1f%%)\n',...
     size(INSIDE,1),NrTotalCases,100*size(INSIDE,1)/NrTotalCases); 
PrintStats(INSIDE,vN);
disp('----------------------------------------');
fprintf(1,'Number of total cases lying OUTSIDE the UB and rho/alpha<0.95: %i of %i (%1.1f%%)\n',...
     size(OUTSIDE,1),NrTotalCases,100*size(OUTSIDE,1)/NrTotalCases); 
PrintStats(OUTSIDE,vN);
disp('----------------------------------------')
disp(' ');
 alpha=OUTSIDE(:,5); rho2alpha=OUTSIDE(:,6);
 fprintf('Cases with rho/alpha < 0.95 and alpha>1: %i\n',sum(alpha>1 & rho2alpha<0.95));
 fprintf('Cases with rho/alpha < 0.95 and alpha<1: %i\n',sum(alpha<1 & rho2alpha<0.95));
  fprintf('Alpha range when rho/alpha > 1.05: %1.3f - %1.3f\n\n',...
      min(alpha(rho2alpha>1.05)),max(alpha(rho2alpha>1.05)));


  %%%%  Fig A1   %%%%   
  fprintf('Figure A1 KGE* vs NSE\n');
  IntX=[0 1]; IntY=IntX;
   figPos(1)=100; 
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);
      serColor=[0 .4 0; 0.5 0.4 0;0.5 0.5 0]; 
      hold all
      hA(1)=plot(INSIDE(:,1),INSIDE(:,2),'.','Markersize',Mksz+4,'linewidth',1.2,'Color',serColor(1,:));
      hA(2)=plot(OUTSIDE_95(:,1),OUTSIDE_95(:,2),'Color',serColor(2,:),...
        'Marker','d','Markersize',Mksz,'linestyle','none');
      hA(1)=plot(OUTSIDE_105(:,1),OUTSIDE_105(:,2),'Color',serColor(3,:),...
        'Marker','o','Markersize',Mksz,'linestyle','none');
      eval(['fplot(@(x) ' feq ',[0 1],' char(39) 'k-' char(39) ');']);
      legTxt{1}='Inside UB';
      legTxt{2}='Outside UB & \rho/\alpha < 0.95';
      legTxt{3}='Outside UB & \rho/\alpha > 1.05';
      legTxt{4}='Eq. 7';

      x=(0.15:0.01:UBmaxVal); eval(['yLB= ' func_L '; yUB= ' func_U ';']);
      
      hBand=area(x,[yUB;yLB-yUB]'); 
       set(hBand,'LineStyle','none','HandleVisibility','off');
       set(hBand(1),'FaceColor','none'); 
       set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
       uistack(hBand,'bottom');
       set(gca,'Layer','top');
       box on
       
   figName='FigA1_CAMELS_KGEm_vs_NSE';
    plot(CAMELS(SeExId,colNSE),CAMELS(SeExId,colKGEm),'k^','Markersize',Mksz,'Markerfacecolor','k');
    formatFig(gcf,axPos,IntX,IntY,'NSE^{ }','KGE*',' ',legTxt,'SE');

   %%%%  Fig A2   %%%%
   fprintf('Figure A2 rho vs alpha\n');
   figPos(1)=600;   
   figure('Position',figPos,'Color','w');
   axes('Position',axPos);

    rho=OUTSIDE_105(:,4); alpha=OUTSIDE_105(:,5); rho2alpha=OUTSIDE_105(:,6);
    scatter(alpha,rho,40,rho2alpha,'filled','o'); hold all
    rho=OUTSIDE_95(:,4); alpha=OUTSIDE_95(:,5); rho2alpha=OUTSIDE_95(:,6);
    scatter(alpha,rho,40,rho2alpha, 'filled','d');
    plot([0.3 1],[0.3 1],'k--');
    yLbl='Correlation component, \rho';
    xLbl='Variability component, \alpha';
    ylim(IntY); xlim(IntX); grid on; box on
      cb=colorbar; 
      cb.Label.FontWeight='bold';  cb.Label.FontSize=fntSz;
      cb.Position=[0.88 0.136 0.044 0.785];
   text(1.38,0.98,'\rho/\alpha','Horizontalalignment','center','Fontsize',fntSz+2,'fontweight','bold');
   text(1.15,0.97,'\rho/\alpha < 0.95','Horizontalalignment','center','Fontsize',fntSz-3,'fontweight','normal');
   text(0.45,0.97,'\rho/\alpha > 1.05','Horizontalalignment','center','Fontsize',fntSz-3,'fontweight','normal');

  figName='FigA2_alpha_rho_rho2alpha';
   formatFig(gcf,[0.125 0.135 0.72 0.85],[0.3 1.3],[0.3 1],xLbl,yLbl,' ','','SE');
   
 
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
  if ~isempty(legendTxt)
    leg1=legend(legendTxt,'Location',LOC,'box','off','Fontsize',fntSz-3);
    leg1.ItemTokenSize = [13,12];
  end
  ylabel(Eq6LbTxt,'Fontsize',fntSz,'fontweight','bold');
  xlabel(xLbTxt,'Fontsize',fntSz,'fontweight','bold');
  AddText(gcf,letter,fntSz+2,-.15,1,'left','bold');
  if ~isempty(figName)
   print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str); figName='';
 end
end%--

function [nse,kge_mod,rmse]=GoF_indices(YoYp)
   noNaNs=~isnan(sum(YoYp,2));
   Yprd=YoYp(noNaNs,2);
   Yobs=YoYp(noNaNs,1);
   r=corrcoef(Yprd,Yobs); r=r(2,1);
   
   alpha=std(Yprd)/std(Yobs);
   beta_n=(mean(Yprd)-mean(Yobs))/std(Yobs); % see Clark et al. 2021 eq. 11
   
   kge_mod=rounddec(1-sqrt((1-r).^2+(1-alpha).^2 + (beta_n).^2),3);
   
   SSQ=sum((Yobs-Yprd).^2);
   SS=sum((Yobs-mean(Yobs)).^2);
   nse=1-SSQ/SS; 
   rmse=sqrt(SSQ/length(Yobs));
end%---

function [medV,MedADV,iqrV,minV,maxV]=ComputeStats(x)
  medV=nanmedian(x);
  MedADV=mad(x(~isnan(x)),1);
  iqrV=iqr(x(~isnan(x)));
  minV=nanmin(x);
  maxV=nanmax(x);  
end%---

function PrintStats(X,varName)
 for cc=1:size(X,2)
     [vmed,vMAD,viqr,vmin,vmax]=ComputeStats(X(:,cc));
     fprintf(1,'  %s: Med %sMAD= %1.2f %s%1.2f; IQR= %1.2f; [%1.3f - %1.3f]\n',varName{cc},...
         char(177),vmed,char(177),vMAD,viqr,vmin,vmax);
 end
end%---