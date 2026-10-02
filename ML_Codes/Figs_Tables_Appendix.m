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
 col_n=12;
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
   
   % Correction du to using sample variance
    CAMELS_n=CAMELS(:,col_n); % number of paired values (n)
    f=(CAMELS_n - 1)./CAMELS_n;
    FF=(1-f).*(CAMELS_alpha.^2-2.*CAMELS_alpha.*CAMELS_rho+CAMELS(:,colNSE));
   
   KGEm_prd=rounddec(1-sqrt(max(0,(CAMELS_rho-1).^2+2*CAMELS_alpha.*(CAMELS_rho-1)+1-CAMELS(:,colNSE)+FF)),3);
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

vN={'NSE','KGE*','beta_n','rho','alpha','rho/alpha',...
    'Catchment area','Mean Elevation','Mean Slope',...
    'Runoff ratio','Baseflow index','Streamflow elasticity',...
    'KGEm-UB_upp','UB_low-KGEm'};
TEMP=[CAMELS(:,colNSE) CAMELS(:,colKGEm) CAMELS_beta-1 ...
      CAMELS_rho CAMELS_alpha CAMELS_Rho2Alpha ...
      CAMELS(:,colArea) CAMELS(:,colElev) CAMELS(:,colSlope) ...
      CAMELS(:,colRO_ratio) CAMELS(:,colBF_index) CAMELS(:,colStreamElas) ...
     ];

INSIDE=TEMP(CAMELS_InsideId,:);
OUTSIDE=TEMP(~CAMELS_InsideId,:);
OUTSIDE_95=TEMP(CAMELS_OutsideId_95,:);
OUTSIDE_105=TEMP(CAMELS_OutsideId_105,:);
disp('--------------------------------------------------------------------------------');
fprintf(1,'Number of total cases lying INSIDE the UB: %i of %i (%1.1f%%)\n',...
     size(INSIDE,1),NrTotalCases,100*size(INSIDE,1)/NrTotalCases); 
PrintStats(INSIDE,vN);
disp('--------------------------------------------------------------------------------');
fprintf(1,'Number of total cases lying OUTSIDE the UB and rho/alpha<0.95: %i of %i (%1.1f%%)\n',...
     size(OUTSIDE,1),NrTotalCases,100*size(OUTSIDE,1)/NrTotalCases); 
PrintStats(OUTSIDE,vN);
disp('--------------------------------------------------------------------------------')
disp(' ');
 alpha=OUTSIDE(:,5); rho2alpha=OUTSIDE(:,6);
 fprintf('Cases with rho/alpha < 0.95 and alpha>1: %i\n',sum(alpha>1 & rho2alpha<0.95));
 fprintf('Cases with rho/alpha < 0.95 and alpha<1: %i\n',sum(alpha<1 & rho2alpha<0.95));
  fprintf('Alpha range when rho/alpha > 1.05: %1.3f - %1.3f\n\n',...
      min(alpha(rho2alpha>1.05)),max(alpha(rho2alpha>1.05)));

%Normalize |beta_n| to a reasonable marker-size range
   size_from=@(b) 15 + 120 * mat2gray(abs(b), [0 0.35]); % 15 - 120
%Convert them to marker sizes with the same function for the data
   beta_legend_vals=[min(abs(OUTSIDE(:,3))) 0.05 0.1 0.2 max(abs(OUTSIDE(:,3)))];  
   s=size_from(beta_legend_vals);

  %%%%  Fig A1   %%%%   
  fprintf('Figure A1 KGE* vs NSE\n');
  IntX=[0 1]; IntY=IntX;
   figPos(1)=60; 
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);
      serColor=[0 .4 0; 0.5 0.4 0;0.5 0.5 0]; 
      hold all
      plot(INSIDE(:,1),INSIDE(:,2),'.','Markersize',Mksz+4,'linewidth',1.2,'Color',serColor(1,:));
      scatter(OUTSIDE_95(:,1),OUTSIDE_95(:,2),size_from(OUTSIDE_95(:,3)),'d',...
              'MarkerEdgeColor',serColor(2,:));
      scatter(OUTSIDE_105(:,1),OUTSIDE_105(:,2),size_from(OUTSIDE_105(:,3)),'o',...
              'MarkerEdgeColor',serColor(3,:));      
      eval(['fplot(@(x) ' feq ',[0 1],' char(39) 'k-' char(39) ');']);
      plot(CAMELS(SeExId,colNSE),CAMELS(SeExId,colKGEm),'k^','Markersize',Mksz,'Markerfacecolor','k');
      legTxt{1}='Inside UB';
      legTxt{2}='Outside UB & \rho/\alpha < 0.95';
      legTxt{3}='Outside UB & \rho/\alpha > 1.05';
      legTxt{4}='Eq. 7';
      legTxt{5}='Selected cases';

      x=(0.15:0.01:UBmaxVal); eval(['yLB= ' func_L '; yUB= ' func_U ';']);
      
      hBand=area(x,[yUB;yLB-yUB]'); 
       set(hBand,'LineStyle','none','HandleVisibility','off');
       set(hBand(1),'FaceColor','none'); 
       set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
       uistack(hBand,'bottom');
       set(gca,'Layer','top');
       box on
       
   figName='FigA1_CAMELS_KGEm_vs_NSE';
 
   y=0.27;
    for i=1:numel(beta_legend_vals)
        x=1.04-.3; 
        scatter(x,y,s(i),'d','markeredgecolor',serColor(2,:),'linewidth',.6); 
        scatter(x-0.05,y,s(i),'o','markeredgecolor',serColor(2,:),'linewidth',.6); 

        text(gca,x+0.05,y, ...
             sprintf('|\\beta_n|= %.2f', beta_legend_vals(i)), ...
             'FontSize', fntSz-4, 'VerticalAlignment','middle');
         y=y-.055;
    end
     formatFig(gcf,axPos,IntX,IntY,'NSE^{ }','KGE*',' ',legTxt,'NW');

   %%%%  Fig A3   %%%%
   fprintf('Figure A3 rho vs alpha\n');
   figPos(1)=520;   
   figure('Position',figPos,'Color','w');
   axes('Position',axPos);

    rho=OUTSIDE_105(:,4); alpha=OUTSIDE_105(:,5); rho2alpha=OUTSIDE_105(:,6); beta_n=OUTSIDE_105(:,3);
    scatter(alpha,rho,size_from(beta_n),rho2alpha,'filled','o','markeredgecolor','k'); hold all
    rho=OUTSIDE_95(:,4); alpha=OUTSIDE_95(:,5); rho2alpha=OUTSIDE_95(:,6); beta_n=OUTSIDE_95(:,3);
    scatter(alpha,rho,size_from(beta_n),rho2alpha,'filled','d','markeredgecolor','k');
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

   y=0.5;
   rectangle(gca,'Position',[0.955 0.315 .325 .215], ...
          'FaceColor',[1 1 1 0.85], 'EdgeColor',[0.6 0.6 0.6]);
    for i=1:numel(beta_legend_vals)
        x=1.04; 
        scatter(x,y,s(i),'d','markeredgecolor','k','linewidth',.6); 
        scatter(x-0.05,y,s(i),'o','markeredgecolor','k','linewidth',.6); 

        text(gca,x+0.05,y, ...
             sprintf('|\\beta_n|= %.2f', beta_legend_vals(i)), ...
             'FontSize', fntSz-4, 'VerticalAlignment','middle');
         y=y-.04;
    end
   
  figName='FigA3_alpha_rho_rho2alpha';
   formatFig(gcf,[0.125 0.135 0.72 0.85],[0.3 1.3],[0.3 1],xLbl,yLbl,' ','','SE');
   
% Table of performance metrics stratified by NSE intervals
  AllData=[INSIDE; OUTSIDE]; 
  isOutside=[false(size(INSIDE,1),1); true(size(OUTSIDE,1),1)];
  isInside=[true(size(INSIDE,1),1); false(size(OUTSIDE,1),1)];

  NSE=AllData(:,1); 
  KGEm= AllData(:,2);
  x=NSE; eval(['KGEmEq7=' feq ';']);
  beta_n= AllData(:,3);
  rho=AllData(:,4);
  alpha=AllData(:,5);
  r2a=AllData(:,6); % rho/alpha
  abs1mr2a=abs(1-r2a); % abs(1 - rho/alpha)

  edges=[0 0.65 0.80 1.01]; 
  labels={['0.2 ' char(45) ' 0.65)'],['[0.65 ' char(45) ' 0.80)'],['[0.80 ' char(45) ' 1.00)'], 'All cases'};

  nBins=numel(edges)-1;
  N=zeros(nBins+1,1);
  pctOut=zeros(nBins+1,1);
  rmse_Eq7=zeros(nBins+1,1);
  medAbsB=zeros(nBins+1,1); iqrAbsB=zeros(nBins+1,2);
  medRho=zeros(nBins+1,1); iqrRho=zeros(nBins+1,2);
  medAlpha=zeros(nBins+1,1); iqrAlpha= zeros(nBins+1,2);
  medR2A=zeros(nBins+1,1); iqrR2A=zeros(nBins+1,2);
  med1mR2A=zeros(nBins+1,1); iqr1mR2A=zeros(nBins+1,2);

  for i=1:nBins+1
    if i<=nBins
      idx=NSE>=edges(i) & NSE<edges(i+1);
    else
      idx=true(size(NSE)); 
    end

   N(i)=sum(idx);
   [~,~,rmse]=GoF_indices([KGEm(idx) KGEmEq7(idx)]);
   rmse_Eq7(i)=rmse;
   pctOut(i)=100 * sum(isOutside(idx)) / N(i);

   getMedIQR=@(x) deal(median(x), quantile(x,[0.25 0.75]));

   [medAbsB(i), iqrAbsB(i,:)]=getMedIQR(abs(beta_n(idx)));
   [medRho(i),  iqrRho(i,:)]=getMedIQR(rho(idx));
   [medAlpha(i),iqrAlpha(i,:)]= getMedIQR(alpha(idx));
   [medR2A(i),  iqrR2A(i,:)]=getMedIQR(r2a(idx));
   [med1mR2A(i),  iqr1mR2A(i,:)]=getMedIQR(abs1mr2a(idx));
  end

  fprintf('\nTable A1. Summary statistics of the CAMELS simulations grouped by NSE class.\n');
  fprintf('%s\n', repmat('=',1,110));
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'NSE class','N','% outside UB','RMSE Eq.(7)','|beta_n|','rho','alpha','ks','|1-ks|');
  fprintf('%s\n', repmat('-',1,110));

  for i=1:nBins+1
    fprintf('%s\t%5d\t%1.0f%%\t%1.3f\t%1.3f (%4.3f %s %4.3f)\t%1.2f (%4.2f %s %4.2f)\t%1.2f (%4.2f %s %4.2f)\t%1.2f (%4.2f %s %4.2f)\t%1.2f (%4.2f %s %4.2f)\n', ...
      labels{i}, N(i), pctOut(i), rmse_Eq7(i),...  
      medAbsB(i), iqrAbsB(i,1), char(45), iqrAbsB(i,2),...
      medRho(i),  iqrRho(i,1), char(45),  iqrRho(i,2),...
      medAlpha(i),  iqrAlpha(i,1), char(45),  iqrAlpha(i,2),...
      medR2A(i),  iqrR2A(i,1), char(45),  iqrR2A(i,2),...
      med1mR2A(i),  iqr1mR2A(i,1), char(45),  iqr1mR2A(i,2));
  end
  fprintf('%s\n\n', repmat('=',1,110)); 
   

  binEdges=[0 0.05 0.10 Inf];
    Nmat=zeros(3,3); NOUTmat= zeros(3,3); pct =zeros(3,3);
    for i=1:3  
      for j=1:3 
         idx=abs(beta_n)  >= binEdges(j) & abs(beta_n)  < binEdges(j+1) & ...
             abs1mr2a >= binEdges(i) & abs1mr2a < binEdges(i+1);
          Nmat(i,j)=sum(idx);
          if Nmat(i,j) > 0
              NOUTmat(i,j)=sum(isOutside(idx));
              pct(i,j)=100*NOUTmat(i,j)/Nmat(i,j);
          end
      end
    end

    Nrow=sum(Nmat,2);
    Ncol=sum(Nmat,1);
    Nall=sum(Nmat(:));

    pctRow=zeros(3,1);
    for i=1:3
      idx=abs1mr2a >= binEdges(i) & abs1mr2a < binEdges(i+1);
      pctRow(i)=100 * sum(isOutside(idx))/sum(idx);
    end

    pctCol=zeros(1,3);
    for j=1:3
      idx=abs(beta_n) >= binEdges(j) & abs(beta_n) < binEdges(j+1);
      pctCol(j)=100 * sum(isOutside(idx))/sum(idx);
    end
    pctAll=100 * sum(isOutside)/numel(isOutside);

    fprintf('\nTable A2. Number of cases outside the uncertainty band / total number of cases (percentage),\n');
    fprintf('jointly stratified by |1 %s ks| and |beta_n|.\n',char(45));
    fprintf('%s\n',repmat('=',1,100));
    fprintf('  %s     %s\t%s\t%s\t%s\n',...
            ['|1 ' char(45) ' ks|'],'|beta_n| <= 0.05','0.05 < |beta_n| <= 0.10','|beta_n| > 0.10','All');
    fprintf('%s\n',repmat('-',1,90));

    rowLabels={'  <= 0.05    ',['(0.05 ' char(45) ' 0.10]'],'   > 0.10    '};
    for i=1:3
        fprintf('%s\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\n',...
             rowLabels{i},...
            NOUTmat(i,1),Nmat(i,1),pct(i,1),...
            NOUTmat(i,2),Nmat(i,2),pct(i,2),...
            NOUTmat(i,3),Nmat(i,3),pct(i,3),...
            sum(NOUTmat(i,:)),  Nrow(i),  pctRow(i));
    end
    CASESOUT=Nmat.*pct;
    
    fprintf('%s\n',repmat('-',1,90));
     fprintf('%s\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\t%d / %d (%1.0f%%)\n',...
            '    All    ',...
            sum(NOUTmat(:,1)),Ncol(1),pctCol(1),...
            sum(NOUTmat(:,2)),Ncol(2),pctCol(2),...
            sum(NOUTmat(:,3)),Ncol(3),pctCol(3),...
            sum(sum(NOUTmat)),Nall,pctAll);
     fprintf('%s\n\n',repmat('=',1,100));
   
   %%%%  Fig A2   %%%%
   fprintf('Figure A2 beta_n sensitivity\n');
   figPos(1)=975;   
   figure('Position',figPos,'Color','w');
   axes('Position',axPos);
   hold on; box on; grid on;

   eval(['LB=@(x) ' func_L ';']);  
   eval(['UB=@(x) ' func_U ';']);

   nse_vec=linspace(0.50,0.95,91);
   ks_list=[0.975 1.000 1.025];
   cols=[0.165 0.471 0.839; 0.322 0.322 0.306; 0.890 0.286 0.282];
   lspec={'--','-','-.'};
   legStr={'k_s= 0.975','k_s= 1.000','k_s= 1.025'};

   X=[];
   for ik=1:3
     ks=ks_list(ik);
     beta_tol=nan(size(nse_vec));
     
     for i=1:numel(nse_vec)
        nse=nse_vec(i);
        L=LB(nse);       
        for b=0:0.001:0.80
           a=sqrt((nse+b^2)/(2*ks - 1));
           r=ks*a;
           trueK=1-sqrt(b^2+(a-1)^2+(r-1)^2);
           if trueK<L,beta_tol(i)=b; break; end
        end
     end
     plot(nse_vec,beta_tol,'Color',cols(ik,:),'LineStyle',lspec{ik},'LineWidth',2);
   end

   for v=[0.50 0.65 0.80 0.90]
      L=LB(v);
      t=NaN;
      for b=0:0.001:0.80
         a=sqrt((v + b^2));
         r=a;
         trueK=1-sqrt(b^2+(a-1)^2+(r-1)^2);
         if trueK<L; t=b; break; end
      end
      plot(v,t,'o','MarkerSize',7,'LineWidth',1.5,...
         'MarkerFaceColor',[0.32 0.32 0.31], 'MarkerEdgeColor','w');
      if (v==.5 || v==.65), offs=0.012; else, offs=-0.012; end
      text(v-0.014, t+offs, sprintf('%.2f',t), 'FontSize',fntSz-3);
   end

  Qth=quantile(abs(beta_n(isOutside)),[0.25 0.5 0.75]);
    patch([0.45 1 1 0.4],[Qth(1) Qth(1) Qth(3) Qth(3)],[0.85 0.85 0.85],...
    'EdgeColor', 'none','FaceAlpha', 0.35);

   xp=xlim;
   plot([.45 1],Qth(2)*[1 1],':','Color',0.1*[1 1 1],'LineWidth',2);
   text(0.482, Qth(2),sprintf('%.3f ',Qth(2)),'FontSize',fntSz-3,'Color',.3*[1 1 1],...
        'BackgroundColor','none','Horizontalalignment','right');


   axPos(1)=axPos(1)+.01;axPos(3)=axPos(3)+.01;
    
   figName='FigA2_bias_tolerance';
   legA3=formatFig(gcf,axPos,[0.5 0.96],[0 0.35],'NSE^{ }',...
    '|\beta_n| at which KGE* leaves the \pm5% band',' ',legStr,'NE');
    legA3.ItemTokenSize=[25,12]; xlim([.48 0.96]);
    set(legA3,'box','on','Color','w','Edgecolor','w');
     ax=gca; set(ax,'YMinorGrid','on','YMinorTick','on');
     ax.YRuler.MinorTickValues=(0:0.05:0.35);
     figName='FigA2_bias_tolerance';
     print(gcf,figName,['-d' figExt],['-r' num2str(figRes)]);
 

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


function leg1=formatFig(hFig,axPos,rngX,rngY,xLbTxt,Eq6LbTxt,varargin)
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
    leg1.ItemTokenSize=[13,12];
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
     if cc>3
     fprintf(1,'  %s: Med %sMAD= %1.2f %s%1.2f; IQR= %1.2f; [%1.2f - %1.2f]\n',varName{cc},...
         char(177),vmed,char(177),vMAD,viqr,vmin,vmax);
     else
     fprintf(1,'  %s: Med %sMAD= %1.3f %s%1.3f; IQR= %1.3f; [%1.3f - %1.3f]\n',varName{cc},...
         char(177),vmed,char(177),vMAD,viqr,vmin,vmax);
     end
     if cc==6,disp(' '); end
 end
end%---