close all; clear 
global  Mksz fntSz ser_color ser_symbol figExt figRes figName
ser_color=['b' 'm' 'r' 'k' 'g' 'r' 'b' 'k' 'm' 'g' 'b' 'm']; % Define color for the series to plot in each panel
ser_symbol=['o' '^' 'd' 'v' '<' '>' 'p' 'h' 's' '+' '*' 'x']; % Define symbol for the series to plot in each panel

clear DATA AUXDATA CAMELS

% Figure resolution and extension, marker size, font size
figRes=300; figExt='png'; Mksz=6; fntSz=16; 
 colNSE=4;  colKGE=5; 
 colmYo=7; 
 colmYp=8; 
 colsdYo=9; 
 colsdYp=10;
 col_r=11;
 col_n=12;
 colRBias=13;
 colKGEm=14;

 feq='1 - sqrt(2)*(1-sqrt(x))';
    func_L='1-sqrt(1.911*(x)-3.909*sqrt(x)+2)';
    func_U='1-sqrt(2.114*(x)-4.111*sqrt(x)+2)';

nseth=0.2;
biasth=5; 
ksth=0.05;

%%%% CAMELS DATA SET %%%%%%% 
fprintf('Reading CAMELS dataset ... ');
 fnm='NSEvsKGE_datasets_camels_48.txt'; 
   A=importdata(fnm); CDATA=A.data; clear A fnm
   fprintf(1,'%i rows and %i columns data available\n',size(CDATA)); 
 
    auxind=CDATA(:,colNSE)>nseth;
    disp(['Nr. of NSE values <= ' num2str(nseth) ' discarded: ' num2str(length(auxind)-sum(auxind))]);
    CAMELS=CDATA(auxind,:); clear auxind CDATA
    fprintf(1,'Number of total cases after discarding NSE<%1.1f: %i\n',nseth,sum(~isnan(CAMELS(:,colNSE))));
 
   CAMELS(:,colNSE)=rounddec(CAMELS(:,colNSE),5); 
   CAMELS(:,colKGEm)=rounddec(CAMELS(:,colKGEm),3); 

   CAMELS_rho=rounddec(CAMELS(:,col_r),5); % Pearson correlation coefficient
   CAMELS_alpha=CAMELS(:,colsdYp)./CAMELS(:,colsdYo); % Variability component
   CAMELS_beta=1+((CAMELS(:,colmYp)-CAMELS(:,colmYo))./CAMELS(:,colsdYo));
   
   % Correction du to using sample variance
    CAMELS_n=CAMELS(:,col_n); % number of paired values (n)
    f=(CAMELS_n - 1)./CAMELS_n;
    FF=(1-f).*(CAMELS_alpha.^2-2.*CAMELS_alpha.*CAMELS_rho+CAMELS(:,colNSE));
   
   KGEm_prd=rounddec(1-sqrt(max(0,(CAMELS_rho-1).^2+2*CAMELS_alpha.*(CAMELS_rho-1)+1-CAMELS(:,colNSE)+FF)),3);
   CAMELS(:,end+1)=KGEm_prd; colKGEm_prd=size(CAMELS,2); % KGE* from NSE

   CAMELS_Rho2Alpha=CAMELS_rho./CAMELS_alpha;

    
 CAMELS_BiasedId=find(abs(CAMELS(:,colRBias))>biasth); 
 disp(['Cases with bias > ' num2str(biasth) '%: ' num2str(length(CAMELS_BiasedId))]);
 

 % Discard cases if outside uncertaintybound 
   clear CAMELS_DiscardedId
   KGEm=CAMELS(:,colKGEm);
   x=CAMELS(:,colNSE); eval(['yLB= ' func_L '; yUB= ' func_U ';']);
   CAMELS_DiscardedId=and(KGEm<0.95,or(KGEm>yUB,KGEm<yLB));
   clear x yLB yUB KGEm
   fprintf(1,'CAMELS Cases (KGE*<0.95) lying outside de Uncertainty band: %i\n\n',sum(CAMELS_DiscardedId));
 
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

    
   %%%%  Fig 5a  %%%%  
   disp('KGE* vs Eq.6');
   figPos(1)=950; figPos(2)=550; 
   fig5a=figure('Position',figPos,'Color','w');
     serColor=[0 .5 0]; MrkrCol='none';
      plot(CAMELS(:,colKGEm_prd),CAMELS(:,colKGEm),'Color',serColor,...
        'Marker',ser_symbol(1),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
      plot(IntX,IntY,'k--');

      figName='Fig5a_KGE_vs_Eq6';
      fprintf(1,' Number of total cases considered in %s: %i\n',...
      figName(1:5),sum(~isnan(CAMELS(:,colKGEm))));
      legTxt={'Set 10'};

      formatFig(gcf,axPos,IntX,IntY,Eq6LbTxt,Eq4LbTxt,'a)',legTxt);
     
     disp(' ');

     fprintf(1,' Set %i: %s [%1.3f - %1.3f]\n',10,'NSE',min(CAMELS(:,colNSE)),max(CAMELS(:,colNSE)));
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]\n',10,'KGE',min(CAMELS(:,colKGE)),max(CAMELS(:,colKGE)));
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]\n\n',10,'KGE*',min(CAMELS(:,colKGEm)),max(CAMELS(:,colKGEm)));

     
NSE=CAMELS(:,colNSE); NSE(CAMELS_DiscardedId)=nan;
NrSelCases=sum(~isnan(NSE)); NrTotalCases=sum(~isnan(CAMELS(:,colNSE)));
fprintf(1,' Number of total cases considered (lying inside UB): %i of %i (%1.1f%%)\n\n',...
     NrSelCases,NrTotalCases,100*NrSelCases/NrTotalCases); 


%%%%  Fig 5b  %%%%
   disp('KGE* vs. Eq.7');
   figPos(1)=950; figPos(2)=50; 
   fig5b=figure('Position',figPos,'Color','w');
       serColor=[0 .5 0]; MrkrCol='none';
       x=NSE; eval(['y=' feq ';']);
       plot(y,CAMELS(:,colKGEm),'Color',serColor,...
        'Marker',ser_symbol(1),'Markerfacecolor',MrkrCol,'linestyle','none');hold all
     legTxt{1}=sprintf('Set 10*');
     plot(IntX,IntY,'k--');
     
     hBand=area(y1,[yUB;yLB-yUB]'); 
       set(hBand,'LineStyle','none','HandleVisibility','off');
       set(hBand(1),'FaceColor','none'); 
       set(hBand(2),'FaceColor',[.87 0.92 0.98]); 
       uistack(hBand,'bottom');
       set(gca,'Layer','top');
   
     figName='Fig5b_KGEm_vs_Eq6';
     fprintf(1,' Number of total cases considered in %s: %i (%1.1f%%)\n',...
     figName(1:5),sum(~isnan(NSE)),100*NrSelCases/NrTotalCases);

     formatFig(gcf,axPos,IntX,IntY,Eq7LbTxt,Eq4LbTxt,'b)',legTxt);

     [vmed,viqr,vmin,vmax]=ComputeStats(NSE);
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'NSE',...
         vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS(~isnan(NSE),colKGEm));
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'KGE*',...
         vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_beta(~isnan(NSE))-1);
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'beta_n',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_rho(~isnan(NSE)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'rho',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_alpha(~isnan(NSE)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'alpha',...
        vmin,vmax,vmed,viqr);
     disp(' '); 
     
     NSED=CAMELS(:,colNSE); NSED(~isnan(NSE))=nan; [vmed,viqr,vmin,vmax]=ComputeStats(NSED);
     fprintf(1,' Number of total cases DISCARDED: %i (%1.1f%%)\n',sum(~isnan(NSED)),100*sum(~isnan(NSED))/NrTotalCases);
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'NSE',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS(~isnan(NSED),colKGEm));
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'KGE*',...
         vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_beta(~isnan(NSED))-1);
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'beta_n',...
       vmin,vmax,vmed,viqr);
        [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_rho(~isnan(NSED)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'rho',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_alpha(~isnan(NSED)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'alpha',...
        vmin,vmax,vmed,viqr);

     disp(' '); clear dCAMELS NSED dbeta_n
     
     
[nse,kge_mod,rmse]=GoF_indices([CAMELS(:,colKGEm) y]);
 fprintf(1,' NSE= %1.3f\n KGE*= %1.3f\n RMSE= %1.3f\n\n',nse,kge_mod,rmse);

      NSE=CAMELS(:,colNSE); [vmed,viqr,vmin,vmax]=ComputeStats(NSE);
     fprintf(1,' Number of total cases: %i (%1.1f%%)\n',sum(~isnan(NSE)),100*sum(~isnan(NSE))/NrTotalCases);
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'NSE',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS(~isnan(NSE),colKGEm));
     fprintf(1,' Set %i: %s [%1.3f - %1.3f]; Median= %1.3f; IQR= %1.3f\n',10,'KGE*',...
         vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_beta(~isnan(NSE))-1);
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'beta_n',...
       vmin,vmax,vmed,viqr);
   
           [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_rho(~isnan(NSE)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'rho',...
        vmin,vmax,vmed,viqr);
     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_alpha(~isnan(NSE)));
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'alpha',...
       vmin,vmax,vmed,viqr);

     [vmed,viqr,vmin,vmax]=ComputeStats(CAMELS_Rho2Alpha);
     fprintf(1,' Set %i: %s [%1.4f - %1.4f]; Median= %1.4f; IQR= %1.4f\n',10,'ks=rho/alpha',...
        vmin,vmax,vmed,viqr);

     disp(' '); clear dCAMELS NSE dbeta_n

     
 
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


function formatFig(hFig,axPos,rngX,rngY,xLbTxt,yLbTxt,varargin)
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
  ylabel(yLbTxt,'Fontsize',fntSz,'fontweight','bold');
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

function [medV,iqrV,minV,maxV]=ComputeStats(x)
  medV=nanmedian(x);
  iqrV=iqr(x(~isnan(x)));
  minV=nanmin(x);
  maxV=nanmax(x);
  
end%---
