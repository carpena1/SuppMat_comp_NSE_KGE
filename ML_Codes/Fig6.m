close all; clear 
global  fntSz  figExt figRes figName
figRes=300; figExt='png'; fntSz=16;
    
%Full space of possible (NSE, KGE*) combinations for all alpha [0.1, 2.5] and rho [0,1]  
  NSE_sp=(0:0.01:1); rho_sp=(0:0.1:1); alpha_sp=(0.1:0.1:2.5);
  ss=1;
  for ii=1:length(NSE_sp)
    NSE=NSE_sp(ii);
    for jj=1:length(rho_sp)
       rho=rho_sp(jj);
       for kk=1:length(alpha_sp)
        alpha=alpha_sp(kk);
        KGEm_sp(ss)=rounddec(1-sqrt((rho-1).^2+2*alpha.*(rho-1)+1-NSE),3);
        NSE_sp_n(ss)=NSE;
        ss=ss+1;
        end; clear kk 
    end; clear jj
  end; clear ii ss 

% --- Eq. (6) scatter: beta*=1, wide (unconstrained) ranges ---
% rho   in [0, 1]: full Pearson correlation range (positive)
% alpha in [0.1, 2.5]: wide variability ratio range
rng(42);   % Set random seed for reproducibility
n = 6000;
rho_s= rand(n,1);                 % Uniform in [0,1]
alpha_s = 0.1 + (2.5-0.1)*rand(n,1); % Uniform in [0.1,2.5]
nse_s = -alpha_s.^2 + 2.*alpha_s.*rho_s; % Eq. (5) with beta* = 1
kge_s = 1 - sqrt((rho_s - 1).^2 + (alpha_s - 1).^2); % Eq. (4) with beta* = 1
 
 
%%% -----------------------------------------
   axPos=[0.137 0.135 0.835 0.85];
   figPos=[10 558 440 390];

   feq7='1 - sqrt(2)*(1-sqrt(x))';
   func_L='1-sqrt(1.9114*(x)-3.9092*sqrt(x)+2)';
   func_U='1-sqrt(2.1139*(x)-4.111*sqrt(x)+2)';
   
   UBmaxVal=1;
   x=(0.15:0.01:UBmaxVal); eval(['yLB= ' func_L '; yUB= ' func_U ';']); 
   eval(['bnd_wdht=' feq7 ';']); clear x

%%% -----------------------------------------
   x=(0:0.01:1); 
   eval(['KGE_eq7= ' feq7 ';']);

   epsilon=0.05; % 1 - ks = 0.05 % 0.95
   phi=sqrt(x/(1-2*epsilon));
   KGE_eq7_UB=1-sqrt(((1-epsilon).*phi-1).^2+(phi-1).^2);

   epsilon=-0.05; % 1 - ks = -0.05 % 1.05
   phi=sqrt(x/(1-2*epsilon));
   KGE_eq7_LB=1-sqrt(((1-epsilon).*phi-1).^2+(phi-1).^2);
 
   disp('KGE* vs NSE');IntX=[0 1]; IntY=[-0.75 1.05];
   figPos(1)=50; figPos(2)=50;  
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);

   scatter(nse_s,kge_s,10,[0.4 0.6 0.8],'filled','MarkerFaceAlpha',0.2);
   hold all
   j=1; legTxt{j}='Eq. (6) [\beta*=1, \rho\in[0,1], \alpha\in[0.1,2.5]]'; j=j+1;

   fill([x fliplr(x)],[KGE_eq7_LB fliplr(KGE_eq7_UB)],[0.93 0.78 0.45],...
    'FaceAlpha',0.45,'EdgeColor','none'); % Uncertainty band
    legTxt{j}=['\pm5% slope uncertainty band' char(10) ' (0.95 \leq k_s \leq 1.05)']; j=j+1;
    hold all
  
   plot(x,KGE_eq7_UB,'--','Color',[0.8 0.4 0.1],'LineWidth',1.5); % ks=0.95
   plot(x,KGE_eq7_LB,':','Color',[0.8 0.4 0.1],'LineWidth',2); % ks=1.05
   plot(x,KGE_eq7,'k','LineWidth',2); % Eq. (7)
   legTxt{j}=['k_s= 0.95 (' char(949) '= +0.05, slope below 1)']; j=j+1;
   legTxt{j}=['k_s= 1.05 (' char(949) '= -0.05, slope above 1)']; j=j+1;
   legTxt{j}='Eq. (7) (k_s= 1)'; 
   
% Reference lines
   yp=ylim; xp=xlim; 
   plot(xp,-0.414*[1 1],'--','Color',[0.7 0.7 0.7]);
   plot(xp,[0 0],'Color',[0.7 0.7 0.7]);
   text(0.05,-0.46,'NSE= 0  \rightarrow  KGE*\approx -0.414',...
      'FontSize',fntSz-4,'Color',[0.45 0.45 0.45]),
   box on
   grid on
   figName='';
   leg1= formatFig(gcf,axPos,IntX,IntY,'NSE^{ }','KGE*','a)',legTxt,'E');
   set(leg1,'FontSize',fntSz-4,'Position',[0.45 0.35 0.51 0.239]);
   set(leg1,'box','on');
   set(gca,'xtick',0:0.2:1,'ytick',-.6:0.2:1);
   leg1.ItemTokenSize = [25,1];
   
figName='Fig6a_KGEm_vs_NSE';
print_str=['print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str);

 
   disp('DKGE* vs NSE');IntX=[0 1]; IntY=[-0.107 0.107];
   figPos(1)=500; figPos(2)=50;  
   figure('Position',figPos,'Color','w');
    axes('Position',axPos);
    serColor=[0 .5 0]; MrkrCol='none';
       
    bnd_wdht=KGE_eq7_UB-KGE_eq7_LB;% band width   
    y2=KGE_eq7_UB - KGE_eq7; % upper deviation KGE* (ks=0.95) - KGE* (Eq. 7)
    y3=KGE_eq7_LB - KGE_eq7; % lower deviation KGE* (ks=1.05) - KGE* (Eq. 7)
    
    plot([x(bnd_wdht>=0) 0.992],[bnd_wdht(bnd_wdht>=0) 0],'-','Color',[1 0.5 0],'LineWidth',2.5); hold all
    plot(x,y2,'--','Color',[0.8 0.4 0.1],'LineWidth',2.2);
    plot(x,y3,':','Color',[0.8 0.4 0.1],'LineWidth',2.5);
    j=1;
    clear legTxt
    legTxt{j}=['Band width [KGE^*(k_s=0.95)   ' char(822) ' KGE^*(k_s=1.05)]']; j=j+1;
    legTxt{j}=['KGE^*(k_s=0.95)   ' char(822) ' KGE^*(Eq. 7) [upper deviation]']; j=j+1;
    legTxt{j}=['KGE^*(k_s=1.05)   ' char(822) ' KGE^*(Eq. 7) [lower deviation]']; j=j+1;
   fill([x(bnd_wdht>=0) 0.992 0],[bnd_wdht(bnd_wdht>=0) 0 0],[0.93 0.78 0.45],...
    'FaceAlpha',0.45,'EdgeColor','none'); hold all % Uncertainty band

    yp=ylim; xp=xlim; 
    plot(xp,[0 0],'Color',[0.7 0.7 0.7]);

   box on
   grid on
   figName='';
   leg1= formatFig(gcf,axPos,IntX,IntY,'NSE^{ }','\DeltaKGE*','b)',legTxt,'SW');
   set(leg1,'FontSize',fntSz-4,'Position',[0.16 0.18 0.65 0.15]);
   set(leg1,'box','on');
   ha=gca;
   set(ha,'xtick',0:0.2:1,'ytick',-.1:0.02:0.1);
   ha.YRuler.TickLabelFormat='%1.2f'; 
   
    pts=[0.50 0.65 0.80 0.90];
    for i=1:length(pts)
      [~,idx] = min(abs(x-pts(i)));
      plot(x(idx),bnd_wdht(idx),'o','MarkerSize',6,'MarkerFaceColor',...
         [1 0.5 0],'MarkerEdgeColor',[1 0.5 0]);
    
      txtPos(i,:)=[x(idx)-0.0,bnd_wdht(idx)+0.0075];
      if x(idx)>=0.875, txtPos(i,:)=[0.99 0.067]; end
     
      text(txtPos(i,1),txtPos(i,2),sprintf('\\Delta=%.3f',bnd_wdht(idx)),...
         'Color',[1 0.5 0],'FontSize',fntSz-4,'Horizontalalignment','right');
      plot(x(idx),0,'k+','MarkerSize',5,'MarkerEdgeColor',[0.7 0.7 0.7]);
      text(x(idx),-.0075,sprintf('%.2f',x(idx)),'Horizontalalignment','center',...
         'Color',[0.7 0.7 0.7],'FontSize',fntSz-5);
    end; clear i
        
figName='Fig6b_DeltaKGEm_vs_NSE';
print_str=['print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str);
    


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
  if ~isempty(legendTxt); leg1=legend(legendTxt,'Location',LOC,'box','off'); end
  ylabel(Eq6LbTxt,'Fontsize',fntSz,'fontweight','bold');
  xlabel(xLbTxt,'Fontsize',fntSz,'fontweight','bold');
  AddText(gcf,letter,fntSz+2,-.15,1,'left','bold');
  if ~isempty(figName)
   print_str=[ 'print(gcf,' char(39) figName char(39) ',' char(39) '-d' figExt char(39) ',' char(39) '-r' num2str(figRes) char(39) ');'];
    eval(print_str); figName='';
 end
end%--
