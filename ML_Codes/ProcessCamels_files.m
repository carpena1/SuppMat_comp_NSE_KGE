% Extracts observed and simulated values from CAMELS_model_output_daymet/model_output/flow_timeseries/daymet
% and generates the corresponding .in file to be used in FITEVAL

Folders=cell(18,1);
for ff=1:18
   if ff<10, Folders{ff}=['0' num2str(ff)]; else,  Folders{ff}=num2str(ff); end
end; clear ff


for j=4:length(Folders)
disp(['Procesando ' Folders{j} '...']);
folder=[cd '/' Folders{j} '_YoYp']; mkdir(folder);
eval(['cd ' Folders{j} '; files=dir(' char(39) '*_model_output*.txt' char(39) ');']); % Create list of .doc files
T=struct2table(files); % convert files to table
fNm=T.name; % extract column with file names
clear T files
n=length(fNm); headerlines=1;
 for i=1:n
   selNm=fNm{i}; outfNm=strrep(selNm,'_model_output.txt','.in'); 
   strg=['! sed ' char(39) '1,' num2str(headerlines) 'd' char(39) ' ' selNm ' >temp.in'];
  eval(strg);
  AUXDATA=load('temp.in'); delete('temp.in');
  if ~isempty(AUXDATA)
 % writeArray([AUXDATA(:,end) AUXDATA(:,end-1)],8,outfNm,1,{'Yo','Yp'});
   writeArray([AUXDATA(:,end) AUXDATA(:,end-1)],8,outfNm);
  end
  clear AUXDATA
 % movefile(outfNm,[folder '/' outfNm]);
 end; clear i
 movefile('*.in',folder);
 cd ..
end
clear headerlines folder Folders selNm outfNm

