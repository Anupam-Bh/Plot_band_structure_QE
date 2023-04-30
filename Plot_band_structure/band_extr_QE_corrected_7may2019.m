%% Band_extraction_QE
clear all
close all;
clc

%% For multiple files (comment this section if single file)
% vec1=[0.02 0.019 0.018 0.017 0.016];
% vec2=[-0.02 -0.019 -0.018 -0.017 -0.016];
% vec_1s=strsplit(sprintf('%.6f ',vec1),{' ','[',']'});
% vec_2s=strsplit(sprintf('%.6f ',vec2),{' ','[',']'});
% vec_1ss=string(vec_1s);
% vec_2ss=string(vec_2s);
% vec_1ss(end)=[];vec_2ss(end)=[];
% for v1=1:length(vec_1ss)
%     for v2=1:length(vec_2ss)
%         filename_scf=strcat('scf_',vec_1ss(v1),'_',vec_2ss(v2),'.out');
%       lsfilename_band=strcat('bands_',vec_1ss(v1),'_',vec_2ss(v2),'.out');
%         scf=fopen(filename_scf,'r');
%         while feof(scf)==0
%             tline = fgetl(scf);
%             str=sprintf(tline);
%             if strncmp(str,'     the Fermi energy is',21)==1
%                 X=strsplit(str,{' '});
%                 fermi=str2double(X(6))
%             end
%         end
%         fclose(scf);

%% For single file 
        
        
% % Fermi energy (find and enter from SCF/NSCF calculation) (uncomment
% following 2 lines for single file)
fermi=input('Enter the value of Fermi energy: ');
filename_band=input('Enter filename: ','s');

cc=fopen(filename_band,'r');
i=0;
tictoc=0;
frac_lines_start=1;
kpts=0;
while feof(cc)==0
    i=i+1;
    tline = fgetl(cc);
    str=sprintf(tline);
    if strncmp(str,'     number of Kohn-Sham states',25)==1
        X=strsplit(str,{'='});
        num_bands=str2double(X(2));
    elseif strncmp(str,'     number of k points',20)==1
        X=strsplit(str,{'=',' '});
        kpts=str2double(X(6));
    elseif strncmp(str,'     End of band structure calculation',25)==1
        starting_line=i;
    elseif strncmp(str,'                       cryst. coord.',25)==1
        frac_lines_start=i;
        tictoc=1;
    elseif i> frac_lines_start && i<=(frac_lines_start+kpts)&& tictoc==1
        X=strsplit(str,{'(',')',' '});
        kpoints(str2double(X(3)),1:3)=[str2double(X(5)) str2double(X(6)) str2double(X(7))];
    end
end
fclose(cc);
kblock_len=ceil(num_bands/8)+3;
fprintf('starting bands line number : %d\n',starting_line);
fprintf('number of K points : %d\n',kpts);
fprintf('Number of bands : %d\n',num_bands);
fprintf('K block length : %d\n',kblock_len);

%% Read Klines
clear X cc i
n=1;
for K=starting_line+2:kblock_len:((kpts*kblock_len)+starting_line+2)
    Klines(n)=K;
    n=n+1;
end
kpts_array=linspace(1,kpts,kpts);

%% Reading Kpath

% kpath_filename=input('Enter kpath filename: ','s'); (uncomment these)
% ver=fopen(kpath_filename,'r');
ver=fopen('kpath','r');
j=1;
while feof(ver)==0
    tline = fgetl(ver);
    str=sprintf(tline);
    X=strsplit(str,{' '});
    vertices_char(j)=X(1);
    path(j,:)=[str2double(X(2)) str2double(X(3)) str2double(X(4))];
    j=j+1;
end
fclose(ver);
kstart=1;
for j=1:length(vertices_char)
    for jk=kstart:kpts
        if(isequal(path(j,:),kpoints(jk,:))==1)
            vertices_xval(j)=jk;
            kstart=jk+1;
            break;
        end
    end
end

%% read eigenvalues
clear X
eigenvalues=zeros(kpts,num_bands);
for K=1:kpts
    i=0;
    cc=fopen(filename_band,'r');
    while feof(cc)==0
        i=i+1;
        tline = fgetl(cc);
        if i>=(Klines(K)+2) && i<(Klines(K)+kblock_len-1)
            str=sprintf(tline);
            formatspecific=[repmat('%9s',1,9)];
            CCC=textscan(str,formatspecific);   
            for xcount=1:length(CCC)-1
                X(xcount)=str2double(cell2mat(CCC{xcount}))-fermi;
            end
                
            len=length(X);
            line=i-(Klines(K)+2);
            eigenvalues(K,(8*line+1):len+(8*line))=X; % pw.x output is already in eV
            (1);
            
        end
    end
    fclose(cc);
end
                
                
%% Plot
figure(2)

for i=1:num_bands
    plot(kpts_array,transpose(eigenvalues(:,i)),'r','LineWidth',1);
    hold on
end
ylabel('E-E_f (eV)');
set(gca,'xlim',[1 kpts],'ylim',[-1 1],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','on',...
    'Fontweight','normal','Fontsize',26,'Fontname','Heveltica');
pbaspect([1 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
% text(153,0.08,'\Gamma_8','FontSize',15,'FontWeight','bold')
% text(153,-0.48,'\Gamma_6','FontSize',15,'FontWeight','bold')
% figname=strcat('band_',vec_1ss(v1),'_',vec_2ss(v2),'.png');
% print(char(figname),'-dpng','-r400');

%close all;
% clear filename_scf filename_band eigenvalues X figname tline
%     end
% end
