%% Band_extraction_QE
clear all
close all;
clc
fprintf('###### Remember to keep verbosity = ''high'' for the calculation\n')
fprintf('This code works for nspin = 2 \n')
fprintf('works for both calculation = ''bands'' or ''nscf''\n')
fprintf('-----------------------------------------------------\n')
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
    elseif strncmp(str,' ------ SPIN UP ------------',25)==1
        starting_line_up=i;
    elseif strncmp(str,' ------ SPIN DOWN ------------',25)==1
        starting_line_down=i;        
    elseif strncmp(str,'                       cryst. coord.',25)==1
        frac_lines_start=i;
        tictoc=1;
    elseif i> frac_lines_start && i<=(frac_lines_start+kpts)&& tictoc==1
        X=strsplit(str,{'(',')',' '});
        kpoints(str2double(X(3)),1:3)=[str2double(X(5)) str2double(X(6)) str2double(X(7))];
    end
end
fclose(cc);
calculation=input('Enter the calculation (bands/nscf): ','s');

kblock_len=ceil(num_bands/8)+3;
fprintf('starting bands (up) line number : %d\n',starting_line_up);
fprintf('starting bands (down) line number : %d\n',starting_line_down);
fprintf('number of K points : %d\n',kpts);
fprintf('Number of bands : %d\n',num_bands);
fprintf('K block length : %d\n',kblock_len);

%% Read Klines (for 'nscf' and 'bands' calculation)
clear X cc i
n=1;
if strcmp(calculation,'bands')
    for K=starting_line_up+3:kblock_len:((kpts*kblock_len)+starting_line_up+3)
        Klines_up(n)=K;
        n=n+1;
    end
    clear K
    n=1;
    for K=starting_line_down+3:kblock_len:((kpts*kblock_len)+starting_line_down+3)
        Klines_down(n)=K;
        n=n+1;
    end
    clear K
elseif strcmp(calculation,'nscf')
    for K=starting_line_up+3:(kblock_len*2)-1:(((kpts-1)*kblock_len*2)+starting_line_up+3)
        Klines_up(n)=K;
        n=n+1;
    end
    clear K
    n=1;
    for K=starting_line_down+3:(kblock_len*2)-1:(((kpts-1)*kblock_len*2)+starting_line_down+3)
        Klines_down(n)=K;
        n=n+1;
    end
    clear K
end


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
%upspin
clear X
eigenvalues_up=zeros(kpts,num_bands);
for K=1:kpts
    i=0;
    cc=fopen(filename_band,'r');
    while feof(cc)==0
        i=i+1;
        tline = fgetl(cc);
        if i>=(Klines_up(K)+2) && i<(Klines_up(K)+kblock_len-1)
            str=sprintf(tline);
            formatspecific=[repmat('%9s',1,9)];
            CCC=textscan(str,formatspecific);   
            for xcount=1:length(CCC)-1
                X(xcount)=str2double(cell2mat(CCC{xcount}))-fermi;
            end
                
            len=length(X);
            line=i-(Klines_up(K)+2);
            eigenvalues_up(K,(8*line+1):len+(8*line))=X; % pw.x output is already in eV
            (1);
            
        end
    end
    fclose(cc);
end         
%downspin
clear X
eigenvalues_down=zeros(kpts,num_bands);
for K=1:kpts
    i=0;
    cc=fopen(filename_band,'r');
    while feof(cc)==0
        i=i+1;
        tline = fgetl(cc);
        if i>=(Klines_down(K)+2) && i<(Klines_down(K)+kblock_len-1)
            str=sprintf(tline);
            formatspecific=[repmat('%9s',1,9)];
            CCC=textscan(str,formatspecific);   
            for xcount=1:length(CCC)-1
                X(xcount)=str2double(cell2mat(CCC{xcount}))-fermi;
            end
                
            len=length(X);
            line=i-(Klines_down(K)+2);
            eigenvalues_down(K,(8*line+1):len+(8*line))=X; % pw.x output is already in eV
            (1);
            
        end
    end
    fclose(cc);
end 
%% Plot
figure(2)
kpts_array=linspace(1,kpts,kpts);
for i=1:num_bands
    up=plot(kpts_array,transpose(eigenvalues_up(:,i)),'b','LineWidth',1.0);
    hold on
end
hold on
for i=1:num_bands
    down=plot(kpts_array,transpose(eigenvalues_down(:,i)),'--r','LineWidth',1.0);
    hold on
end
legend([up(1) down(1)],'spin up','spin down ','box','off');
ylabel('E-E_f (eV)');
set(gca,'xlim',[1 kpts],'ylim',[-3 2],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','off',...
    'Fontweight','bold','Fontsize',16,'Fontname','Heveltica');
pbaspect([1.5 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
% figname=strcat('band_',vec_1ss(v1),'_',vec_2ss(v2),'.png');
% print(char(figname),'-dpng','-r400');
print('band_structure.png','-dpng','-r400')
%close all;
% clear filename_scf filename_band eigenvalues X figname tline
%     end
% end