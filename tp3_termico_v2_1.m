clear all;close all; clc
% SPT - TP3 - Antonio Almeida - 64960
%% --------------------------------------------------
% - Dados - Criar Geometria
%----------------------------------------------------
clear all;close all;format compact; clc

%--------------------------------------------------------------------------
disp('  Criar Geometria')
disp('  -- Termica 2D  --')
%--------------------------------------------------------------------------
nome = input('Nome do ficheiro de resultados -> ','s');
tipo = 'termico2dv2';
nlados = str2double(input('Enter -> Nº lados: ','s'));
%--------------------------------------------------------------------------
fig1 = figure();
%--------------------------------------------------------------------------
k = 1;
keycheck = '';
ptslim = [];
clc;
while k < nlados+1
    keycheck='';
    ok1 = false;
    disp('Coordenadas sequenciais dos vertices')
    disp(strcat('------------------------- Ponto: ',num2str(k)))
    ptslim(1,k) = str2double(input('Enter -> Coordenada x: ','s'));
    ptslim(2,k) = str2double(input('Enter -> Coordenada y: ','s'));
    plot(ptslim(1,:),ptslim(2,:),'bx--')
    disp('-----------------------------------------')
    disp(strcat(['x:';'y:'],num2str(ptslim)))
    keycheck = input('Continue[y]/repeat[n]: ','s');
    clc;
    disp(strcat(['x:';'y:'],num2str(ptslim)))
    disp(strcat('------------------------ Face: ',num2str(k))); 
    axis equal
    if (keycheck == sprintf('y')) 
        k = k + 1;
    else
        clc
        disp('repetir/erro')
    end
end
ptslim(:,nlados+1) = ptslim(:,1);
plot(ptslim(1,:),ptslim(2,:))
%axis([(min(ptslim(1,:))-20),(max(ptslim(1,:))+20),(min(ptslim(2,:))-20),(max(ptslim(1,:))+20)]);
for kk = 1:nlados
    text(((ptslim(1,kk)+ptslim(1,kk+1))/2)+2,((ptslim(2,kk)+ptslim(2,kk+1))/2)+4,sprintf('F%d',kk));
end
axis auto
axis equal
grid on
pause(1)
clc;
%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'nome','tipo','ptslim','nlados');
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% - Dados - Criar Con. fronteira
%--------------------------------------------------------------------------
clear all;close all;format compact; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'termico2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------
tface = zeros(1,nlados);
tcond = '';
tlado = '';
%--------------------------------------------------------------------------
plot(ptslim(1,:),ptslim(2,:))
axis auto
%axis([(min(ptslim(1,:))-20),(max(ptslim(1,:))+20),(min(ptslim(2,:))-20),(max(ptslim(1,:))+20)]);
for kk = 1:nlados
    text(((ptslim(1,kk)+ptslim(1,kk+1))/2)+2,((ptslim(2,kk)+ptslim(2,kk+1))/2)+4,sprintf('F%d',kk));
end
axis equal
grid on

k = 1;
while k <= nlados
    keycheck='';
    clc;
    disp(strcat('------------------------ Face: ',num2str(k))); 
    disp('Cond. fronteira: TºC -> ''t'' / Q -> ''q''');
    tcond(k) = input('Enter ->  ','s');
    tface(k) = str2double(input('Valor: ','s'));
    disp('Lado da fronteira: cima -> ''c'' / baixo -> ''b'' / direita -> ''d'' / esquerda -> ''e''');
    tlado(k) = input('Dir: ','s');
    keycheck = input('Continue[y]/repeat[n]: ','s');
    if (keycheck == sprintf('y')) && ((tcond(k) == 't') || (tcond(k) == 'q')) %&& (((tlado(k) == 'c') || (tlado(k) == 'b') || (tlado(k) == 'e') || (tlado(k) == 'd'))   )
        text(((ptslim(1,k)+ptslim(1,k+1))/2)+2,((ptslim(2,k)+ptslim(2,k+1))/2)+4,sprintf('F%d: %s = %d',k,tcond(k),tface(k)));
        k = k + 1;
    else
        clc
        disp('repetir/erro')
    end
end
%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'nome','tipo','ptslim','nlados','tcond','tface','tlado');
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------
% - Dados - Nº pontos
%--------------------------------------------------------------------------
clear all;close all;format compact; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'termico2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------

matriz1 = [];
disp('Resolução')
ok2 = 'n';
a = 0;
b = 0;
while (a < 3) || (b < 3)
    while ok2 ~= 'y'
        disp('Dim de intervalos (dx)e(dy)')
        dx  = str2double(input('dxy: ','s'));
    %     disp('Dim de intervalos (dy)')
    %     dy  = str2double(input('dy: ','s'));
        dy=dx;
        ok2 = input('Continue/repeat [y]/[n]: ','s');
    end
    dimxx = round(max(ptslim(1,:))-min(ptslim(1,:)));
    dimyy = round(max(ptslim(2,:))-min(ptslim(2,:)));
    npx = round(dimxx/dx);
    npy = round(dimyy/dy);

    ptslimno(1,:) = round((ptslim(1,:)+dx)/dx);
    ptslimno(2,:) = round((ptslim(2,:)+dy)/dy);

    a = abs(ptslimno(1,1:end-1)-ptslimno(1,2:end));
    b = abs(ptslimno(2,1:end-1)-ptslimno(2,2:end));
    a = min(a(a~=0)); if a < 3;disp('!----------------!'); disp('dxy Muito grande não há pontos suficientes');ok2 = 'n';disp('!----------------!');end;
    b = min(b(b~=0)); if b < 3;disp('!----------------!'); disp('dxy Muito grande não há pontos suficientes');ok2 = 'n';disp('!----------------!');end;
end
h  = zeros(npy,npx);
h  = poly2mask(ptslimno(1,:), ptslimno(2,:),npy,npx);  %criar matrizes de valores 1 para pontos de inconita     %max(ptslimno(2,:))-
h  = flipud(h); 
h  = imdilate(h,[0,0,0;1,1,0;1,1,0]);
figure
spy(h)

%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'nome','tipo','ptslim','nlados','tcond','tface','tlado','dx','dy','npx','npy');
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% - Processamento - 
%--------------------------------------------------------------------------
clear all;close all;format compact; clc
%--------------------------------------------------------------------------
tic
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'termico2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------
ptslimno(1,:) = round((ptslim(1,:)+dx)/dx);
ptslimno(2,:) = round((ptslim(2,:)+dy)/dy);
h  = zeros(npy,npx);
H  = h; % matriz com valor dos pontos; (pi+exp(1))*10^10 para inconita, outro valor para t 
h  = poly2mask(ptslimno(1,:),ptslimno(2,:),npy,npx);  %criar matrizes de valores 1 para pontos de inconita  % 
h  = flipud(h);
h  = imdilate(h,[0,0,0;1,1,0;1,1,0]);
%--------------------------------------------------------------------------
% Mostrar figura representativa da malha de pontos
%--------------------------------------------------------------------------
figure
spy(h)
%--------------------------------------------------------------------------
h  = flipud(h);  % corrigir modelo para as coordenadas de matrix
h1 = h;
h  = im2double(h);
%--------------------------------------------------------------------------
% Corrigir coordenadas para não execederem o tamanho do modelo 
%--------------------------------------------------------------------------
ptxx = ptslimno(1,:);
ptyy = ptslimno(2,:);
ptxx(ptxx>npx/2) = ptxx(ptxx>npx/2)-1;
ptyy(ptyy>npy/2) = ptyy(ptyy>npy/2)-1;
%--------------------------------------------------------------------------
% Construir as condições de fronteira
%--------------------------------------------------------------------------
h(h==1) = (pi+exp(1))*10^10;
hh = zeros(npy,npx);
for k = 1:nlados
    xx = (min(ptxx(k),ptxx(k+1)):max(ptxx(k),ptxx(k+1)));
    yy = (min(ptyy(k),ptyy(k+1)):max(ptyy(k),ptyy(k+1)));
    h(yy,xx) = tface(k);
    hh(yy,xx) = k;
end
for k = 1:nlados
    h(ptyy(k),ptxx(k))  = 2*(pi+exp(1))*10^10;
    hh(ptyy(k),ptxx(k)) = 0;
end
H(h==(pi+exp(1))*10^10) = 1:sum(sum((h==(pi+exp(1))*10^10))); % matriz de posições das inconitas
s = max(max(H)); %numero de elementos
h(h==(pi+exp(1))*10^10) = 1;
h(h==2*(pi+exp(1))*10^10) = 1;
%--------------------------------------------------------------------------
% Criar vector b e matrix A ; de Ax=b
%--------------------------------------------------------------------------
b = (zeros(s,1));
A = (zeros(s));
pos = find(H~=0);  %posições das inconitas na mariz
mat1 =  [0  1 0
         1 -4 1
         0  1 0]; %molecula para dx==dy
for k = 1:numel(pos)
    [J,I] = ind2sub(size(h),pos(k));
    matpos = [0        H(J-1,I) 0
              H(J,I-1)  H(J,I)  H(J,I+1)   
              0        H(J+1,I) 0];

    matval = [0         h(J-1,I)  0
              h(J,I-1)   h(J,I)   h(J,I+1)
              0         h(J+1,I)  0];

    matbvl = [0         hh(J-1,I)  0
              hh(J,I-1)  hh(J,I)   hh(J,I+1)
              0         hh(J+1,I)  0];      

          pos1 = (find(matpos>0));
          
    for kkk = 1:numel(pos1) 
        kk = pos1(kkk);
        A(matpos(5),matpos(kk)) = A(matpos(5),matpos(kk))+ mat1(kk)*matval(kk);
    end
       
    for kkk = 1:9
        if matbvl(kkk) ~= 0 &&matpos(kkk)== 0
            switch tcond(matbvl(kkk))
                case 't'
                    bcond = matval(kkk);
                case 'q'
                    bcond = matval(kkk)*(dx^2);
                    A(matpos(5),matpos(5)) = A(matpos(5),matpos(5))+ 1;
            end
            b(matpos(5)) = b(matpos(5)) - bcond;
        end
    end
%      pause
end
disp('Resultados Cálculados - RR')
RR = A\b;
disp('Max(RR)')
R1 = max(abs(RR))
disp('Min(RR)')
R2 = min(abs(RR))
clear A b
%--------------------------------------------------------------------------
% Corrigir valores dos cantos e das faces com condução
%--------------------------------------------------------------------------
chapa = h;
for k = 1:numel(RR)
    chapa(H==k) = RR(k);
end
clear H 
for k = 1:nlados
    if tcond(k) == 'q'
        pos2 = find(hh==k);
        for kk = 1:numel(pos2)
            [J,I] = ind2sub(size(hh),pos2(kk));
            switch tlado(k)
                case 'c'
                    chapa(J,I) = chapa(J-1,I)+(tface(k)*(dy^2));
                case 'b'
                    chapa(J,I) = chapa(J+1,I)+(tface(k)*(dy^2));
                case 'd'
                    chapa(J,I) = chapa(J,I-1)+(tface(k)*(dx^2));
                case 'e'
                    chapa(J,I) = chapa(J,I+1)+(tface(k)*(dx^2));
            end
        end
    end
end
chapa2 = zeros(size(chapa)+2);
chapa2(2:end-1,2:end-1)=chapa;
for k = 1:nlados
    mo =[0, chapa2(ptyy(k),ptxx(k)+1), 0; chapa2(ptyy(k)+1,ptxx(k)), 0, chapa2(ptyy(k)+1,ptxx(k)+2); 0, chapa2(ptyy(k)+2,ptxx(k)+1), 0];
    chapa2(ptyy(k)+1,ptxx(k)+1) = (mo(4)+mo(6)+mo(2)+mo(8))/(sum(sum(mo~=0)));
end

chapa = chapa2(2:end-1,2:end-1);
clear chapa2
disp('Resultados Globais - chapa')
disp('Max(chapa)')
Rmx = max(max(abs(RR)))
disp('Min(chapa)')
Rmn = min(min(abs(RR)))
chapa(h1==0) = nan;
%--------------------------------------------------------------------------
toc
%--------------------------------------------------------------------------
tipo = 'termico2dv2_res';
save(strcat(nome,'_res_t','.mat'),'tipo','chapa');
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Mostrar resultados
%--------------------------------------------------------------------------
clear all;close all;format compact; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'termico2dv2_res'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------

figure
pcolor(chapa)
colorbar
if ((size(chapa,1)*size(chapa,2))>1000)
shading interp
end



