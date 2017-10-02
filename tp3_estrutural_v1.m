clear all;close all; clc
% SPT - TP3 - Antonio Almeida - 64960
%% ------------------------------------------------------------------------
% - Dados - Criar Geometria
%--------------------------------------------------------------------------
clear all;close all; clc
disp('           MDF Estrutural')
disp('    Função para chapa Rectangular ')
disp('       & dx==dy & q distribuido')  
format compact;
%--------------------------------------------------------------------------
nome = input('Nome do ficheiro de resultados -> ','s');
tipo = 'estrutural2dv2';
%--------------------------------------------------------------------------
disp('Dimensões da chapa:')
lx = str2double(input('Comp. x (m) -> ','s'));
ly = str2double(input('Comp. y (m) -> ','s'));
t  = str2double(input('Esp.  t (m) -> ','s'));
disp('Propriedades da chapa:')
E  = str2double(input('E (GPa)-> ','s'));
v  = str2double(input('v -> ','s'));
qq = str2double(input('Carga q (N)-> ','s'));
D  = (E*(10^9)*t^3)/(12*(1-v^2))
%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'tipo','nome','lx','ly','E','t','v','qq','D')
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% - Dados - Criar Con. fronteira
%--------------------------------------------------------------------------
clear all;close all; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'estrutural2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------
ok  = 0;
ok1 = 1;
ta = '----';
while ok ~= 1
    while ok1 <= 4
        disp('      F1    ')
        disp('   -------- ')
        disp('   |      | ')
        disp('F4 |      | F2')
        disp('   |      | ')
        disp('   -------- ')
        disp('      F3    ')
        disp('--------------------')
        disp(sprintf('Face %d -> %s',1,ta(1)))
        disp(sprintf('Face %d -> %s',2,ta(2)))
        disp(sprintf('Face %d -> %s',3,ta(3)))
        disp(sprintf('Face %d -> %s',4,ta(4)))
        disp('--------------------')
        disp(sprintf('Face %d',ok1))
        disp('Tipo de apoio (Encastrado -> e | Apoio simples -> a  | Livre -> l) ')
        disp(ok1)
        ta(ok1) = input('Apoio -> ','s');
        if ta(ok1) == 'e' ||ta(ok1) == 'a' ||ta(ok1) == 'l' 
            disp(sprintf('Face %d -> %s',ok1,ta(ok1)))
            ok1 = ok1 + 1;
            pause(0.8)
            clc
            disp('--------------------')
        else
            disp('Repetir cond. não válida')
            ta(ok1) = '-';
        end
    end
    if ta(1) == ta(2) && ta(1) == ta(3) && ta(1) == ta(4)
        simetria = 4; % simetrica em 4 dir
        ok =1;
        simetria2 = 'Cond. fronteira simetricas em x e y';
        disp('Cond. fronteira simetricas em x e y')
    elseif ta(4) == ta(2) 
        simetria = 3; % simetria horizontal
         ok =1;
         simetria2 = 'Cond. fronteira simetricas em y';
         disp('Cond. fronteira simetricas em y')
    elseif ta(1) == ta(3)
        simetria = 2; % simetria vertucal
         ok =1;
         simetria2 = 'Cond. fronteira simetricas em x';
         disp('Cond. fronteira simetricas em x')
    else
        simetria = 1; % nao simetrico
         ok =1;
         simetria2 = 'Con. fronteira não simetricas';
         disp('Con. fronteira não simetricas')
    end
    pconf = 0;  % variável para confirmar que os apoios estão ok
    for k = 1:4
        if ta(k)== 'e'; pconf = pconf + 2; elseif ta(k) == 'a'; pconf = pconf + 1; end
    end
    if pconf <= 1; ok = 0; 
        simetria = 0;
        clc
        disp('As extremidades não podem estar todas livres')
        disp('Pelo menos uma tem que estar encastrada ou duas apoiadas')
        disp('Repetir')
    end
end
%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'tipo','nome','lx','ly','E','t','v','qq','D','ta','simetria','simetria2')
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% - Dados - Nº pontos
%--------------------------------------------------------------------------
clear all;close all; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'estrutural2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------
disp(sprintf('%s com: dimx = %d ; dimy = %d',simetria2,lx,ly))

dxy = str2double(input('dxy (m) -> ','s'));
if simetria == 4;
    npx = round((lx/dxy)/2);
    npy = round((ly/dxy)/2);
    npt = npx*npy;
    disp(sprintf('Nº elementos = %d',npt*4))
elseif simetria == 3;
    npx = round((lx/dxy)/2);
    npy = round((ly/dxy));
    npt = npx*npy;
    disp(sprintf('Nº elementos = %d',npt*2))
elseif simetria == 2;
    npx = round((lx/dxy));
    npy = round((ly/dxy)/2);
    npt = npx*npy;
    disp(sprintf('Nº elementos = %d',npt*2))
elseif simetria == 1;
    npx = round((lx/dxy));
    npy = round((ly/dxy));
    npt = npx*npy;
    disp(sprintf('Nº elementos = %d',npt))
end
%--------------------------------------------------------------------------
save((strcat(nome,'.mat')),'tipo','nome','lx','ly','E','t','v','qq','D','ta','simetria','dxy','npx','npy','npt','D','qq','ta','simetria2')
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% - Processamento - 
%--------------------------------------------------------------------------
clear all;close all; clc
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'estrutural2dv2'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------
tic
%--------------------------------------------------------------------------
% criar condições de fronteira e simetrias
tta = ones(1,4);
%--------------------------------------------------------------------------
for k = 1:4;
    if ta(k) == 'e'
        tta(k) = 1;
    elseif ta(k) == 'a'
        tta(k) = -1;
    elseif ta(k) == 'l'
        tta(k) = 0;
    end
end

switch simetria
    case 4                         %simetria 4 dir (horizontal e vertical)
        h = zeros(npy+2,npx+2);
        h(3:end,3:end) = 1;
        h1 = h;
        h0 = h;
        h(h==1)=1:sum(sum(h));
        
        h(3:end,1) = h(3:end,3);    %fronteira 4
        h(1,3:end) = h(3,3:end);    %fronteira 1
        
        h1(:,1) = tta(4);
        h1(1,:) = tta(1);

        h2 = fliplr(h(:,1:end-1));
        h3 = flipud(h(1:end-1,:));
        h4 = flipud(h2(1:end-1,:));

        h12 = fliplr(h1(:,1:end-1));
        h13 = flipud(h1(1:end-1,:));
        h14 = flipud(h12(1:end-1,:));

        h02 = fliplr(h0(:,1:end-1));
        h03 = flipud(h0(1:end-1,:));
        h04 = flipud(h02(1:end-1,:));

        H0 = [h0 h02; h03 h04];
        H = [h h2; h3 h4];
        hh = [h1 h12; h13 h14];   
        
        s = sum(sum(h0));
        pos = (find(h0==1)');
        
    case 3                             %simetria horizontal 
        
        h = zeros(npy+4,npx+2);
        h(3:end,3:end) = 1;
        h1 = h;
        h0 = h;
        h(h==1)=1:sum(sum(h));
        
        h(3:end-2,1) = h(3:end-2,3);   %fronteira 4
        h(1,3:end)   = h(3,3:end);     %fronteira 1
        h(end,3:end) = h(end-2,3:end); %fronteira 3

        h1(:,1)   = tta(4);
        h1(1,:)   = tta(1);
        h1(end,:) = tta(3);

        h2 = fliplr(h(:,1:end-1));

        h12 = fliplr(h1(:,1:end-1));
        
        h02 = fliplr(h0(:,1:end-1));

        H0 = [h0 h02];
        H = [h h2];
        hh = [h1 h12];   
        
        s = sum(sum(H0));
        pos = (find(H0==1)');
        
    case 2                            % simetria vertical
        h = zeros(npy+2,npx+4);
        h(3:end,3:end) = 1;
        h1 = h;
        h0 = h;
        h(h==1)=1:sum(sum(h));
        
        h(3:end,1)   = h(3:end,3);     %fronteira 4
        h(1,3:end-2) = h(3,3:end-2);   %fronteira 1
        h(end,3:end) = h(end-2,3:end); %fronteira 2
        
        h1(:,1)   = tta(4);
        h1(1,:)   = tta(1);
        h1(:,end) = tta(2);

        h3 = flipud(h(1:end-1,:));

        h13 = flipud(h1(1:end-1,:));

        h03 = flipud(h0(1:end-1,:));

        H0 = [h0; h03];
        H = [h; h3];
        hh = [h1; h13]; 
        
        s = sum(sum(H0));
        pos = (find(H0==1)');
        
    case 1                       %sem simetria
        h = zeros(npy+4,npx+4);
        h(3:end-2,3:end-2) = 1;
        if tta(1)==0; h(2,2:end-2) = 1; end
        if tta(4)==0; h(3:end-1,end-1) = 1; end
        if tta(3)==0; h(end-1,3:end-1) = 1; end
        if tta(2)==0; h(3:end-1,2) = 1; end
        h1 = h;
        h0 = h;
        h(h==1)=1:sum(sum(h));
        
        if tta(2)~=0; h(3:end,1)   = h(3:end,3);     end; %fronteira 4
        if tta(1)~=0; h(1,3:end-1) = h(3,3:end-1);   end; %fronteira 1
        if tta(3)~=0; h(end,3:end) = h(end-2,3:end); end; %fronteira 2
        if tta(4)~=0; h(3:end,end) = h(3:end,end-2); end; %fronteira 3
        
        h1(:,1)   = tta(2);
        h1(1,:)   = tta(1);
        h1(:,end) = tta(4);
        h1(end,:) = tta(3);  
        
        H  = zeros(npx+6,npy+6);
        hh = zeros(npx+6,npy+6);
        hhh = zeros(npx+6,npy+6);
        H0 = h0;
        H(2:end-1,2:end-1)  = h;
        hh(2:end-1,2:end-1) = h1;
        
        hhh(3:end-2,3:end-2) = hh(3:end-2,3:end-2)~=0;
        
        pos = (find(hhh'));
        s = sum(sum(hhh));
end
%--------------------------------------------------------------------------
figure
spy(H0)
clear H0 h0 h02 h03 h04 h1 h12 h13 h14 hhh
%--------------------------------------------------------------------------
mat1 =  [0  0  1  0 0
        0  2 -8  2 0
        1 -8 20 -8 1
        0  2 -8  2 0
        0  0  1  0 0];
%--------------------------------------------------------------------------
res = zeros(s);
rrr = (qq*(dxy^4)/D);
resb = repmat(rrr,s,1);
%--------------------------------------------------------------------------
A = zeros(s);
%--------------------------------------------------------------------------
for k = 1:numel(pos)
    [I,J] = ind2sub(size(h),pos(k))
    matpos = [0        0          H(J-2,I) 0          0
              0        H(J-1,I-1) H(J-1,I) H(J-1,I+1) 0
              H(J,I-2) H(J,I-1)   H(J,I)   H(J,I+1)   H(J,I+2)
              0        H(J+1,I-1) H(J+1,I) H(J+1,I+1) 0
              0        0          H(J+2,I) 0          0];
          
    matval = [0         0           hh(J-2,I) 0           0
              0         hh(J-1,I-1) hh(J-1,I) hh(J-1,I+1) 0
              hh(J,I-2) hh(J,I-1)   hh(J,I)   hh(J,I+1)   hh(J,I+2)
              0         hh(J+1,I-1) hh(J+1,I) hh(J+1,I+1) 0
              0         0           hh(J+2,I) 0           0];
          pos1 = (find(matpos>0)');
    for kkk = 1:numel(pos1())
        kk = pos1(kkk);
        A(matpos(13),matpos(kk)) = A(matpos(13),matpos(kk))+mat1(kk)*matval(kk);
    end
end
disp('Resultados - RR')
RR = A\resb;
disp('Max(RR)')
R = max(abs(RR))
disp('Min(RR)')
R = min(abs(RR))
%--------------------------------------------------------------------------
chapa = H(2:end-1,2:end-1);
for k = 1:numel(RR)
    chapa(chapa==k) = RR(k);
end
toc
%--------------------------------------------------------------------------
tipo = 'estrutural2dv2_res';
save(strcat(nome,'_res_e','.mat'),'tipo','chapa');
disp('------------------')
disp('Resultado guardado')
%--------------------------------------------------------------------------

%% Mostrar resultados
%--------------------------------------------------------------------------
load(strcat(input('Nome do ficheiro de resultados -> ','s'),'.mat'))
if exist('tipo','var')~= 1; warning('Ficheiro inválido!'); break; end;
if ~strcmp(tipo,'estrutural2dv2_res'); warning('Ficheiro inválido!'); break; end
%--------------------------------------------------------------------------

figure
contour(chapa); 
colorbar;
axis equal
figure
surf(chapa)
colorbar;
shading interp
axis([0 size(chapa,2) 0 size(chapa,1) -0.00002 0.00001])
%--------------------------------------------------------------------------

%% agra.
% Agradecemos ao Mario por ter reparado o codigo, sem ele as condições de
% fronterira nao seriam tidas em conta
disp(' Agradecemos ao Mario por ter reparado o codigo, sem ele as condições de fronterira nao seriam tidas em conta')