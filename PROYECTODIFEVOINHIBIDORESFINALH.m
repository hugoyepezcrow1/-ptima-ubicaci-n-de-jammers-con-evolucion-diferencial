% Determina la cobertura de inhibidores de se�al instalados en la mejor
% posici�n para tener la mayor covertura, el software fue desarrollado en
% Matlab 2018 b, simulado en una computadora con procesador i7 y 12 Gb de
% ramde

clc
a=2;     %Ingrese el largo del centro de rehabilitaci�n social
b=1;   %Ingrese el ancho del fonde de rehabilitaci�n social

Nb=20;   % numero de inhibidores
F= 0.4+0.6*rand;   %factor escalar
Cr=0.7+0.3*rand;  % tasa de cruce
Fp=0.95;   %factor de protecci�n perimetro
% tama�os de radio de cobertura del inhibidor:
mab=mean([a b]);
aa=0.15*mab+0.3*mab*rand(1,Nb);
bb=aa; %Variable auxiliar para el radio
m2=min([aa bb]/2); 
AA=(aa/2.*bb/2)*pi; % areas de los inhibidores de se�al

penalty=1e-6*a*b;
nac=0.8; % coheficiente de �rea negativa



N=100; % Tama�o de la poblacion
ng=10000; % numero de generaciones
pmpe=0.5; % probabilidad de plazas de intercambio de mutaci�n
pmbj=0.5; % Saltos grandes
pmsj=0.5; % Saltos peque�os
pmrr=0.5; % traslacion
pmvi=0.5; % visible/invisible
pmne=0.5; % se mueve al limite cercano

%limites de la grafica

figure;
ha1=axes;
ha1=subplot(2,1,1);
plot([0 a a 0 0], [0 0 b b 0],'b-');
xlim([-0.1*a 1.1*a]);
ylim([-0.1*b 1.1*b]);
set(ha1,'NextPlot','add');
ht=title(ha1,'start');
ha2=subplot(2,1,2);
drawnow;


set_cl; % pone los clores de los inhibidores de diferentes colores


% Iniciando el proceso de cubrimiento 
tic % Para el control de tiempo
 
% inicio de poblaci�n randomicamente:
G=zeros(N,4*Nb);
Gch=zeros(N,4*Nb); % hijo
for Nc=1:N % para cada individuo
    G1=zeros(4,Nb); % 1 individuao
    % G1(1,i)=1 si i-inhibidor es visible
    % G1(2,i)=1 si i-inhibidor es movido
    % G1(3,i) - coordenada x del centro del inhibidor
    % G1(4,i) - coordenada y del centro del inhibidor

    G1(1,:)=double(rand(1,Nb)<0.2);
    G1(2,:)=double(rand(1,Nb)<0.5);
    G1(3,:)=m2+(a-m2)*rand(1,Nb);
    G1(4,:)=m2+(b-m2)*rand(1,Nb);
    G(Nc,:)=(G1(:))'; % (G1(:))' convierte matriz a vector fila
end



hi=imagesc(G,'parent',ha2);
ylabel('n�mero de individuos');
title('genes');
colorbar;
drawnow;

Gpr1=zeros(4,Nb);
Gpr2=zeros(4,Nb); % 2 padres
Gch1=zeros(4,Nb);
Gch2=zeros(4,Nb); % hijos

delta=1e-6;

 % Mutaci�n: Esquema de mutaci�n adaptativo
 if delta<1e-6
    % intercambio de genes
    for Nc=1:N % para cada individuo
        if rand>pmpe
            G1(:)=(G(Nc,:))';
            ir1=ceil(Nb*rand);
            ir2=ceil(Nb*rand);
            ir3=ceil(Nb*rand);
            ir4=ceil(Nb*rand);
            G1(3:4,ir4)=G1(3:4,ir1)+F*(G1(3:4,ir3)-G1(3:4,ir2));
            G(Nc,:)=(G1(:))';
            
        else
             if rand<pmpe
            G1(:)=(G(Nc,:))';
            ir1=ceil(Nb*rand);
            ir2=ceil(Nb*rand);
            ir3=ceil(Nb*rand);
            ir4=ceil(Nb*rand);
            ir5=ceil(Nb*rand);
            ir6=ceil(Nb*rand);
            G1(3:4,ir6)=G1(3:4,ir1)+F*(G1(3:4,ir2)-G1(3:4,ir3))+F*(G1(3:4,ir5)-G1(3:4,ir4));
            G(Nc,:)=(G1(:))';
             end
             
        end
    end
    
 else
        if delta>1e-6
    % intercambio de genes
    for Nc=1:N % para cada individuo
        if rand>pmpe
            G1(:)=(G(Nc,:))';
            ir1=ceil(Nb*rand);
            ir2=ceil(Nb*rand);
            ir3=ceil(Nb*rand);
            G1(3:4,ir3)=Gb+F*(G1(3:4,ir1)-G1(3:4,ir2));
            G(Nc,:)=(G1(:))';
            
        else
             if rand<pmpe
            G1(:)=(G(Nc,:))';
            ir1=ceil(Nb*rand);
            ir2=ceil(Nb*rand);
            ir3=ceil(Nb*rand);
            ir4=ceil(Nb*rand);
            ir5=ceil(Nb*rand);
            G1(3:4,ir5)=Gb+F*(G1(3:4,ir1)-G1(3:4,ir2))+F*(G1(3:4,ir3)-G1(3:4,ir4));
            G(Nc,:)=(G1(:))';
        end
             
        end
    end
        end
 end
       
for ngc=1:ng % conteo de generaciones
    % encuentra la funci�n objetivo:
    fitnesses=zeros(N,1);
    for Nc=1:N % para cada individuo
        G1(:)=(G(Nc,:))';
        vis=G1(1,:);
        ind=find(vis);
        L=length(ind);
        if L>0
            % solo inhibidores visibles:
            rot=G1(2,ind);
            x=G1(3,ind);
            y=G1(4,ind);
            if L==1
                aaa=aa(ind);
                bbb=bb(ind);
                if rot
                    tmp=aaa;
                    aaa=bbb;
                    bbb=tmp;
                end
                A0=AA(ind); % area del inhibidor
                x1=max([x-aaa/2  0]);
                y1=max([y-bbb/2  0]);
                x2=min([x+aaa/2  a]);
                y2=min([y+bbb/2  b]);
                % x1 - x2,  y1 - y2 coordenadas dentro del �rea de
                % cobertura
                if (x1>=x2)||(y1>=y2)||(x1>=a*Fp)||(y1>=b*Fp)
                    A=0; % inhibidor dentro del area principal
                else
                    A=((x2-x1)/2*(y2-y1)/2)*pi; % calculo del area del inhibidor
                end
                %if A<A0 % si no esta totalmente adentro del area principal
               if (aaa/2<=x*Fp)&&(x<=a-aaa/2)&&(bbb/2<=y*Fp)&&(y<=b-bbb/2) % si esta totalmente adentro
                    fitness=A;
               else
                    nac=nac+(0.0081/log(2*nac));
                  
                    fitness=A-nac*(A0-A)-penalty;
                end
               
                    
            else
                fitness=0;
                ispen=false; % es 1 si hay penalidad
                
                % chequeo de cruce con el area de cobertura principal:
                % anade inhibidores y sustrae a los que estan afuera:
                for n=1:L % por cada inhibidor
                    ind1=ind(n);
                    aaa=aa(ind1);
                    bbb=bb(ind1);
                    if rot(n)
                        tmp=aaa;
                        aaa=bbb;
                        bbb=tmp;
                    end
                    A0=AA(ind1); %  area inhibidor
                    x1=max([x(n)-aaa/2  0]);
                    y1=max([y(n)-bbb/2  0]);
                    x2=min([x(n)+aaa/2  a]);
                    y2=min([y(n)+bbb/2  b]);
                    % x1 - x2,  y1 - y2 esta el inhibidor dentro del area
                    % de cobertura
                    if (x1>=x2)||(y1>=y2)||(x1>=a*Fp)||(y1>=b*Fp)
                        A=0; % inhibidor dentro del area principal
                    else
                        A=((x2-x1)/2*(y2-y1)/2)*pi; % inhibidor dentro del area principal
                    end
                    
                     %para probar
                    if A<A0 % if not fully inside main box
                       fitness=fitness + A-nac*(A0-A);
                        ispen=true; % penality
                    else
                        fitness=fitness + A;
                   end
                    
                    if (aaa/2<=x(n)*Fp)&&(x(n)<=a*Fp-aaa/2)&&(bbb/2<=y(n)*Fp)&&(y(n)<=b*Fp-bbb/2) % si esta totalmente adentro
                        fitness=fitness + A;
                    else
                        fitness=fitness + A-nac*(A0-A);
                        ispen=true; % penalidad
                    end
                    
                end
                
                % por cada par de inhibidores:
                for n1=1:L-1
                    ind1=ind(n1);
                    aaa1=aa(ind1);
                    bbb1=bb(ind1);
                    if rot(n1)
                        tmp=aaa1;
                        aaa1=bbb1;
                        bbb1=tmp;
                    end
                    A1=AA(ind1);
                    x1=x(n1);
                    y1=y(n1); % posicion del inhibidor
                    for n2=n1+1:L
                        ind2=ind(n2);
                        aaa2=aa(ind2);
                        bbb2=bb(ind2);
                        if rot(n2)
                            tmp=aaa2;
                            aaa2=bbb2;
                            bbb2=tmp;
                        end
                        A2=AA(ind2);
                        x2=x(n2);
                        y2=y(n2); % posicion de segundo inhibidor
                        dx=abs(x1-x2);
                        dy=abs(y1-y2); % distancias
                        a12=(aaa1/2+aaa2/2);
                        b12=(bbb1/2+bbb2/2);
                        if (dx<a12)&&(dy<b12) % si cruza
                            ispen=true;
                            Ac=((a12-dx)/2*(b12-dy)/2)*pi; % area de cruce
                            fitness=fitness-Ac-Ac; % porque area n1 y n2 fueron totalmente a�adidas
                            fitness=fitness-2*nac*Ac;
                        end

                    end
                end
                
                if ispen
                    fitness=fitness-penalty;
                end
        
            end
        else
            fitness=0;
        end
        fitnesses(Nc)=fitness;
    end
    
    [fb bi]=max(fitnesses); % el mejor
    
   %para acelerar la tasa de convergencia del algoritmo
   
         delta = abs(fitness/abs(mean(fitnesses))-1);
        
        
    
  % grafica del mejor:
    G1(:)=(G(bi,:))';
    Gb=G(bi,:); % best
    if mod(ngc,10)==0
        cla(ha1);
        Atmp=0;
        for Nbc=1:Nb
            vis1=G1(1,Nbc);
            if vis1
                rot1=G1(2,Nbc);
                aaa=aa(Nbc);
                bbb=bb(Nbc);
                
                if rot1
                    tmp=aaa;
                    aaa=bbb;
                    bbb=tmp;
                end
                x=G1(3,Nbc);
                y=G1(4,Nbc);
                
                 theta = linspace(0, 2*pi);
        xcentro1 = G1(3,Nbc)-aaa/2*cos(theta);
        xcentro2 = G1(3,Nbc)+aaa/2*cos(theta);
        
        ycentro1 = G1(4,Nbc)-bbb/2*sin(theta);
        ycentro2 = G1(4,Nbc)+bbb/2*sin(theta);
     
        
           plot([xcentro1  xcentro2  xcentro2  xcentro1  xcentro1],...
                 [ycentro1  ycentro1  ycentro2  ycentro2  ycentro1],...      
            '-','color',cl(Nbc,:),...
                     'parent',ha1);
                 hold on;
                Atmp=Atmp+aaa*bbb;
                
   
       
        %%%%%%%%%
                 
                
                
            end
        end
        plot([0 a a 0 0], [0 0 b b 0],'b-','parent',ha1);
        xlim(ha1,[-0.1*a 1.1*a]);
        ylim(ha1,[-0.1*b 1.1*b]);
        
        
        %%%%%%%%%
        
        
        set(hi,'Cdata',G);
        
        nvb=length(find(G1(1,:))); % numero de inhibidores visibles
        pcubierto=num2str(fb)/a*b;
        
        set(ht,'string',[' generation: ' num2str(ngc)  ', inhibidores: ' num2str(nvb) ', area: ' num2str(fb)]);
        
        drawnow;
    end
     
    
    % preparacion para el cruce, seleccion:
    fmn=min(fitnesses);
    fst=std(fitnesses);
    if fst<1e-7
        fst=1e-7;
    end
    fmn1=fmn-0.01*fst; % bajo entonces m�nimo
    P=fitnesses-fmn1; % valores positivos
    p=P/sum(P); % probabilidad
    ii=roulette_wheel_indexes(N,p); %uso del algoritmo de ruleta de casino para generar el cruce
    Gp=G(ii,:); % padres
    
    % cruzamiento:
    for n=1:2:N
        pr1=Gp(n,:);
        pr2=Gp(n+1,:); % 2 padres
        % en forma de matriz:
        Gpr1(:)=pr1'; 
        Gpr2(:)=pr2';
        
        for Nbc=1:Nb
            
            % visibilidad
            if rand<Cr
                Gch1(1,Nbc)=Gpr1(1,Nbc);
            else
                Gch1(1,Nbc)=Gpr2(1,Nbc);
            end
            if rand<Cr
                Gch2(1,Nbc)=Gpr1(1,Nbc);
            else
                Gch2(1,Nbc)=Gpr2(1,Nbc);
            end
            
            % movimiento:
            if rand<Cr
                Gch1(2,Nbc)=Gpr1(2,Nbc);
            else
                Gch1(2,Nbc)=Gpr2(2,Nbc);
            end
            if rand<Cr
                Gch2(2,Nbc)=Gpr1(2,Nbc);
            else
                Gch2(2,Nbc)=Gpr2(2,Nbc);
            end
            
            % posicion:
            % hijo 1:
           
            i3=1+ceil(2*rand);
            switch i3
                case 1 % obtiene la posicion media
                    Gch1(3,Nbc)=(Gpr1(3,Nbc)+Gpr2(3,Nbc))/2;
                    Gch1(4,Nbc)=(Gpr1(4,Nbc)+Gpr2(4,Nbc))/2;
                case 2 %obtiene la posicion del padre 1
                    Gch1(3,Nbc)=Gpr1(3,Nbc);
                    Gch1(4,Nbc)=Gpr1(4,Nbc);
                case 3 %obtiene la posicion del padre 2
                    Gch1(3,Nbc)=Gpr2(3,Nbc);
                    Gch1(4,Nbc)=Gpr2(4,Nbc);
            end
            
            % hijo 2:
          
            i3=1+ceil(2*rand);
            switch i3
                case 1 % obtiene la posicion media
                    Gch2(3,Nbc)=(Gpr1(3,Nbc)+Gpr2(3,Nbc))/2;
                    Gch2(4,Nbc)=(Gpr1(4,Nbc)+Gpr2(4,Nbc))/2;
                case 2 %obtiene la posicion del padre 1
                    Gch2(3,Nbc)=Gpr1(3,Nbc);
                    Gch2(4,Nbc)=Gpr1(4,Nbc);
                case 3 %obtiene la posicion del padre 2
                    Gch2(3,Nbc)=Gpr2(3,Nbc);
                    Gch2(4,Nbc)=Gpr2(4,Nbc);
            end
            
            
        end
        ch1=(Gch1(:))';
        ch2=(Gch2(:))';
        Gch(n,:)=ch1;
        Gch(n+1,:)=ch2;
        
        
    end
    G=Gch; % nuevo hijo
    
   
    % salto grande:
    for Nc=1:N % para cada individuo
        if rand<pmbj
            G1(:)=(G(Nc,:))';
            ir=ceil(Nb*rand);
            G1(3:4,ir)=G1(3:4,ir)+[0.05*a*Fp*randn;
                                   0.05*b*Fp*randn];
            G(Nc,:)=(G1(:))';
        end
    end
    
    % pequeno salto gauss:
    for Nc=1:N % para cada individuo
        if rand<pmsj
            G1(:)=(G(Nc,:))';
            ir=ceil(Nb*rand);
            G1(3:4,ir)=G1(3:4,ir)+[0.005*a*Fp*randn;
                                   0.005*b*Fp*randn];
            G(Nc,:)=(G1(:))';
        end
    end
    
    % movimiento aleatorio:
    for Nc=1:N % para cada individuo
        if rand<pmrr
            G1(:)=(G(Nc,:))';
            ir=ceil(Nb*rand);
            G1(2,ir)=double(rand<Cr);
            G(Nc,:)=(G1(:))';
        end
    end
    
    % aleatorio visible o no:
    for Nc=1:N % para cada individuo
        if rand<pmvi
            G1(:)=(G(Nc,:))';
            ir=ceil(Nb*rand);
            G1(1,ir)=double(rand<Cr);
            G(Nc,:)=(G1(:))';
        end
    end
    
    % se mueve al mas cercano limite:
    for Nc=1:N % para cada individuo
        if rand<pmne
            G1(:)=(G(Nc,:))';
            ir=ceil(Nb*rand); % inhibidor aleatorio pequeno
            rv=find((G1(1,:))&((1:Nb)~=Nc)); % encuentra el resto visible
            if rand<Cr
                % al borde vertical
                eax=[G1(3,rv)-aa(rv)/2  G1(3,rv)+aa(rv)/2  0  Fp*a]; % borde xs
                deax=[(G1(3,ir)-aa(ir)/2) - eax  (G1(3,ir)+aa(ir)/2) - eax]; % distancies
                [dmn indm]=min(abs(deax));
                G1(3,ir)=G1(3,ir)-deax(indm);
            else
                % al borde horizontal
                eay=[G1(4,rv)-bb(rv)/2  G1(4,rv)+bb(rv)/2  0  Fp*b]; % borde ys
                deay=[(G1(4,ir)-bb(ir)/2) - eay  (G1(4,ir)+bb(ir)/2) - eay]; % distancias
                
                [dmn indm]=min(abs(deay));
                G1(4,ir)=G1(4,ir)-deay(indm);
            end
        end
    end
    
  
    
    toc % Tiempo transcurrido
    
    % elite:
    G(1,:)=Gb;
   
        
 
end

