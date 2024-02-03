clc; clear all; close all;

% UCITAVANJE SLIKE
fname=uigetfile('*.gif');% otvara ui za odabranu sliku
I=imread(fname); % cita odabranu sliku
I=imresize(I,[256 256]);% zbog lakše komputacije, smanjuje sliku na 256*256
figure,imshow(I);title('Originalna slika');drawnow; % prikaz slike
R=8; % blokovi kodomena (4x4, 8x8)
D=2*R; % blokovi domena duplo veci
step=8; % korak pri dijeljenju slike na blokove domena (4, 8)
n=R*R; % broj piksela u jednom bloku slike

% DFINISANJE MATRICA I PRORAcUN AFFINE TRANSFORMACIJA
tx = 0;            % translacija po x osi
ty = 0;            % translacija po y osi

T{1} = [1  0  0; 0  1  0; tx  ty  1]; 
T{2} = [-1  0  0; 0  1  0; tx  ty  1]; 
T{3} = [1  0  0; 0  -1  0; tx  ty  1]; 
T{4} = [-1  0  0; 0  -1  0; tx  ty  1]; 
T{5} = [0  1  0; 1  0  0; tx  ty  1]; 
T{6} = [0  -1  0; 1  0  0; tx  ty  1]; 
T{7} = [0  1  0; -1  0  0; tx  ty  1]; 
T{8} = [0  -1  0; -1  0  0; tx  ty  1]; 

for i=1:length(T)
    aff{i}=affine2d(T{i});
end

% PODJELA SLIKE NA BLOKOVE KODOMENE R

br1=1;
for i=1:R:256
    for j=1:R:256
    blok=I(i:i+R-1,j:j+R-1);
    Range{br1}=blok;
    br1=br1+1;
    end
end


% ODREDjIVANJE PRAGOVA Ts i Td
l=1;
for j=1:length(Range)  
    dctkoef=dct2(Range{j}); % proracun dct koeficijenata
    f10=abs(dctkoef(2,1));
    f01=abs(dctkoef(1,2));
    fovi{j}=[f10 f01]; %spremanje vrijednosti dct koef u niz
    fz=max(f10,f01); % odredjivanje vrijednosti Fmax za prvi prag
    fzzz(j)=fz;
end
[fzzz,r]=sort(fzzz);
c=round(length(Range)/3); % vrijednost na prvoj trecini odgovara pragu TS
Ts=fzzz(c);
l=1;

for b=c+1:length(Range)
    fk=abs(fovi{r(b)}(1)-fovi{r(b)}(2)); %proracun apsolutne vrijednosti za odredjivanje drugog praga
    fkkk(l)=fk; 
    l=l+1;
end    
fkkk=sort(fkkk);
Td=fkkk(round(length(fkkk)/2)); %vrijednost na polovini preostalih blokova odgovara pragu TD

sr=1;
dr=1;
hr=1;

% klasifikacija blokova R na klase S,D,H
Newrange=Range;
for qq=1:length(Range)  
    dctkoef=dct2(Range{qq});
    f10=abs(dctkoef(2,1));
    f01=abs(dctkoef(1,2));
    if f10<Ts && f01<Ts
        typeSr{sr}=Range{qq};
        sr=sr+1;
        Newrange{qq}(1:R,1:R)=0; % blokovi S obojeni u crno
    else if abs(f10-f01)<Td
            typeDr{dr}=Range{qq};
            dr=dr+1;
            Newrange{qq}(1:R,1:R)=500; % blokovi D obojeni u bijelo
        else
            typeHr{hr}=Range{qq};
            hr=hr+1;
            Newrange{qq}(1:R,1:R)=100; % blokovi H obojeni u sivo
        end
    end
end

% prikaz slike pomocu klasificiranih blokova
pr=1;
for qr=1:256/R
  for vr=1:256/R
     oldimgr{qr,vr}=Newrange{pr};
     pr=pr+1;
     end
end
Pr=cell2mat(oldimgr);
figure; imshow(Pr); title('Slika blokova S, D i H');drawnow;


% PODJELA SLIKE NA BLOKOVE DOMENE D

br2=1;
br3=1;
for k=1:step:256-D+1
    for h=1:step:256-D+1
            blok2=I(k:k+D-1,h:h+D-1);
            Domain{br2}=imresize(blok2,[R R]); % smanjenje velicine bloka D na velicinu R (poduzorkovanje)
        
        % PRIMJENA AFFINE TRANSFORMACIJA NA BLOKOVE DOMENA D
        for pp=1:length(T)
            Newdomain{br3}=imwarp(Domain{br2},aff{pp});
            aff_odg{br3}=aff{pp}; %spremanje affinih transf.
            br3=br3+1;
        end
        br2=br2+1;
    end
end

% RASPODJELA BLOKOVA DOMENE U 3 KLASE: S, D I H
s=1;
h=1;
d=1;
for q=1:length(Newdomain)  
    dctkoef=dct2(Newdomain{q});
    f10=abs(dctkoef(2,1));
    f01=abs(dctkoef(1,2));
    f{q}=[f10 f01];
    % ispitivanje uslova 
    if f10<Ts && f01<Ts
        typeS{s}=Newdomain{q};
        affS{s}=aff_odg{q};
        koorS(s)=q;
        s=s+1;
    else if abs(f10-f01)<Td
            typeD{d}=Newdomain{q};
            affD{d}=aff_odg{q};
            koorD(d)=q;
            d=d+1;
        else
            typeH{h}=Newdomain{q};
            affH{h}=aff_odg{q};
            koorH(h)=q;
            h=h+1;
        end
    end
end

% PRORACUN KOEFICIJENATA BRIGHTNESS I CONTRAST (si, oi)
% Nakon primjene affinih transformacija i podjele blokova u klase
% izracunati su koeficijenti kontrasta si i osvjetljenja oi. 
% Dakle, u ovom slucaju koeficijenti su izracunati za one blokove R i
% transformisane blokove D koji pripadaju istoj klasi (S, D ili H).
% Na osnovu toga odredjeni su si i oi koef. za blokove klasa S, D i H.

% KODIRANJE FIC-DCT: POREDjENJE BLOKOVA IZ ODGOVARAJUCIH KLASA, TE
% PRONALAZAK NAJSLICNIJEG DOMENSKOG BLOKA SA OPTIMALNIM KOEF. si i oi i
% afinom transform., tako da zadovoljava min vrijednost MSE.

brmse=1;
brmse1=0;
brmse2=0;
brmse3=0;
tic;
for a=1:length(Range)
    % Poredjenje blokova R sa blokovima D iz odgovarajuce klase
    % te traženje minimalne MSE na osnovu primjenjenih tranformacija i si
    % oi koef.
    if fovi{a}(1)<Ts && fovi{a}(2)<Ts
        for sp=1:length(typeS)
            Ra=double(Range{a});
            Db=double(typeS{sp});
            % sume u formulama za si i oi su izracunate pojedinacno radi
            % jednostavnosti
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula za si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula za oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err1(sp)=mse;
            koefS{sp}=[s o];
            brmse1=brmse1+1;
        end
        [M,In]=min(err1);
        aff_rekons{a}=affS{In}; % spremanje odgovarajuce affine transf.
        koordinate(a)=ceil(koorS(In)/8); % pozicija bloka D prije affinih transf.
        s_o{a}=koefS{In}; % spremanje optimalnih blokova si oi
    else if abs(fovi{a}(1)-fovi{a}(2))<Td
        for dp=1:length(typeD)
            Ra=double(Range{a});
            Db=double(typeD{dp});
            % sume u formulama za si i oi su izracunate pojedinacno radi
            % jednostavnosti
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula za si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula za oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err2(dp)=mse;
            koefD{dp}=[s o];
            brmse2=brmse2+1;
        end
        [M,In]=min(err2);
        aff_rekons{a}=affD{In}; % spremanje odgovarajuce affine transf.
        koordinate(a)=ceil(koorD(In)/8); % pozicija bloka D prije affinih transf.
        s_o{a}=koefD{In}; % spremanje optimalnih blokova si oi
        else
        for hp=1:length(typeH)
            Ra=double(Range{a});
            Db=double(typeH{hp});
            % sume u formulama za si i oi su izracunate pojedinacno radi
            % jednostavnosti
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula za si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula za oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err3(hp)=mse;
            koefH{hp}=[s o];
            brmse3=brmse3+1;
        end
        [M,In]=min(err3);
        aff_rekons{a}=affH{In}; % spremanje odgovarajuce affine transf.
        koordinate(a)=ceil(koorH(In)/8); % pozicija bloka D prije affinih transf.
        s_o{a}=koefH{In}; % spremanje optimalnih blokova si oi
        end
    end
            
end
brmse=brmse1+brmse2+brmse3;
Vrijeme_kodiranja=toc
Broj_poredjenja=brmse

% % DEKODIRANJE

fname=uigetfile('*.gif');% otvara ui za inicijalnu sliku
P=imread(fname); % cita odabranu sliku
RekonsI=P; % inicijalna slika


for it=1:15
    broj2=1;
    for kk=1:step:256-D+1
        for hh=1:step:256-D+1
                blokk=RekonsI(kk:kk+D-1,hh:hh+D-1);
                Domain_dekod{broj2}=imresize(blokk, [R R]);
                broj2=broj2+1;
        end
    end
    for f=1:length(Range)
    RekonsIm=imwarp(Domain_dekod{koordinate(f)},aff_rekons{f}); % primjena affine transformacije
    mat1=ones(size(RekonsIm)).*s_o{f}(2); 
    % optimizacija kontrasta i osvjetljenja na osnovu si i oi
    ReI{f}=imresize(uint8(double(RekonsIm).*s_o{f}(1)+mat1), [R R]);
    end
    % Rekonstrukcija slike
    p=1;
    for q=1:length(I)/R
     for v=1:length(I)/R
         oldimg{q,v}=ReI{p};
         p=p+1;
     end
    end
    it
    RekonsI=cell2mat(oldimg); % slika u sljedecoj iteraciji
end



% PRIKAZ REKONSTRUISANE SLIKE

figure, imshow(RekonsI); title('Rekonstruisana slika');drawnow;
[peaksnr,~] = psnr(RekonsI,I);
PSNR=peaksnr