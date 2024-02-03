clc; clear all; close all;

% UCITAVANJE SLIKE
fname=uigetfile('*.gif');% otvara ui za odabranu sliku
I=imread(fname); % ?ita odabranu sliku
I=imresize(I,[256 256]);% zbog lakše komputacije, smanjuje sliku na 256*256
figure,imshow(I);title('Originalna slika');drawnow; % prikaz slike
vel=size(I); % veli?ina slike
R=8; % blokovi kodomena (4x4, 8x8)
D=2*R; % blokovi domena duplo veci
step=8; % korak pri dijeljenju slike na blokove domena (4, 8)
n=R*R; % broj piksela u jednom bloku slike

% DFINISANJE MATRICA I PRORACUN AFFINE TRANSFORMACIJA
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

% PODJELA SLIKE NA BLOKOVE DOMENE D

br2=1;
br3=1;
for k=1:step:256-D+1
    for h=1:step:256-D+1
            blok2=I(k:k+D-1,h:h+D-1);
            %Domain{br2}=skaliranje(blok2, 2, 2); % skaliranje bloka D na velicinu R
            Domain{br2}=imresize(blok2, [R R]);
        % PRIMJENA AFFINE TRANSFORMACIJA NA BLOKOVE DOMENA D
        % Na blokove domena nakon decimiranja izvršeno je 8 geometrijskih
        % afinih transformacija, te je broj blokova domena 8*D_blokova.
        for pp=1:length(T)
            Newdomain{br3}=imwarp(Domain{br2},aff{pp});
            af_odg{br3}=aff{pp}; %spremanje affinih transf.
            br3=br3+1;
        end
            br2=br2+1;
    end
end

% PRORACUN KOEFICIJENATA BRIGHTNESS I CONTRAST (si, oi)
% Nakon primjene affinih transformacija izracunati su koeficijenti
% kontrasta si i osvjetljenja oi. Dakle, za svaki blok R i skup blokova
% domena nakon afinih transformacija, odredio se par koeficijenata.

% KODIRANJE
tic;
brmse=1;
br1=0;
for a=1:length(Range)
    for b=1:length(Newdomain)
        Ra=double(Range{a});
        Db=double(Newdomain{b});
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
        Dom{b}=uint8(Db.*s+mat);
        mse=immse(Range{a},Dom{b}); % proracun mse blokova R i transformisanih
        % blokova D, gdje su uracunati i parametri kontrasta i osvjetljenja
        err(b)=mse; % spremanje MSE
        koef{b}=[s o]; %spremanje para koeficijenata [si oi] u niz
        brmse=brmse+1; % broj poredjenja blokova
    end
    [M,In]=min(err); % pronalazenje minimalne MSE, te odgovarajuceg bloka D
    % na osnovu indexa In pronasao se odgovarajuci transformisan blok D,
    % koji je najslicniji bloku R.
    aff_rekons{a}=af_odg{In}; % spremanje odgovarajuce affine transf.
    koordinate(a)=ceil(In/8); % pozicija bloka D prije affinih transf.
    s_o{a}=koef{In}; % spremanje optimalnih blokova si oi
end
% 
Vrijeme_kodiranja=toc
Broj_poredjenja=brmse-1

% DEKODIRANJE
% Na osnovu transformisanih blokova D pronasle su se afine transformacije
% i koeficijentima si i oi koji daju najmanju MSE za odg. blokove R i D.
% Spremljeni optimalni koef. i afine transformacije primijenjeni su na 
% inicijalnu sliku kako bi se formirala rekonstruisana slika.

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
    if (it==1 || it==2 || it==3 || it==5)
        figure, imshow(RekonsI);
    end
end



% PRIKAZ REKONSTRUISANE SLIKE
figure, imshow(RekonsI); title('Rekonstruisana slika');drawnow;
[peaksnr,~] = psnr(RekonsI,I);
PSNR=peaksnr