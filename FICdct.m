clc; clear all; close all;

% IMAGE LOADING
fname=uigetfile('*.gif');% opens ui for the chosen image
I=imread(fname); % reads chosen image
I=imresize(I,[256 256]);% resizing image to 256*256 for better computation
figure,imshow(I);title('Originalna slika');drawnow; % image display
R=8; % range blocks (4x4, 8x8)
D=2*R; % domain blocks (2x range)
step=8; % step when dividing the image into domain blocks (4, 8)
n=R*R; % number of pixels in one image block

% DEFINITION OF MATRICES AND CALCULATION OF AFFINE TRANSFORMATIONS
tx = 0;            % x-axis translation
ty = 0;            % y-axis translation
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

% DIVIDING IMAGE INTO RANGE BLOCKS R
num1=1;
for i=1:R:256
    for j=1:R:256
    block=I(i:i+R-1,j:j+R-1);
    Range{num1}=block;
    num1=num1+1;
    end
end


% THRESHOLDS DETERMINATION Ts & Td
l=1;
for j=1:length(Range)  
    dctcoeff=dct2(Range{j}); % DCT coefficients calculation
    f10=abs(dctcoeff(2,1));
    f01=abs(dctcoeff(1,2));
    fovi{j}=[f10 f01]; % saving DCT values to an array 
    fz=max(f10,f01); % value Fmax for the first threshold
    fzzz(j)=fz;
end
[fzzz,r]=sort(fzzz);
% the value in the first third corresponds to the TS threshold
c=round(length(Range)/3); 
Ts=fzzz(c);
l=1;

for b=c+1:length(Range)
    % calculation of the absolute value for determining the second threshold
    fk=abs(fovi{r(b)}(1)-fovi{r(b)}(2));
    fkkk(l)=fk; 
    l=l+1;
end    
fkkk=sort(fkkk);
% the value on half of the remaining blocks corresponds to the TD threshold
Td=fkkk(round(length(fkkk)/2));

sr=1;
dr=1;
hr=1;

% Range blocks classification into classes S,D,H
Newrange=Range;
for qq=1:length(Range)  
    dctcoeff=dct2(Range{qq});
    f10=abs(dctcoeff(2,1));
    f01=abs(dctcoeff(1,2));
    if f10<Ts && f01<Ts
        typeSr{sr}=Range{qq};
        sr=sr+1;
        Newrange{qq}(1:R,1:R)=0; % S blocks - black
    else if abs(f10-f01)<Td
            typeDr{dr}=Range{qq};
            dr=dr+1;
            Newrange{qq}(1:R,1:R)=500; % D blocks - white
        else
            typeHr{hr}=Range{qq};
            hr=hr+1;
            Newrange{qq}(1:R,1:R)=100; % H blocks - gray
        end
    end
end

% image display using classified blocks
pr=1;
for qr=1:256/R
  for vr=1:256/R
     oldimgr{qr,vr}=Newrange{pr};
     pr=pr+1;
     end
end
Pr=cell2mat(oldimgr);
figure; imshow(Pr); title('Slika blockova S, D i H');drawnow;


% DIVIDING IMAGE INTO DOMAIN BLOCKS D
num2=1;
num3=1;
for k=1:step:256-D+1
    for h=1:step:256-D+1
            block2=I(k:k+D-1,h:h+D-1);
        % reducing the block size D to the size R (subsampling)
            Domain{num2}=imresize(block2,[R R]);
        % Applying AFFINE Transformations on domain blocks D
        for pp=1:length(T)
            Newdomain{num3}=imwarp(Domain{num2},aff{pp});
            aff_corr{num3}=aff{pp}; % Saving Affine Transformations
            num3=num3+1;
        end
        num2=num2+1;
    end
end

% DISTRIBUTION OF DOMAIN BLOCKS INTO 3 CLASSES: S, D & H
s=1;
h=1;
d=1;
for q=1:length(Newdomain)  
    dctcoeff=dct2(Newdomain{q});
    f10=abs(dctcoeff(2,1));
    f01=abs(dctcoeff(1,2));
    f{q}=[f10 f01];
    % testing of conditions
    if f10<Ts && f01<Ts
        typeS{s}=Newdomain{q};
        affS{s}=aff_corr{q};
        coorS(s)=q;
        s=s+1;
    else if abs(f10-f01)<Td
            typeD{d}=Newdomain{q};
            affD{d}=aff_corr{q};
            coorD(d)=q;
            d=d+1;
        else
            typeH{h}=Newdomain{q};
            affH{h}=aff_corr{q};
            coorH(h)=q;
            h=h+1;
        end
    end
end

% CALCULATION OF BRIGHTNESS AND CONTRAST COEFFICIENTS (oi, si)
% After applying affine transformations and dividing blocks into classes
% the coefficients of contrast si and brightness oi were calculated.
% So, in this case, the coefficients were calculated for those blocks R and
% transformed blocks D belonging to the same class (S, D or H).
% Based on this, si and oi coefficients were determined for blocks of classes S, D and H.

% FIC-DCT CODING: COMPARISON OF BLOCKS FROM CORRESPONDING CLASSES, AND
% FINDING THE MOST SIMILAR DOMAIN BLOCK WITH THE OPTIMUM COEF. si & oi &
% affine transform., so that it meets the min value of MSE.

nummse=1;
nummse1=0;
nummse2=0;
nummse3=0;
tic;
for a=1:length(Range)
     % Comparison of blocks R with blocks D from the corresponding class
     % and searching for the minimum MSE based on the applied transformations and si
     % & oi coef.
    if fovi{a}(1)<Ts && fovi{a}(2)<Ts
        for sp=1:length(typeS)
            Ra=double(Range{a});
            Db=double(typeS{sp});
            % sums in the formulas for si & oi are calculated individually for simplicity
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula for si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula for oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err1(sp)=mse;
            coeffS{sp}=[s o];
            nummse1=nummse1+1;
        end
        [M,In]=min(err1);
        aff_recons{a}=affS{In}; % saving corresponding affine transf.
        coordinates(a)=ceil(coorS(In)/8); % D block position before affine transf.
        s_o{a}=coeffS{In}; % saving optimal blocks si & oi
    else if abs(fovi{a}(1)-fovi{a}(2))<Td
        for dp=1:length(typeD)
            Ra=double(Range{a});
            Db=double(typeD{dp});
            % sums in the formulas for si & oi are calculated individually for simplicity
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula for si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula for oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err2(dp)=mse;
            coeffD{dp}=[s o];
            nummse2=nummse2+1;
        end
        [M,In]=min(err2);
        aff_recons{a}=affD{In}; % saving corresponding affine transf.
        coordinates(a)=ceil(coorD(In)/8); % D block position before affine transf.
        s_o{a}=coeffD{In}; % saving optimal blocks si & oi
        else
        for hp=1:length(typeH)
            Ra=double(Range{a});
            Db=double(typeH{hp});
            % sums in the formulas for si & oi are calculated individually for simplicity
            sum1=sum(sum(Ra.*Db));
            sum2=sum(sum(Db))*sum(sum(Ra));
            sum3=sum(sum(Db.^2));
            sum4=(sum(sum(Db)))^2;
            % formula for si:
            s=((n)*(sum1)-sum2)/((n)*sum3-sum4);
            % formula for oi:
            o=(sum(sum(Ra))-s*sum(sum(Db)))/((n));
            mat=ones(size(Ra)).*o;    
            mse=immse(Range{a},uint8(Db.*s+mat));
            err3(hp)=mse;
            coeffH{hp}=[s o];
            nummse3=nummse3+1;
        end
        [M,In]=min(err3);
        aff_recons{a}=affH{In}; % saving corresponding affine transf.
        coordinates(a)=ceil(coorH(In)/8); % D block position before affine transf.
        s_o{a}=coeffH{In}; % saving optimal blocks si & oi
        end
    end
            
end
nummse=nummse1+nummse2+nummse3;
Coding_time=toc
Comparisons=nummse

% DECODING

fname=uigetfile('*.gif');% opens ui for initial image
P=imread(fname); % reads chosen image
reconsI=P; % initial image


for it=1:15
    numb2=1;
    for kk=1:step:256-D+1
        for hh=1:step:256-D+1
                blockk=reconsI(kk:kk+D-1,hh:hh+D-1);
                Domain_decod{numb2}=imresize(blockk, [R R]);
                numb2=numb2+1;
        end
    end
    for f=1:length(Range)
    % applying affine transformation
    reconsIm=imwarp(Domain_decod{coordinates(f)},aff_recons{f});
    mat1=ones(size(reconsIm)).*s_o{f}(2); 
    % optimization of contrast and brightness based on si and oi
    ReI{f}=imresize(uint8(double(reconsIm).*s_o{f}(1)+mat1), [R R]);
    end
    % Image reconstruction
    p=1;
    for q=1:length(I)/R
     for v=1:length(I)/R
         oldimg{q,v}=ReI{p};
         p=p+1;
     end
    end
    it
    reconsI=cell2mat(oldimg); % image in the next iteration
end



% DISPLAY OF RECONSTRUCTED IMAGE
figure, imshow(reconsI); title('reconstruisana slika');drawnow;
[peaksnr,~] = psnr(reconsI,I);
PSNR=peaksnr