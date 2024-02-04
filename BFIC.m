clc; clear all; close all;

% IMAGE LOADING
fname=uigetfile('*.gif');% opens ui for the chosen image
I=imread(fname); % % reads chosen image
I=imresize(I,[256 256]);% resizing image to 256*256 for better computation
figure,imshow(I);title('Originalna slika');drawnow; % image display
R=16; % range blocks (4x4, 8x8)
D=2*R; % domain blocks (2x range)
step=8; % step when dividing the image into domain blocks (4, 8)
n=R*R; %  number of pixels in one image block

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
n1=1;
for i=1:R:256
    for j=1:R:256
        block=I(i:i+R-1,j:j+R-1);
        Range{n1}=block;
        n1=n1+1;
    end
end

% DIVIDING IMAGE INTO DOMAIN BLOCKS D
n2=1;
n3=1;
for k=1:step:256-D+1
    for h=1:step:256-D+1
            block2=I(k:k+D-1,h:h+D-1);
            Domain{n2}=imresize(block2, [R R]); % scaling of D blocks to size of R
        % Applying AFFINE Transformations on domain blocks D
        % 8 geometric affine transformations were performed on the domain blocks 
        % after scaling, so the number of domain blocks is 8*D_blocks.
        for pp=1:length(T)
            Newdomain{n3}=imwarp(Domain{n2},aff{pp});
            af_corr{n3}=aff{pp}; % saving affine transf.
            n3=n3+1;
        end
            n2=n2+1;
    end
end

% CALCULATION OF BRIGHTNESS AND CONTRAST COEFFICIENTS (oi, si)
% After applying affine transformations and dividing blocks into classes
% the coefficients of contrast si and brightness oi were calculated. 
% Thus, for each block R and the set of domain blocks after affine
% transformations, a pair of coefficients was determined.

% CODING
tic;
nummse=1;
n1=0;
for a=1:length(Range)
    for b=1:length(Newdomain)
        Ra=double(Range{a});
        Db=double(Newdomain{b});
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
        Dom{b}=uint8(Db.*s+mat);
        % calculation of MSE for blocks R and transformed blocks D, where 
        % the parameters of contrast and brightness are taken into account
        mse=immse(Range{a},Dom{b}); 
        err(b)=mse; % saving MSE
        coeff{b}=[s o]; % saving coefficients [si oi] into an array
        nummse=nummse+1; % number of block comparisons
    end
    % finding the minimum MSE, and the corresponding block D based on the
    % index In, the corresponding transformed block D, which is the most 
    % similar to the block R, was found.
    [M,In]=min(err); 
    aff_recons{a}=af_corr{In}; % saving corresponding affine transf.
    coordinates(a)=ceil(In/8); % position of D block before affine transf.
    s_o{a}=coeff{In}; % saving optimal blocks si & oi
end
% 
Coding_time=toc
Comparisons=nummse-1

% DEKODIRANJE
% Based on the transformed blocks D, affine transformations were found with
% the coefficients si & oi, which give the smallest MSE for corresponding 
% blocks R and D. Saved optimal coeff. and affine transformations were applied
% to the initial image to form the reconstructed image.

fname=uigetfile('*.gif');% opens ui for initial image
P=imread(fname); % reads chosen image
ReconsI=P; % initial image


for it=1:15 
    num2=1;
    for kk=1:step:256-D+1
        for hh=1:step:256-D+1
                blockk=ReconsI(kk:kk+D-1,hh:hh+D-1);
                Domain_dekod{num2}=imresize(blockk, [R R]);
                num2=num2+1;
        end
    end
    for f=1:length(Range)
    % applying affine transformation
    ReconsIm=imwarp(Domain_dekod{coordinates(f)},aff_recons{f});
    mat1=ones(size(ReconsIm)).*s_o{f}(2); 
    % optimization of contrast and brightness based on si and oi
    ReI{f}=imresize(uint8(double(ReconsIm).*s_o{f}(1)+mat1), [R R]);
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
    ReconsI=cell2mat(oldimg); % image in the next iteration
    if (it==1 || it==2 || it==3 || it==5)
        figure, imshow(ReconsI);
    end
end



% DISPLAY OF RECONSTRUCTED IMAGE
figure, imshow(ReconsI); title('Rekonstruisana slika');drawnow;
[peaksnr,~] = psnr(ReconsI,I);
PSNR=peaksnr