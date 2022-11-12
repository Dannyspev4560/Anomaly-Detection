clear global;
close all;
imtool close all;

% setup paths
addpath('util');
addpath('Bullwinkle');
addpath('signature_detectors');

% load one of the campus images
load muufl_gulfport_campus_3;


%crop hsi.Data
[row,col,dim] = size(hsi.Data);
hcube = hypercube(hsi.Data,hsi.info.wavelength);
Row = 145:180;
Column = 150:180;
newhcube = cropData(hcube,':',':',1:dim);
newhcube = cropData(newhcube,Row,Column,':');
hsiMat=newhcube.DataCube;
[row,col,dim] = size(hsiMat);
coloredImg = colorize(newhcube,"Method","rgb","ContrastStretching",true);
imtool(coloredImg)
%create 4 spectral signatures
%arr=[i for each norm(pixel-spectral signature)
%minimum identified as outliers - mark those pixeles run ssm on them


%%
%------------------------------------------------mapping outliers------
%red spectral band
red_spectral_signature=[];
red_spectral_signature(:,1) = squeeze(newhcube.DataCube(20,17,:));
red_spectral_signature(:,2) = squeeze(newhcube.DataCube(20,18,:));
red_spectral_signature(:,3) = squeeze(newhcube.DataCube(20,16,:));
red_spectral_signature(:,4) = squeeze(newhcube.DataCube(21,17,:));
red_spectral_signature(:,5) = squeeze(newhcube.DataCube(22,17,:));

red_spectral_signature_mean=[];
for i=1:72
    red_spectral_signature_mean(i,1) = sum(red_spectral_signature(i,:))/5;
end

red_outliers_map=[];
for i=1:row
    for j=1:col
        pixel_vec = squeeze(newhcube.DataCube(i,j,:));
        red_outliers_map(i,j) = norm(red_spectral_signature_mean(:)-pixel_vec);
    end
end

%blue spectral band
blue_spectral_signature=[];
blue_spectral_signature(:,1) = squeeze(newhcube.DataCube(12,8,:));
blue_spectral_signature(:,2) = squeeze(newhcube.DataCube(13,8,:));
blue_spectral_signature(:,3) = squeeze(newhcube.DataCube(12,9,:));
blue_spectral_signature(:,4) = squeeze(newhcube.DataCube(11,8,:));
blue_spectral_signature(:,5) = squeeze(newhcube.DataCube(12,7,:));

blue_spectral_signature_mean=[];
for i=1:72
    blue_spectral_signature_mean(i,1) = sum(blue_spectral_signature(i,:))/5;
end

blue_outliers_map = [];
for i=1:row
    for j=1:col
        pixel_vec = squeeze(newhcube.DataCube(i,j,:));
        blue_outliers_map(i,j) = norm(blue_spectral_signature_mean(:)-pixel_vec);
    end
end

%green spectral band
green_spectral_signature=[];
green_spectral_signature(:,1) = squeeze(newhcube.DataCube(9,14,:));
green_spectral_signature(:,2) = squeeze(newhcube.DataCube(9,13,:));
green_spectral_signature(:,3) = squeeze(newhcube.DataCube(9,15,:));
green_spectral_signature(:,4) = squeeze(newhcube.DataCube(10,14,:));
green_spectral_signature(:,5) = squeeze(newhcube.DataCube(9,13,:));

green_spectral_signature_mean=[];
for i=1:72
    green_spectral_signature_mean(i,1) = sum(green_spectral_signature(i,:))/5;
end

green_outliers_map = [];
for i=1:row
    for j=1:col
        pixel_vec = squeeze(newhcube.DataCube(i,j,:));
        green_outliers_map(i,j) = norm(green_spectral_signature_mean(:)-pixel_vec);
    end
end

%black spectral band
black_spectral_signature=[];
black_spectral_signature(:,1) = squeeze(newhcube.DataCube(25,10,:));
black_spectral_signature(:,2) = squeeze(newhcube.DataCube(25,9,:));
black_spectral_signature(:,3) = squeeze(newhcube.DataCube(25,11,:));
black_spectral_signature(:,4) = squeeze(newhcube.DataCube(26,10,:));
black_spectral_signature(:,5) = squeeze(newhcube.DataCube(24,10,:));

green_spectral_signature_mean=[];
for i=1:72
    black_spectral_signature_mean(i,1) = sum(black_spectral_signature(i,:))/5;
end

black_outliers_map = [];
for i=1:row
    for j=1:col
        pixel_vec = squeeze(newhcube.DataCube(i,j,:));
        black_outliers_map(i,j) = norm(black_spectral_signature_mean(:)-pixel_vec);
    end
end




outliers_red_map = red_outliers_map < 0.25;
%imshow(outliers_red_map)

outliers_blue_map = blue_outliers_map < 0.25;
%imshow(outliers_blue_map)

outliers_green_map = green_outliers_map < 0.25;
%imshow(outliers_green_map)

outliers_black_map = black_outliers_map < 0.25;
%imshow(outliers_black_map)


subplot(2,2,1), imshow(outliers_blue_map),title('blue');
subplot(2,2,2), imshow(outliers_green_map),title('green');
subplot(2,2,3), imshow(outliers_black_map),title('black');
subplot(2,2,4), imshow(outliers_red_map),title('red');
figure;
outliers_map = outliers_red_map+outliers_blue_map+outliers_green_map+ outliers_black_map;
imtool(outliers_map), title('outliers');


%%
Y_copy = hsiMat;
Y_copy(:,:,dim+1)=outliers_map;


Y_copy = permute(Y_copy, [3 2 1]); %% change order for correct reshape
Y_copy=reshape(Y_copy,[dim+1,row*col]);% each col represents pixel+ map
Y = permute(hsiMat, [3 2 1]);
Y=reshape(Y,[dim,row*col]);%to rows:(36*31)*72==>1116(cols)*72(rows)
%%return back: 
   %Y_copy = reshape(Y_copy,[dim+1,col,row]);
   %Y_copy = permute(Y_copy, [3 2 1]);
%problem- values ??? needed to be in range of or 0-1 ???? negative values
%Y=rescale(Y)
[rows,cols] = size(Y);
SSM1 = cell(1,cols);
SSM1_tild = cell(1,cols);
delta = cell(1,cols);
Pt = (Y'*Y);
%Pt = rescale(Pt);
a=1e-5;%increased value due to singularity issue
diag_load = a*trace(Y*Y');
[rows,cols] = size(Pt);
Py_tilda=Y*inv(Y'*Y +diag_load*eye(size(Pt)))*Y';%values too big
[none,iters]=size(Y);
for k=1:iters % k iterator like i
    Yk = double(Y(:,k));
    SSM1{k} = norm(Py_tilda*(Yk/norm(Yk)))^2;
    
end

SSM_MAT = cell2mat(SSM1);
%SSM_MAT = rescale(SSM_MAT);
stem(SSM_MAT);

%---------Q calc--------------------
M = mean(Y,2);%%addidtion by Amir
Y = Y-M;
[U,S,V] = svd(im2double(Y),'econ');
s = diag(S);
cumsumS = cumsum(s);
stem(cumsumS);
cumsumTH = 0.7*cumsumS(length(cumsumS));%play here
Q = length(find(cumsumS<cumsumTH));%low rank dim
disp("Q: "+Q);




Yi = [];

 
Yo = [];

ii=1;
oo=1;

for i=1:rows
    if Y_copy(73,i)==1
        Yo(:,oo)=Y(:,i);
        oo=oo+1;
    else
        Yi(:,ii)=Y(:,i);
        ii=ii+1;
    end
end
%----Q inliers
Mi = mean(Yi,2);%%addidtion by Amir
Yi = Yi-Mi;
[Ui,Si,Vi] = svd(im2double(Yi),'econ');
si = diag(Si);
cumsumSi = cumsum(si);
stem(cumsumSi);
cumsumTHi = 0.7*cumsumSi(length(cumsumSi));%play here
Qi = length(find(cumsumSi<cumsumTHi));%low rank dim
disp("Qi: "+Qi);




%----Q outliers
Mo = mean(Yo,2);%%addidtion by Amir
Yo = Yo-Mo;
[Uo,So,Vo] = svd(im2double(Yo),'econ');
so = diag(So);
cumsumSo = cumsum(so);
stem(cumsumSo);
cumsumTHo = 0.7*cumsumSo(length(cumsumSo));%play here
Qo = length(find(cumsumSo<cumsumTHo));%low rank dim
disp("Qo: "+Qo);


%--------------

ssm_copy = sort(SSM_MAT,'descend') ;
Y_tilda = ssm_copy(1:Qi);

Y_tild = [];
for q=1:Qi
    [m,i]= max(ssm_copy);
    yq = Y(:,q);
    Y_tild = [yq Y_tild];
    ssm_copy(i) = 0;
end
PY_tild2 = Y_tild*inv(Y_tild'*Y_tild + diag_load*eye(size(Y_tild'*Y_tild)))*Y_tild';

for k=1:iters % k iterator like i
    Yk = double(Y(:,k));
    SSM1_tild{k} = norm(PY_tild2*(Yk/norm(Yk)))^2;
    
end
SSM_tild_MAT = cell2mat(SSM1_tild);
%SSM_tild_MAT = rescale(SSM_tild_MAT);
stem(SSM_tild_MAT);
ylabel("SSM value");
xlabel("# of pixel");
%yline(0.4664,'-','Threshold')

title("SSM evaluation ");

sorted_SSM_tild = sort(SSM_tild_MAT,'descend');

N = length(Y);
std_div_SSM=std(SSM_tild_MAT);

Th = sorted_SSM_tild(uint8(N/9)) - std_div_SSM;

outliers_indexes = [];
for k=1:iters % k iterator like i
   if SSM_tild_MAT(k)< Th% if lower than Threshold it means outlier
       outliers_indexes = [outliers_indexes k] ;
   end 
end

%disp(outliers_indexes);
disp("threshold")
disp(Th);


cubeSSM = reshape(SSM_tild_MAT,[1,col,row]);
cubeSSM = permute(cubeSSM, [3 2 1]);%for visualization-removing the 3rd(72) dim

for i=1:row
    for j=1:col
        if cubeSSM(i,j)>= Th
            cubeSSM(i,j) = 0;
        else
            cubeSSM(i,j) = 1;%outlier painted in white
        end
    end
end

imtool(cubeSSM)

Tn=0;
Tp=0;
Fn=0;
Fp=0;
outliers=0;
cubeSSM = double(cubeSSM);

for i=1:row
    for j=1:col
        if outliers_map(i,j)==1
            outliers= outliers+1;
        end
        if cubeSSM(i,j) == outliers_map(i,j) &&  outliers_map(i,j) == 1
            Tp=Tp+1;
        end
        if cubeSSM(i,j) == outliers_map(i,j) &&  outliers_map(i,j) == 0
            Fp=Fp+1;
        end
        if cubeSSM(i,j) ~= outliers_map(i,j) &&  outliers_map(i,j) == 1
            Fn=Fn+1;
        end
        if cubeSSM(i,j) ~= outliers_map(i,j) &&  outliers_map(i,j) == 0
            Tn=Tn+1;
        end
    end
end

disp((Tp+Tn)/(Tp+Tn+Fn+Fp));