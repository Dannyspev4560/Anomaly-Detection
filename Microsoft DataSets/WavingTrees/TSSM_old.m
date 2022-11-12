imagefiles = dir('*.bmp');      
%nfiles = length(imagefiles);    % Number of files found
nfiles = 287;%inliers
graycells = cell(1,nfiles-1) ;
rgbcells = cell(1,nfiles-1) ;

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   %resize image by 1/8
   currentimage = imresize(currentimage,0.2);
   [rows,cols,dim] = size(currentimage);
   bwimage = reshape(rgb2gray(currentimage),[rows*cols,1]);
   graycells {ii} = bwimage;
   images{ii} = currentimage;
end
InGray = cell2mat(graycells );%B&W representation for 1 image in each row
No = 17;
Qi = 4;
Y=InGray;% both inliers and outliers

%[rows,cols]=size(Yi);
SSM1 = cell(1,cols);
SSM1_tild =  cell(1,cols);
delta = cell(1,cols);
a=1e-7;
[U,S,V] = svd(im2double(Y),'econ');
s=diag(S);
%stem(s);
Y = rescale(Y);
Pt = (Y'*Y);
Pt = rescale(Pt);
diag_load = a*trace(Y*Y');%should do rescale here?????????????????????????
diag_load_mat = diag_load*eye(size(Pt));
[rows,cols] = size(Pt);
PY_tild = Y*inv(Pt + diag_load_mat)*Y';
%PY_tild2 = 

[none,iters]=size(Y);
for k=1:iters % k iterator like i
    Yk = double(Y(:,k));
    SSM1{k} = norm(PY_tild*(Yk/norm(Yk)))^2;
    
end

SSM_MAT = cell2mat(SSM1);
SSM_MAT = rescale(SSM_MAT);
SSM_MAT = 1-SSM_MAT;%check why it flips
ssm_copy = SSM_MAT;
[rows,cols] = size(Y);
Y_tild = [];
for q=1:Qi
    [m,i]= max(ssm_copy);
    yq = Y(:,q);
    Y_tild = [yq Y_tild];
    ssm_copy(i) = 0;
end
Pt_tild = rescale(Y_tild'*Y_tild);
diag_load = a*trace(Y_tild*Y_tild');%should do rescale here?????????????????????????
diag_load_mat = diag_load*eye(size(Pt_tild));
PY_tild2 = Y_tild*inv(Pt_tild + diag_load_mat)*Y_tild';


for k=1:iters % k iterator like i
    Yk = double(Y(:,k));
    SSM1_tild{k} = norm(PY_tild2*(Yk/norm(Yk)))^2;
    
end
SSM_tild_MAT = cell2mat(SSM1_tild);
SSM_tild_MAT = rescale(SSM_MAT);
SSM_tild_MAT = 1 - SSM_tild_MAT;%check why it flips- same here
%stem(SSM_tild_MAT);
sorted_SSM_tild = sort(SSM_tild_MAT,'descend');
SSM_tilda_sum = sum(SSM_tild_MAT);
N = nfiles;
SSM_tilda_sum_fifth = sum(sorted_SSM_tild(1:uint8(N/5)));
quin_SSM = SSM_tilda_sum_fifth/SSM_tilda_sum;  %fix !!!!! too low if i take the lowest in the highest 20%
N = nfiles;
%std_div_SSM=std(SSM_tild_MAT);
%std_div_SSM = (1/sqrt(N))*sqrt((SSM_tilda_sum-(1/N)*(SSM_tilda_sum))^2);
Th = quin_SSM - std(SSM_tild_MAT);


outliers_indexes = [];
for k=1:N
    if SSM_tild_MAT(k)> Th
        outliers_indexes = [outliers_indexes k]
    end
end

disp(outliers_indexes);


