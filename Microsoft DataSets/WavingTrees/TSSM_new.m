imagefiles = dir('*.bmp');      
nfiles = 287;
graycells = cell(1,nfiles-1) ;
rgbcells = cell(1,nfiles-1) ;

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   %resize image by 1/8
   currentimage = imresize(currentimage,0.3);
   [rows,cols,dim] = size(currentimage);
   bwimage = reshape(rgb2gray(currentimage),[rows*cols,1]);
   graycells {ii} = bwimage;
   images{ii} = currentimage;
end
InGray = cell2mat(graycells );%B&W representation for 1 image in each row
No = 17;
Qi = 4;%need to calc it automtically
Y=InGray;% both inliers and outliers

SSM1 = cell(1,cols);
SSM1_tild =  cell(1,cols);
delta = cell(1,cols);
a=1e-7;
[U,S,V] = svd(im2double(Y),'econ');

s=diag(S);
Y = double(Y)/255;

%M = mean(Y,2);%%addidtion by Amir
%Y=Y-M;

Pt = (Y'*Y);
%Pt = rescale(Pt);
diag_load = a*trace(Y*Y');%should do rescale here?????????????????????????
[rows,cols] = size(Pt);
Py_tilda=Y*inv(Y'*Y +diag_load*eye(size(Pt)))*Y';
[none,iters]=size(Y);
for k=1:iters % k iterator like i
    Yk = double(Y(:,k));
    SSM1{k} = norm(Py_tilda*(Yk/norm(Yk)))^2;
    
end

SSM_MAT = cell2mat(SSM1);
%SSM_MAT = rescale(SSM_MAT);

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

title("Inliers with high values");
outliers_indexes = [];

sorted_SSM_tild = sort(SSM_tild_MAT,'descend');

N = nfiles;
std_div_SSM=std(SSM_tild_MAT);

Th = sorted_SSM_tild(uint8(N/5)) - std_div_SSM;

outliers_indexes = [];
for k=1:iters % k iterator like i
   if SSM_tild_MAT(k)< Th% if lower than Threshold it means outlier
       outliers_indexes = [outliers_indexes k] ;
   end
    
end

disp(outliers_indexes);


