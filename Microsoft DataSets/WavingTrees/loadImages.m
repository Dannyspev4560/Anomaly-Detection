imagefiles = dir('*.bmp');      
%nfiles = length(imagefiles);    % Number of files found
nfiles = 287;%inliers
graycells = cell(1,nfiles-1) ;
rgbcells = cell(1,nfiles-1) ;
%%
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   %resize image by 1/8
   currentimage = imresize(currentimage,0.2);
   [rows,cols,dim] = size(currentimage);
   bwimage = reshape(rgb2gray(currentimage),[rows*cols,1]);
   graycells {ii} = bwimage;
   %Rim = reshape(currentimage(:,:,1),[rows*cols,1]);
   %Gim = reshape(currentimage(:,:,2),[rows*cols,1]);
   %Bim = reshape(currentimage(:,:,3),[rows*cols,1]);
   %RGB = reshape([Rim,Gim,Bim],[rows*cols*3,1]);
   %rgbcells{ii} = RGB;
   images{ii} = currentimage;
end
%%
InGray = cell2mat(graycells );%B&W representation for 1 image in each row
%InRGB = cell2mat(rgbcells );%RGB representation for 1 image in each row
No = 17;
Qi = 4;
Y=InGray;% both inliers and outliers
Yi = Y(:,1:240);
%YY=double(InGray)/255
%[A,~,~]=svd(YY,'econ')
[rows,cols]=size(Yi);
SSM = cell(1,cols);
%how to calculate those values they are to big- can i divide it by 255
%what alfa stands for- learning rate??
delta = cell(1,cols);
a=1e-7;
[U,S,V] = svd(im2double(Y),'econ');
s=diag(S);
%stem(s);
Pt = (Yk'*Yk);
Y_no=rescale(Y);
Pt2 = (Y_no'*Y_no);
%Pt22 = double(Pt2)/255;
Pt21 = rescale(Pt2);
[rows,cols] = size(Pt);
PY_tild = Yk*inv(Pt + delta_k*eye(rows,cols))*Yk';
%PY_tild2 = 
%%
[none,iters]=size(Y);
for k=1:iters % k iterator like i
     Yk = double(Y(:,k));
%     %P = (Yk'*Yk);
%     Pt = (Yk'*Yk);%hermetian operator is inverse- according to matlab func ishermetian - assure with Amir!!!
%     [U,S,V] = svd(Pt);
%     %Ptinv = inv(Pt);%composed poorly- need other way to inverse- how to apply svd????
%     %Pyk = Yk.*Ptinv.*(Yk');
%     delta_k = a*trace(Pt);
%     [rows,cols] = size(Pt);
%     %Matrix is close to singular or badly scaled - issue!!! - very small
%     %numbers calculated
%     Pytilda = Yk*inv(Pt + delta_k*eye(rows,cols))*Yk';% to correct singularity by adding small diagonal constant
%     %||~PY-PYi||^2=.....
    
    SSM{k} = norm(PY_tild*(Yk/norm(Yk)))^2;
    
end
%%
SSM_MAT = cell2mat(SSM);