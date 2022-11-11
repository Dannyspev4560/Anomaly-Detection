function [outliers_indexes] = TSSM_SSM(datapath,format,NOutliers,flag)
% To run the algorithm from the command window :
% TSSM_SSM('C:\....\','*.tif',1,1)
% important to add "\" at the end of the directory
% works with standard RGB formatts, fr .tif required excplicit adaption


%% Input:
% datapath_format = "C:\......\*.bmp ro *.jpg and so on
%datapath - path to images
%format - format of images- default is RGB , can process gray images too 
%NOutliers - for SSM -number of outliers, in TSSM can pass defeulat number
%flag - 1 ==> TSSM, other ==>SSM
%% Output:
%outliers_indexes - after SSM metric identified as outliers

%% Load Images local storage ( resolution reduction and conversion to gray included)
tic
full_path = [datapath format];
imagefiles = dir(full_path);      
numOfframes = length(imagefiles);
graycells = cell(1,numOfframes-1) ;

for ii=1:numOfframes
   currentfilename = imagefiles(ii).name;
   currentimage = imread([datapath currentfilename]);
   %resize image by 1/8
   currentimage = imresize(currentimage,0.2);%high res : 0.4
   [rows,cols,dim] = size(currentimage);
   if (dim > 1)
     currentimage =  rgb2gray(currentimage);
   end
   bwimage = reshape((currentimage),[rows*cols,1]);
   graycells {ii} = bwimage;
   images{ii} = currentimage;
end
InGray = cell2mat(graycells );%B&W representation for 1 image in each row
if NOutliers ~= 0
    No = NOutliers;
end

Y=InGray;% both inliers and outliers
Qi = calcDataRank(Y);%need to calc it automtically

SSM1 = cell(1,cols);
SSM1_tild = cell(1,cols);


%% first SSM matric calculation
a=1e-7;
Y = double(Y)/255;
Pt = (Y'*Y);
diag_load = a*trace(Y*Y');
[rows,cols] = size(Pt);
Py_tilda=Y*inv(Y'*Y +diag_load*eye(size(Pt)))*Y';
[~,iters]=size(Y);
for k=1:iters 
    Yk = double(Y(:,k));
    SSM1{k} = norm(Py_tilda*(Yk/norm(Yk)))^2;
end

SSM_MAT = cell2mat(SSM1);
ssm_copy = sort(SSM_MAT,'descend') ;
Y_tilda = ssm_copy(1:Qi);

%% second SSM matric calucation
Y_tild = [];
for q=1:Qi
    [~,i]= max(ssm_copy);
    yq = Y(:,i);
    Y_tild = [yq Y_tild];
    ssm_copy(i) = 0;
end
PY_tild2 = Y_tild*inv(Y_tild'*Y_tild + diag_load*eye(size(Y_tild'*Y_tild)))*Y_tild';

for k=1:iters 
    Yk = double(Y(:,k));
    SSM1_tild{k} = norm(PY_tild2*(Yk/norm(Yk)))^2;
    
end

%% print SSM metric for each images
SSM_tild_MAT = cell2mat(SSM1_tild);
stem(SSM_tild_MAT);
ylabel("SSM value");
xlabel("# of frame");
xlim([0 numOfframes]);
title("SSM value per each frame");

%yline(0.875,'-','Threshold');
outliers_indexes = [];

%% TSSM threshold calculation
if flag == 1
    sorted_SSM_tild = sort(SSM_tild_MAT,'descend');

    N = numOfframes;
    std_div_SSM=std(SSM_tild_MAT);

    Th = sorted_SSM_tild(uint8(N/5)) - std_div_SSM;

    outliers_indexes = [];
    for k=1:iters % k iterator like i
       if SSM_tild_MAT(k)< Th% if lower than Threshold it means outlier
           outliers_indexes = [outliers_indexes k] ;
       end 
    end
else 
    for k=1:No
        [~,i]= min(SSM_tild_MAT);%extracting the index of min value
        outliers_indexes = [outliers_indexes i];
        SSM_tild_MAT(i) = 1;
    end
end

disp(outliers_indexes);
disp("Threshold:  "+Th);
toc
end

%% Q rank calcuation 
function rankQ = calcDataRank(data_MAT)
    pctg=0.85;
    Y = double(data_MAT)/255;
    M = mean(Y,2);
    Y = Y-M;
    [U,S,V] = svd(im2double(Y),'econ');
    s = diag(S);
    cumsumS = cumsum(s);
    stem(cumsumS);
    xlabel("# of frames");
    ylabel("cumulative sum");
    title("Cumulative Sum of diag(S) as function of frame number");
    cumsumTH = pctg*cumsumS(length(cumsumS));%play here
    rankQ = length(find(cumsumS<cumsumTH));%low rank dim
    disp("total Q: "+rankQ)
  
    %% outliers and inliers rank calcuation for research (commented)
%      %-Q rank outliers
%     Yo=Y(:,1:48);
%     Mo = mean(Yo,2);%%addidtion by Amir
%     Yo = Yo-Mo;
%     [Uo,So,Vo] = svd(im2double(Yo),'econ');
%     so = diag(So);
%     cumsumSo = cumsum(so);
%     stem(cumsumSo);
%     cumsumTHo = pctg*cumsumSo(length(cumsumSo));%play here
%     rankQo = length(find(cumsumSo<cumsumTHo));%low rank dim
%     disp("Q outliers: "+rankQo)
% 
%     %Q rank inliers
%     %Yi=Y(:,[1:64 140:200]) ;
%     Yi=Y(:,49:200) ;
%     Mi = mean(Yi,2);%%addidtion by Amir
%     Yi = Yi-Mi;
%     [Ui,Si,Vi] = svd(im2double(Yi),'econ');
%     si = diag(Si);
%     cumsumSi = cumsum(si);
%     stem(cumsumSi);
%     cumsumTHi = pctg*cumsumSi(length(cumsumSi));%play here
%     rankQi = length(find(cumsumSi<cumsumTHi));%low rank dim
%     disp("Q inliers: "+rankQi)
    
end