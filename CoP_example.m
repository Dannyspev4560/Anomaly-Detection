function [outliers] = CoP_example(datapath,format,NOutliers)
%% Input:
% datapath_format = "C:\......\*.bmp ro *.jpg and so on
%datapath - path to images
%format - format of images- default is RGB , can process gray images too 
%NOutliers - for SSM -number of outliers, in TSSM can pass defeulat number
%% Output:
%outliers - after SSM metric identified as outliers
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
    
 
    D=InGray;% both inliers and outliers
    r = rank(double(D));
  
    n=numOfframes;
                    
    [Uh,p,th] = CoP(D , n, r) ;   
    
    outliers= find(p <= th);
    disp(outliers);
    toc
    
end

function [U,p,th] = CoP(D,n,r)
    n = fix(n) ;
    r = fix(r) ;
    [N1,~] = size(D) ; 
    T = repmat(sum(D.^2).^0.5 , N1 , 1) ;
    X = double(D)./T ; 
    G = X'*X ; 
    G = G - diag(diag(G)) ;
    p = sum(G.^2) ; 
    p = p/max(p) ; 
    figure ; stem(p); title('The elements of vector p') ; 
    grid on ;
    th=calcThreshold(p);
    [~,b] = sort(p , 'descend') ;
    Y = X(:, b(1:n)) ; 
    [s,~,~] = svd(Y , 'econ') ; 
    U = s(: , 1:r) ;  

end

function rankQ = calcDataRank(data_MAT)
    pctg=0.85;
    Y = double(data_MAT)/255;
    M = mean(Y,2);
    Y = Y-M;
    [~,S,~] = svd(im2double(Y),'econ');
    s = diag(S);
    cumsumS = cumsum(s);
    %stem(cumsumS);
    cumsumTH = pctg*cumsumS(length(cumsumS));%play here
    rankQ = length(find(cumsumS<cumsumTH));%low rank dim
    %disp("total Q: "+rankQ)
end

function th=calcThreshold(p)
    averg=mean(p);
    sortedVec=sort(p);
    th=averg;

end