function [inlierest_ns,inlierest,outlierest,Uhat] = ROMA_N(datapath,ending,flag)
%% Input: 
%M - data matrix of size nxN - n - dimension and N - number of vectors 
%flag - whether to adapt threshold or not (flag=1 - for adaptation )
%% Output:
%inlierest_ns - inlier index set estimate after removing only unstructured
%outliers
%inlierest - inier index set estimate
%outlierest - outlier index set estimate
%Uhat - estimated subsapce bases
%% Find all the angles - and the mean principle angle
    tic
    normc_fcn = @(M) sqrt(M.^2 ./ sum(M.^2));
    full_path = [datapath ending];
    imagefiles = dir(full_path);      
    numOfframes = length(imagefiles);
    graycells = cell(1,numOfframes-1) ;
    rgbcells = cell(1,numOfframes-1) ;
    
    for ii=1:numOfframes
       currentfilename = imagefiles(ii).name;
       currentimage = imread([datapath currentfilename]);
       %resize image by 1/8
       currentimage = imresize(currentimage,0.8);%high res : 0.4
       [rows,cols,dim] = size(currentimage);
       if (dim > 1)
         currentimage =  rgb2gray(currentimage);
       end
       bwimage = reshape((currentimage),[rows*cols,1]);
       graycells {ii} = bwimage;
       images{ii} = currentimage;
    end
    InGray = cell2mat(graycells );%B&W representation for 1 image in each row
    
    M=InGray;% both inliers and outliers
    [n,N] = size(InGray);
    alp=1/(N^2);
    %nm = normc(M);%normalize matrix foe angle computation
    X=[];
    [frameLen,~]=size(M(:,1));
    for i=1:numOfframes
        Xi = double(M(:,i));
        X(:,i) = Xi/(norm(Xi));
    end
    nm = X;
    c = acos(abs(nm'*nm));% angles - phi_ij 's
    c = abs(c);
    c1=c;
    for i=1:N
        c1(i,i) = 10;
    end
    % Calculating mean of all angles for adaptation
    if flag ==1
        cp = acos(nm'*nm); % angles - theta_ij 's 
        cp = abs(cp);
        c2 = zeros(N*(N-1)/2,1);
        k=1;
        for i=1:N-1
            for j=i+1:N
                c2(k) = cp(i,j);
                k=k+1;
            end
        end
        mean_angle = mean(c2);
    else
        mean_angle = pi/2;
    end
    sig_angle = 1/sqrt(n-2);
%% Find minimum angle formed by each point and also the mean of all angles
    cmin = min(c1);
    [~,i_pos]  = min(cmin);% Reference inlier point
    c_ipos = c(i_pos,:);
    [~,o_pos] = max(c_ipos);% Reference outlier point
    %% Threshold computation 
    thrnew  = mean_angle-(norminv(1-(alp/(2*(N-1))))*sig_angle);
    if (thrnew<0)
       thrnew  = pi/2-(norminv(1-(alp/(2*(N-1))))*sig_angle);
    end
    %% Find inlier indices
    inlierest_ns = find(cmin<=thrnew);% After removing unstructured outliers
    %% Find the score of each point - number of angles above threshold
    scr = zeros(1,N);
    for i=1:N
       scr(i) = length(find(c(i,inlierest_ns)>thrnew));
    end
    %% Recompute inliers and outliers using new metric  
    inlierest = [];
    outlierest = find(cmin>thrnew);
    % Find a reference structured outlier from the current inlier set
    while (isempty(find(inlierest_ns==o_pos, 1)))
       c_ipos(o_pos) = 0;
       [~,o_pos] = max(c_ipos);
    end
    %% Calculate inlier indices removing structured outliers
    for i=inlierest_ns
        if abs(scr(i)-scr(i_pos))<abs(scr(i)-scr(o_pos))
            inlierest = [inlierest,i];
        else
            outlierest = [outlierest,i];
        end
    end
    %% Estimate the subspace using SVD
    
%     Y = nm(:,inlierest);
%     r1 = min(size(Y));
%     [U,S1,~] = svd(Y);
%     Uhat=[];
%      for i=1:r1
%          if S1(i,i)>10^-7
%              Uhat = [Uhat,U(:,i)];
%          end
%      end
    toc
end