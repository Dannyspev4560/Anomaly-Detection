function C = TORP(datapath,ending,NOutliers,boolTreshold)    
    full_path = [datapath ending];
    imagefiles = dir(full_path);      
    numOfframes = length(imagefiles);
    graycells = cell(1,numOfframes-1) ;
    rgbcells = cell(1,numOfframes-1) ;
    
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
    
    %D
    M=double(InGray);% both inliers and outliers
    

    [U,S,Vt]=svd(M,'econ');
    V=Vt';
    Vnormal=norm((exp(1)*V),'fro');
    miu=Vnormal*(1/sqrt(0.8));
    ro=1/(128*miu^2*287);

    T=round(log((10*287*norm(M,2)) /(norm((eye(768)-U*U')*M,"fro"))));
    ro2=1/(128*mean(mean(M))^2*287);
    r=90;
    C=[];
    C(:,1)=0;
    L=[];
    
    for t=1:T
        [U,S,V]=svd(M-C(:,t));
        L=U*S*V';
        R=((eye(768)-U*U')*M);
        %Sinv=inv(S(1:287,:));
        E=pinv(S)*U'*M;
        col_norm=norm(M(:,t),2);
        %normE=norm(E())
    end

end