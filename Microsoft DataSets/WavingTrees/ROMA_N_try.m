function [outliers] = ROMA_N(datapath,ending,NOutliers,boolTreshold)
% datapath_format = "C:\......\*.bmp ro *.jpg and so on
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

M=InGray;% both inliers and outliers
X=[];
%Normalzing matrix
[frameLen,~]=size(M(:,1));%n = 768
for i=1:numOfframes
    Xi = double(M(:,i));
    X(:,i) = Xi/(norm(Xi,2));
end

theta=[];
phi=[];
qi=[];
for i=1:numOfframes
    for j=1:numOfframes
        xi=X(:,i);
        xj=X(:,j);
        %CosTheta = max(min(dot(xi,xj)/(norm(xi)*norm(xj)),1),-1);
        %phi(i) = real(acosd(CosTheta));
        if i==j
            phi(i,j)=0;
            continue;
        else
            theta(i,j)= acos(xi' * xj);
            phi(i,j)=  acos(abs(xi' * xj));%
        end

        if length(qi)<i
            qi(i)=phi(i,j);
        else
            qi(i)=min(qi(i),phi(i,j));
        end
    end
end

N = numOfframes;
n = frameLen;
CN = 0;%(1/(1-1/(N^2*(N-1))))*  (1-1/(N^2*(N-1)))  ;%???

threshold = pi/2 -(CN/sqrt(n-2));

ang = 0:0.01:pi/2;
pow=n-2;
h = (2/sqrt(pi))*(gamma(n/2))*(1/(gamma((n-1)/2)))*(sin(ang).^pow);

%then: if (qi > threshold) --> outlier and vice versa

end