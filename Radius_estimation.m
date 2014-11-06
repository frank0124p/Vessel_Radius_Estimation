clear;close all;
path=pwd;
path=strcat(pwd,'\Data');
addpath(path);
% Initial  load centerline_points.mat 
load tubular_response_bigBranch.mat %parameter --smooth_result with big branch
load no_zero.mat
load centerline_gaussian_smooth.mat % parameter-result 
skel=allconnectnode;
result=allconnectnode;
load smooth_Ori.mat %The Dicom image
load connectPairTree.mat %pred for the minimal spanning tree structure


% load Precontrast.mat
% Pre_Ori=Ori;
%ignore the root node due to we can judge the connect node from pred


[m,~]=size(no_zero);
temp=zeros(m,5);
temp(:,1:4)=no_zero;
temp(:,5)=pred;
for i=1:m
  if isnan(pred(i))
    temp(i,5)=0;
  end
end
[k,~]=find(temp(:,5));
[n,~]=size(k);
center_pts=zeros(n,3);
for i=1:n
   center_pts(i,:)=temp(k(i),1:3);
end


%-----------------Calculate the normal_vector---------------------------
%use the connect node to difference the normal vector
temp=ones(m,4);
temp(1,4)=0;
for i=2:m
    if isnan(pred(i))
      temp(i,4)=0; 
    else
      temp(i,1:3)=no_zero(i,1:3)-no_zero(pred(i),1:3);
    end
end
[g,~]=find(temp(:,4));
[n,~]=size(g);

normal_vector=ones(n,3);
for i=1:n
    normal_vector(i,1:3)=temp(g(i),1:3);
end 

normal_vector=abs(normal_vector);


%-------------------Radius Estimation  with two descriptor----------------------
tic
disp('Start Radius Estimation...');


load Vein.mat %load the liver position first slice and last slice

Ori=smooth_Ori;
Ori=Ori(:,:,firstslice:lastslice);


%-------------Using the Precontrast data -------------------------------
% Pre_Ori=Pre_Ori(:,:,firstslice:lastslice);
% Ori=Pre_Ori;
%-------------------------------------------------------------------------------

for i=1:n
 


        %Inner circle Calculation 
          Radius=15;

         [slice, ~,loc_x,loc_y,loc_z] =  extractSlice(Ori,center_pts(i,2),center_pts(i,1),center_pts(i,3),normal_vector(i,2),normal_vector(i,1),normal_vector(i,3),Radius);
       %Calculate the Radius Flux  
          scaleStep =1;
          r = 1:scaleStep:10;
          sigma = 1;
          response = fastflux2(slice, r, sigma);

          MM = zeros(size(r));
         
        for j = 1:length(r)
              MM(j) = max(max(response(:,:,j)));
        end
        
        
        figure;
        plot(r, MM, 'linewidth', 2)
        hold on;
        plot(r(MM == max(MM)), max(MM), 'or', 'markersize', 8)
    
%     disp(['The label is ... ' num2str(label)] );
%     
%     [k,g]=size(outX);
%     outX=reshape(outX,k*g,1);outY=reshape(outY,k*g,1);outZ=reshape(outZ,k*g,1);
%     outX=round(outX);    outY=round(outY);    outZ=round(outZ);
% 
% 
%      [length,~]=size(outX);
% 
%       if label==1
%         for l=1:length
%          result(round(outX(l)),round(outY(l)),round(outZ(l)))=0.17;
%         end
%       end
%       
       disp(i);
      
     
end     


%--------------------------------------------------------------------------------------------------------------------
%Using the local hessian matrix to determine the radius 


% for i=1:m
%     
%      [inner_slice, ~,inner_x,inner_y,inner_z] =extractSlice(Ori,center_pts(i,2),center_pts(i,1),center_pts(i,3),normal_vector(i,2),normal_vector(i,1),normal_vector(i,3),3);
%      [a,b]=size(inner_slice);
%      slicedouble=zeros(a,b,3);
%      slicedouble(:,:,1)=inner_slice;
%      slicedouble(:,:,2)=inner_slice;
%      slicedouble(:,:,3)=inner_slice;
% 
% %    options.ItenValue = MaxValue(:,2);
%      options.FrangiScaleRange =[1 10];
%      options.BlackWhite = false;
%      options.FrangiScaleRatio = 0.5; 
% 
% 
%      [FilterOut,whatScale,Voutx,Vouty,Voutz] = VesselnessFilter(slicedouble,options);
% 
% 
% 
% %      IMout(IMout>0.9)=0;
% 
% 
% 
% 
% 
%      hsigma = 2.5;
%      hsize = round((6*hsigma - 1)/2); 
%      h = fspecial('gaussian', hsize, hsigma);
%      G_Im = imfilter(slicedouble,h,'replicate');
%      [FinalclusterIM,Cluster_Max_Info] = ResponsePostProcessing(FilterOut,Voutx,Vouty,Voutz,G_Im);
%         
%      
%      temp=FinalclusterIM(:,:,1);
%      ind=find(temp>0);
%      if ~isempty(ind)
%       [length,~]=size(ind);
%       [k,g]=size(inner_x);
%       outX=reshape(inner_x,k*g,1);outY=reshape(inner_y,k*g,1);outZ=reshape(inner_z,k*g,1);
%       outX=round(outX);    outY=round(outY);    outZ=round(outZ);
%        for j=1:length
%            result(outX(ind(j)),outY(ind(j)),outZ(ind(j)))=0.12;
%        end
%      end
%       disp(i);
%     
% end



  

disp('End Radius Estimation...');

toc 


%Write the STL file
figure,imshow(max(result,[],3));


fv = isosurface(skel,0.1);  
% fv=smoothpatch(fv,1,5);

stlwrite('Ori_skel.stl', fv);


fv = isosurface(result,0.1);  
fv=smoothpatch(fv,1,5);

stlwrite('Radius_Estimation.stl', fv);



allconnectnode = smooth3(result,'gaussian',[3 3 3]);
figure,imshow(max(allconnectnode,[],3));

fv = isosurface(allconnectnode,0.08);  
fv2=smoothpatch(fv,1,1);

stlwrite('Radius_Estimation_smooth.stl', fv2);
