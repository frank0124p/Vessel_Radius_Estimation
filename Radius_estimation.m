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
temp(:,4)=1:1:7311;
temp(:,5)=pred;
for i=1:m
  if isnan(pred(i))
    temp(i,5)=0;
  end
end
[k,~]=find(temp(:,5));
[n,~]=size(k);
center_pts=zeros(n,5);
for i=1:n
   center_pts(i,1:3)=temp(k(i),1:3);
   center_pts(i,4)=temp(k(i),4); %original position
   center_pts(i,5)=temp(k(i),5); %connect position
end

%Get the connect node
connect_pair=cell(n,1);
for i=1:n
    temp_key=center_pts(i,4);%search the connect element
    [ind,~]=find(center_pts(:,5)==temp_key);
    if size(ind,1)>=2 | isempty(ind)==1
        if size(ind,1)>=2
                    connect_pair(i,1)={center_pts(ind,4)};

        else
                    temp_key=center_pts(i,5);
                    [ind,~]=find(center_pts(:,4)==temp_key);
                    if size(ind,1)>=2 | isempty(ind)==1
                        count=count+1;
                    else
                     connect_pair(i,1)={center_pts(ind,5)};
                    end

            
        end
    else
        connect_pair(i,1)={center_pts(ind,4)};
    end

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



%--------------------------------The smooth term reconstruct ------------------------
 load the_radius_matrix.mat














% Radius Estimation matrix 

% the_radius_matrix=zeros(n,1);

% %-------------------Radius Estimation  with two descriptor----------------------
% tic
% disp('Start Radius Estimation...');
% 
% 
% load Vein.mat %load the liver position first slice and last slice
% 
% Ori=smooth_Ori;
% Ori=Ori(:,:,firstslice:lastslice);
% 
% % %Use for the big branch Radius Estimation 
% % Ori=smooth_result;
% 
% 
% %-------------Using the Precontrast data -------------------------------
% % Pre_Ori=Pre_Ori(:,:,firstslice:lastslice);
% % Ori=Pre_Ori;
% %-------------------------------------------------------------------------------
% 
% for i=1:n
%  
% 
% 
%         %Inner circle Calculation 
%          Radius=8;
% 
%          [slice, ~,loc_x,loc_y,loc_z] =  extractSlice(Ori,center_pts(i,2),center_pts(i,1),center_pts(i,3),normal_vector(i,2),normal_vector(i,1),normal_vector(i,3),Radius);
%        %Processing the slice matrix without nan
%        if ~isempty(isnan(slice))
%         [row, col] = find(isnan(slice));
%         slice(row,col)=0;
%        end
%        %Calculate the Radius Flux  
%           scaleStep =0.5;
%           r = 1:scaleStep:5;
%           sigma = 1;
%           response = fastflux2(slice, r, sigma);
% 
%           MM = zeros(size(r));
%          
%         for j = 1:length(r)
%               MM(j) = max(max(response(:,:,j)));
%         end
% 
%         if max(MM)==0
%         estimate_radius=1;
%         else
%             if r(MM == max(MM))>5
%                         estimate_radius=5;
% 
%             else
%                         estimate_radius=r(MM == max(MM));
%             end
%         end
% 
%        
% %         [slice_plane, ~,outX,outY,outZ] =  extractSlice(Ori,center_pts(i,2),center_pts(i,1),center_pts(i,3),normal_vector(i,2),normal_vector(i,1),normal_vector(i,3),estimate_radius);
%         
%         
% 
%         disp(i);
% 
% 
%         %Add the smooth  term to determine the radius 
%         
%         the_radius_matrix(i)=estimate_radius;
% 
%         
% %         disp(['The label is ... ' num2str(estimate_radius)] );
% %     
% %         [k,g]=size(outX);
% %         outX=reshape(outX,k*g,1);outY=reshape(outY,k*g,1);outZ=reshape(outZ,k*g,1);
% %         outX=round(outX);    outY=round(outY);    outZ=round(outZ);
% % 
% % 
% %          [size_outx,~]=size(outX);
% % 
% %             for l=1:size_outx
% %              result(round(outX(l)),round(outY(l)),round(outZ(l)))=1;
% %             end
% %           
% 
%       
%      
% end     
% 
% 
% %--------------------------------------------------------------------------------------------------------------------
% 
% 
% 
%   
% 
% disp('End Radius Estimation...');
% 
% toc 


%Write the STL file
figure,imshow(max(result,[],3));

%Original Skeleton Data
fv = isosurface(skel,0.1);  
stlwrite('Ori_skel.stl', fv);


fv = isosurface(result,0.8);  
% fv=smoothpatch(fv,1,5);

stlwrite('Radius_Estimation.stl', fv);



allconnectnode = smooth3(result,'gaussian',[3 3 3]);
figure,imshow(max(allconnectnode,[],3));

fv = isosurface(allconnectnode,0.3);  

stlwrite('Radius_Estimation_gaussian.stl', fv);
