function [S,S_old] = segmentation_3Dprocess(img_bac_sharpened_label,img_bac_sharpened,S,Estimated_core_size,cell_size_lower_thresh,erode_switch)
% This file is non-executable 

% Step 4: Core extraction and 3D Data curing

% The weird structue is for parallel calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is the core of the whole program. Please contact Yanlab, Yale(jing.yan@yale.edu)
% or write to haoran lu (hl2396@cornell.edu). We welcome any suggestions! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input: 
% img_bac_sharpened_label, img_bac_sharpened and S (See common parameter)
% Estimated_core_size and cell_size_lower_thresh: parameter controlling the segmentation criteria

% Output: 
% "S". The elements in output S represents the location information of each 'Repaired cells'(See common parameter)
% The old S structure is also recorded for visualization and comparation

S_old=S;
newcelli=zeros(1,size(S,2));% Record newly generated cells. It's a non-negative integer 
num=size(S,2);
parfor i=1:num % Here it can be changed between parfor and for 
    % Sreplace: New coordinates' matrix that raplace the old ones in the same position in S matrix
    % Sadd: new coordinates' matrix that will be list as new elements in S
    
    % Core extraction segmentation algorithm
    [Sreplace,nofnewcell,Sadd] = segmentation_3Dprocess_core_extraction(img_bac_sharpened_label,img_bac_sharpened,S,i,Estimated_core_size,cell_size_lower_thresh,erode_switch);
    % Data curing by mapping the deleted voxels back to the nearest cells
    [Sreplace_repaired,Sadd_repaired]=segmentation_3Dprocess_mapping(S{i},Sreplace,Sadd,nofnewcell);
    Sreplacemiddle{i}=Sreplace_repaired; % "middle" is for further substituting the old S
    Sreplacemiddle_old{i}=Sreplace;
    if nofnewcell>0
        Saddmiddle{i}=Sadd_repaired;
        Saddmiddle_old{i}=Sadd;
    end
    newcelli(i)=nofnewcell;% Numbers of cells that need to be replaced

    % Speed display
    disp('current process of 3D processing')
    disp((i)/num)
end

supplimentS=size(S,2)+1;% Cells that need to be added into S
num=size(S,2);
for i=1:num
        S_old{i}=[];
        S_old{i}=Sreplacemiddle_old{i};% replece those bacteria already in the structure of S
    if newcelli(i)>0
        for n=1:newcelli(i)
        S_old{supplimentS}=Saddmiddle_old{i}{n}(1:(Saddmiddle_old{i}{n}(1,1)+1),1:3);
        if newcelli(i)==1
            S_old{i}(1,2)=supplimentS;
        elseif newcelli(i)==2
            S_old{i}(1,3)=supplimentS;
        end
        supplimentS=supplimentS+1;
        end
    end
end
% Add newly seperated cells into S
supplimentS=size(S,2)+1;
for i=1:num
        S{i}=[];
        S{i}=Sreplacemiddle{i};% Replece those bacteria already in the structure of S
    if newcelli(i)>0 % Number i cell in the old S
        for n=1:newcelli(i)
        S{supplimentS}=Saddmiddle{i}{n}(1:(Saddmiddle{i}{n}(1,1)+1),1:3);
        if newcelli(i)==1
            S{i}(1,2)=supplimentS;
        elseif newcelli(i)==2
            S{i}(1,3)=supplimentS;
        end
        supplimentS=supplimentS+1;
        end
    end
end



