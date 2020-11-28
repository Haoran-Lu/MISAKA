function[Sreplace,nofnewcell,Sadd] = segmentation_3Dprocess_core_extraction(img_bac_sharpened_label,img_bac_sharpened,S,i,Estimated_core_size,cell_size_lower_thresh,erode_switch)
% This program is applied within the 3D process, which accelerate the whole process by making it runs parelled

% Inputs:
% img_bac_sharpened_label,img_bac_sharpened,S:Cluster information
% Estimated_core_size and cell_size_lower_thresh: parameter controlling the
% segmentation criteria

% Outputs:
% For each input cluster, it will output its corresponding cores.
% These cores are smaller than Estimated_core_size.

% Meaning of the parameters can be found in the parameter book
Smodifycontrol=1;% Judging whether the core is the first cell generated
nofnewcell=-1;% Whether this cluster can be split into more

% Variables for output
Sreplace=[];
Sadd=[];


len=S{i}(1,1)+1;
% Extract the small local labeled matrix
onecell=zeros(size(img_bac_sharpened_label));
onecell(img_bac_sharpened_label==i)=1; % pick the spaces occupied by the single cell.


%SS is the small coordinate matrix
SS=zeros(len-1,3);
for j=2:len
    SS(j-1,1:3)=S{i}(j,1:3);
end
xmin=min(SS(:,1))-1;
xmax=max(SS(:,1))+1;
ymin=min(SS(:,2))-1;
ymax=max(SS(:,2))+1;
zmin=min(SS(:,3))-1;
zmax=max(SS(:,3))+1;

% Extracting a small region to avoid operating on the whole image space
onecellnew=onecell(xmin:xmax,ymin:ymax,zmin:zmax);
onecell_repair=img_bac_sharpened(xmin:xmax,ymin:ymax,zmin:zmax).*onecellnew;


% Start modification loop
if size(SS,1)>Estimated_core_size
    % Each loop will increase the threshold slightly
    gap=(max(max(max(onecell_repair)))-min(min(min(onecell_repair))))/100;
    minimumlocal=min(min(min(onecell_repair)));
    for repeat=1:100
        thresholdlocal=gap*repeat+minimumlocal;
        
        % Increasing local intensity threshold 
        onecell_repair(onecell_repair(:,:,:)<thresholdlocal)=0;
        onecell_repair_binaryfilter=onecell_repair;
        onecell_repair_binaryfilter(onecell_repair_binaryfilter(:,:,:)>0)=1;
        
        if erode_switch==1
            if mod(repeat,2)==0 
                for isearch_layer=1:size(onecell_repair_binaryfilter,3)
                   [Coor_list_raw,~]=three_D_erosion(onecell_repair_binaryfilter(:,:,isearch_layer));
                end
                
                onecell_repair=onecell_repair_binaryfilter.*onecell_repair;
            else
                XX=bwconncomp(onecell_repair_binaryfilter(:,:,:),6);% for speed up
                Coor_list_raw=XX.PixelIdxList;%coordinates
            end
        else
            XX=bwconncomp(onecell_repair_binaryfilter(:,:,:),6);% for speed up
            Coor_list_raw=XX.PixelIdxList;%coordinates
        end
        
      
        
        
        onecell_repair_binaryfilter=onecell_repair*0;
        onecell_repair_small=onecell_repair_binaryfilter;
        %onecell_repairbackup=onecell_repair;
        
        mark_correct=0;mark_big=0;mark_small=0;
        onecell_repair_remain=onecell_repair_small;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Function part 2: Multiple new possible cores are generated %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if size(Coor_list_raw,2)>1 
            for innerclusternumber=1:size(Coor_list_raw,2)
                x=zeros(1,size(Coor_list_raw{innerclusternumber},1));y=x;z=x;
                %extract the coordinates
                xlength=size(onecell_repair_binaryfilter,1);
                ylength=size(onecell_repair_binaryfilter,2);
                testlocal=0*onecell_repair;
                for ilabel=1:size(Coor_list_raw{innerclusternumber},1)
                    z(ilabel)=1+floor(Coor_list_raw{innerclusternumber}(ilabel,1)/(xlength*ylength+0.00001));
                    y(ilabel)=1+floor((Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength))/(xlength+0.00001));
                    x(ilabel)=Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength)-(y(ilabel)-1)*(xlength);
                    testlocal(x(ilabel),y(ilabel),z(ilabel))=1;
                end
                
                %[x,y,z] = ind2sub(size(testlocal),find(testlocal == innerclusternumber));
                if size(x,2)<=Estimated_core_size && size(x,2)>cell_size_lower_thresh
                    mark_correct=1;
                    %A: extract out the needed small matrix
                    %A=testlocal;%/innerclusternumber;
                    %in the inner loop. onecell_repair2 serve
                    %as a matrix for deleting the corrected ones
                    onecell_repair_binaryfilter=onecell_repair_binaryfilter+testlocal;
                    %filling the S matrix
                    for newcellnum=1
                        if Smodifycontrol==1
                            Sreplace(newcellnum,1:3)=[size(x,2), 0, 0];
                        else
                            Sadd{nofnewcell+1}(newcellnum,1)=size(x,2);
                            Sadd{nofnewcell+1}(newcellnum,2)=i;
                        end
                    end
                    for newcellnum=2:size(x,2)+1
                        if Smodifycontrol==1
                            Sreplace(newcellnum,1:3)=[x(newcellnum-1)+xmin-1 y(newcellnum-1)+ymin-1 z(newcellnum-1)+zmin-1];
                        else
                            Sadd{nofnewcell+1}(newcellnum,1:3)=[x(newcellnum-1)+xmin-1  y(newcellnum-1)+ymin-1 z(newcellnum-1)+zmin-1];
                            
                        end
                    end
                    
                    %confirm Sreplace always only have one
                    %elements
                    if Smodifycontrol==1
                        Smodifycontrol=0;
                    end
                    if Smodifycontrol==0
                        nofnewcell=nofnewcell+1;%mark that this amount start from -1
                    end
                elseif size(x,2)<cell_size_lower_thresh+1
                    mark_small=1;
                    
                    x=zeros(1,size(Coor_list_raw{innerclusternumber},1));y=x;z=x;
                    %extract the coordinates
                    xlength=size(onecell_repair_binaryfilter,1);
                    ylength=size(onecell_repair_binaryfilter,2);
                    testlocal=0*onecell_repair;
                    for ilabel=1:size(Coor_list_raw{innerclusternumber},1)
                        z(ilabel)=1+floor(Coor_list_raw{innerclusternumber}(ilabel,1)/(xlength*ylength+0.00001));
                        y(ilabel)=1+floor((Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength))/(xlength+0.00001));
                        x(ilabel)=Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength)-(y(ilabel)-1)*(xlength);
                        testlocal(x(ilabel),y(ilabel),z(ilabel))=1;
                    end
                    onecell_repair_small=onecell_repair_small+testlocal;
                elseif size(x,2)>Estimated_core_size
                    mark_big=1;
                    
                    x=zeros(1,size(Coor_list_raw{innerclusternumber},1));y=x;z=x;
                    %extract the coordinates
                    xlength=size(onecell_repair_binaryfilter,1);
                    ylength=size(onecell_repair_binaryfilter,2);
                    testlocal=0*onecell_repair;
                    for ilabel=1:size(Coor_list_raw{innerclusternumber},1)
                        z(ilabel)=1+floor(Coor_list_raw{innerclusternumber}(ilabel,1)/(xlength*ylength+0.00001));
                        y(ilabel)=1+floor((Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength))/(xlength+0.00001));
                        x(ilabel)=Coor_list_raw{innerclusternumber}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength)-(y(ilabel)-1)*(xlength);
                        testlocal(x(ilabel),y(ilabel),z(ilabel))=1;
                    end
                    onecell_repair_remain=onecell_repair_remain+testlocal;%/innerclusternumber+remainfilter;
                end
                
            end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Function part 3: One new possible core is generated %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif size(Coor_list_raw,2)==1
            if size(Coor_list_raw{1},1)<=Estimated_core_size && size(Coor_list_raw{1},1)>cell_size_lower_thresh
                mark_correct=1;
                testlocal=onecell_repair;
                testlocal(testlocal(:,:,:)>0)=1;
                onecell_repair_binaryfilter=testlocal;
                %extract the coordinates
                x=[];y=[];z=[];
                xlength=size(onecell_repair_binaryfilter,1);
                ylength=size(onecell_repair_binaryfilter,2);
                
                for ilabel=1:size(Coor_list_raw{1},1)
                    z(ilabel)=1+floor(Coor_list_raw{1}(ilabel,1)/(xlength*ylength+0.00001));
                    y(ilabel)=1+floor((Coor_list_raw{1}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength))/(xlength+0.00001));
                    x(ilabel)=Coor_list_raw{1}(ilabel,1)-(z(ilabel)-1)*(xlength*ylength)-(y(ilabel)-1)*(xlength);
                    testlocal(x(ilabel),y(ilabel),z(ilabel))=1;
                end
                
                for newcellnum=1
                    if Smodifycontrol==1
                        Sreplace(newcellnum,1:3)=[size(x,2), 0, 0];
                    else
                        Sadd{nofnewcell+1}(newcellnum,1)=size(x,2);
                        Sadd{nofnewcell+1}(newcellnum,2)=i;
                    end
                end
                for newcellnum=2:size(x,2)+1
                    if Smodifycontrol==1
                        Sreplace(newcellnum,1:3)=[x(newcellnum-1)+xmin-1 y(newcellnum-1)+ymin-1 z(newcellnum-1)+zmin-1];
                    else
                        Sadd{nofnewcell+1}(newcellnum,1:3)=[x(newcellnum-1)+xmin-1  y(newcellnum-1)+ymin-1 z(newcellnum-1)+zmin-1];
                    end
                end
                
                if Smodifycontrol==1
                    Smodifycontrol=0;
                end
                if Smodifycontrol==0
                    nofnewcell=nofnewcell+1;%mark that this amount start from -1
                end
            elseif size(Coor_list_raw{1},1)>Estimated_core_size
                mark_big=1;
                testlocal=onecell_repair;
                testlocal(testlocal(:,:,:)>0)=1;
                onecell_repair_remain=onecell_repair_remain+testlocal;
            else
                mark_small=1;
                testlocal=onecell_repair;
                testlocal(testlocal(:,:,:)>0)=1;
                onecell_repair_small=testlocal;
            end
        end
        
        % Using notice variable to tell the program what operation need to
        % be done
        if mark_correct==1 && mark_big==0
            break
        elseif mark_correct==1
            onecell_repair(onecell_repair_binaryfilter(:,:,:)>0)=0;
        end
        if mark_big==1
            onecell_repair=onecell_repair.*onecell_repair_remain;
        end
        if mark_small==1
            onecell_repair(onecell_repair_small(:,:,:)>0)=0;
        end
        %%
    end

else
    Sreplace=zeros(size(SS,1)+1,3);
    Sreplace(1,1)=size(SS,1);
    Sreplace(2:(Sreplace(1,1)+1),1:3)=SS;
    nofnewcell=0;
end
