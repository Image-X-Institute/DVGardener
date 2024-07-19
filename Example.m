% ------------------------------------------
% FILE   : Example.m
% AUTHOR : Owen Dillon, The University of Sydney
% DATE   : 2024-05-20  Created.
% ------------------------------------------
% PURPOSE
%   This script runs through an example of DVGardener.
%   This was tested on Matlab R2021b
%   Assumes you have a working installation of Elastix
%   https://github.com/SuperElastix/elastix/releases/tag/5.1.0
%   And have downloaded images from TCIA
%   https://www.cancerimagingarchive.net/collection/4d-lung/
%   As well as the matlab modules of this workshop added to your path.
%   We use elastix to register 4DCT images around the respiratory cycle
%   First as inhale to exhale, then phase 0 to 1 to 2 ... to 9 to 1.
%   We compare the inhale to exhale DVF to the 0 to 1 ... to 5 DVF (as in
%   thory these will be the same)
%   We show how to impose consistency across the cycle (0 to ... to 0 = 0)
%   and average consistency (0 to .. to 5 = to exhale AND ... to 0 = 0)
%   We generate perturbations of the inhale/exhale DVF, one from the
%   circular registration, one by gaussian sampling.
%
%
% ------------------------------------------
% Initialise
close all; clear all;
% toggles
do_dcmtomha = false;
do_figures = false;
do_elastix_inex = false;
do_elastix_circle = false;
do_consistent = false;
do_compare = true;
do_compress = false;
% parameters
circ_lam = 20000;             %regularisation parameter for circular consistency
mask_mid = [256,110,190];        %rectangular mask middle
mask_rad = [100,35,100];      %mask width
% Define locations
name_elastix = ['C:\elastix\elastix.exe'];
name_transformix = ['C:\elastix\transformix.exe'];
elastix_parfile = 'Owens_Elastix_BSpline_OpenCL';
elastix_parfile_full = ['C:\elastix\parfiles\Owens_Elastix_BSpline_OpenCL.txt'];
%where we keep the data from TCIA. Note this was developed for
% Subject ID 100_HM10395
% Study UID 1.3.6.1.4.1.14519.5.2.1.6834.5010.335014117706890137582032169351
% Series ID 1.3.6.1.4.1.14519.5.2.1.6834.5010.124525254231969293228570633659
% We renamed each dicom folder to g0, g1 for 0% gating, 10% gating etc. to
% avoid an error in windows with directories having spaces.
dir_rawdata = ['C:\DVGardener\Workshop\pat101\pat101\4D-Lung\101_HM10395\10-21-1997-NA-p4-86157'];
dir_data = ['C:\DVGardener\Workshop\input'];
dir_out = ['C:\DVGardener\Workshop\output'];
%% Convert .dcm data from TCIA to .mha files
if do_dcmtomha;
    mkdir(dir_data);
    dirs_raw = ls(dir_rawdata); dirs_raw = cellstr(dirs_raw);
    dirs_raw(1,:) = []; dirs_raw(1,:) = [];
    n_phases = length(dirs_raw);
    for jj = 1:n_phases;
        name_out = [dir_data,'\\CT_',num2str(jj-1),'.mha'];
        names_dicoms = ls([dir_rawdata,'\',dirs_raw{jj,1}]);
        names_dicoms = cellstr(names_dicoms);
        names_dicoms(1,:) = []; names_dicoms(1,:) = [];
        n_slices = length(names_dicoms);
        vol = [];
        for kk = 1:n_slices;
            name_slice = [dir_rawdata,'\',dirs_raw{jj,1},'\',...
                names_dicoms{kk,1}];
            info_raw{jj,kk} = dicominfo(name_slice);
            indx_slice(kk) = info_raw{jj,kk}.InstanceNumber;
            slice = dicomread(name_slice);
            slice = rot90(slice,2);
            vol = cat(3,vol,slice);
        end
        % Reorder (just in case)
        vol = vol(:,:,indx_slice);
        order = [2,3,1];
        vol = permute(vol,order);
        voxel_size = [info_raw{1,1}.PixelSpacing;info_raw{1,1}.SliceThickness]';
        voxel_size = voxel_size(order);
        % Build the .MHA header
        load MhaHeaderTemplate;
        info_mha = mhaHeaderTemplate;
        dims = size(vol);
        info_mha.PixelDimensions = voxel_size;
        info_mha.Dimensions = dims;
        info_mha.Offset = -0.5*(dims.*voxel_size);
        % Write
        MhaWrite(info_mha,vol,name_out);
    end
end
%% plotting some things.
% Note this makes use of imshow3Dfull by Maysam Shahedi
% https://au.mathworks.com/matlabcentral/fileexchange/47463-imshow3dfull
if do_figures
    [info_ct,vol_ct] = MhaRead([dir_data,'\CT_0.mha']);
    mask = logical(vol_ct); mask = false&mask;
    mask(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3)) = true;
    vol_ct(mask) = vol_ct(mask)+300;
    figure; imshow3Dfull(vol_ct);
end
%% Register inhale CT to exhale
if do_elastix_inex;
    %0 is inhale, 5 is exhale
    mkdir([dir_out,'\inex']);
    mkdir([dir_out,'\inex\elastix']);
    mkdir([dir_out,'\inex\dvf']);
    name_elastixpams = [dir_out,'\inex\elastix\CT_0_to_CT_5.txt'];
    name_elastixresult = [dir_out,'\inex\elastix\CT_0_to_CT_5.mha'];
    name_dvf = [dir_out,'\inex\dvf\dvf_0_to_5.mha'];
    % do the registration
    [fail, msg_elastix] = system([name_elastix,' ',...
        ' -f ',dir_data,'\\CT_5.mha',...
        ' -m ',dir_data,'\\CT_0.mha',...
        ' -out ',dir_out,'\inex\elastix ',...
        ' -p ',elastix_parfile_full,...
        ]);
    copyfile([dir_out,'\inex\elastix\result.0.mha'],name_elastixresult);
    copyfile([dir_out,'\inex\elastix\TransformParameters.0.txt'],name_elastixpams);
    delete([dir_out,'\inex\elastix\result.0.mha']);
    delete([dir_out,'\inex\elastix\TransformParameters.0.txt']);
    % convert elastix parameter file to a DVF
    [fail, msg_transformix] = system([name_transformix,' ',...
        ' -tp ',name_elastixpams,...
        ' -out ',dir_out,'\inex\dvf ',...
        ' -def all ',...
        ]);
    copyfile([dir_out,'\inex\dvf\deformationField.mha'],name_dvf);
    delete([dir_out,'\inex\dvf\deformationField.mha']);
end
%% register around the respiratory cycle
if do_elastix_circle
    %0 is inhale, 5 is exhale
    mkdir([dir_out,'\circ']);
    mkdir([dir_out,'\circ\elastix']);
    mkdir([dir_out,'\circ\dvf']);
    for jj = 1:9;
        names_elastixpams{1,jj} = [dir_out,'\circ\elastix\',...
            'CT_',num2str(jj-1),'_to_CT_',num2str(jj),'.txt'];
        names_elastixresults{1,jj} = [dir_out,'\circ\elastix\',...
            'CT_',num2str(jj-1),'_to_CT_',num2str(jj),'.mha'];
        names_dvfs{1,jj} = [dir_out,'\circ\dvf\',...
            'dvf_',num2str(jj-1),'_to_',num2str(jj),'.mha'];
    end
    names_elastixpams{1,10} = [dir_out,'\circ\elastix\',...
        'CT_9_to_CT_0.txt'];
    names_elastixresults{1,10} = [dir_out,'\circ\elastix\',...
        'CT_9_to_CT_0.mha'];
    names_dvfs{1,10} = [dir_out,'\circ\dvf\',...
        'dvf_9_to_0.mha'];
    % do the registrations
    for jj = 1:10
        if jj<10;
            [fail, msg_elastix] = system([name_elastix,' ',...
                ' -f ',dir_data,'\\CT_',num2str(jj),'.mha',...
                ' -m ',dir_data,'\\CT_',num2str(jj-1),'.mha',...
                ' -out ',dir_out,'\circ\elastix ',...
                ' -p ',elastix_parfile_full,...
                ]);
        else
            [fail, msg_elastix] = system([name_elastix,' ',...
                ' -f ',dir_data,'\\CT_0.mha',...
                ' -m ',dir_data,'\\CT_9.mha',...
                ' -out ',dir_out,'\circ\elastix ',...
                ' -p ',elastix_parfile_full,...
                ]);
        end
        copyfile([dir_out,'\circ\elastix\result.0.mha'],names_elastixresults{1,jj});
        copyfile([dir_out,'\circ\elastix\TransformParameters.0.txt'],names_elastixpams{1,jj});
        delete([dir_out,'\circ\elastix\result.0.mha']);
        delete([dir_out,'\circ\elastix\TransformParameters.0.txt']);
        % convert elastix parameter file to a DVF
        [fail, msg_transformix] = system([name_transformix,' ',...
            ' -tp ',names_elastixpams{1,jj},...
            ' -out ',dir_out,'\circ\dvf ',...
            ' -def all ',...
            ]);
        copyfile([dir_out,'\circ\dvf\deformationField.mha'],names_dvfs{1,jj});
        delete([dir_out,'\circ\dvf\deformationField.mha']);
    end
end
%% make the circular registration consistent i.e. should end near 0
%need to add up DVFs along the cycle - do this by applying DVF to itself
%cumulatively
% load dvf j, add to the sum, move sum by dvf j.
% we do regularised end consistency, and inhale/exhale consistency
if do_consistent
    mkdir([dir_out,'\circ\cons']);
    %compose the DVFs
    for jj = 1:9;
        names_elastixpams{1,jj} = [dir_out,'\circ\elastix\',...
            'CT_',num2str(jj-1),'_to_CT_',num2str(jj),'.txt'];
        names_dvfs{1,jj} = [dir_out,'\circ\dvf\',...
            'dvf_',num2str(jj-1),'_to_',num2str(jj),'.mha'];
    end
    names_elastixpams{1,10} = [dir_out,'\circ\elastix\',...
        'CT_9_to_CT_0.txt'];
    names_dvfs{1,10} = [dir_out,'\circ\dvf\',...
        'dvf_9_to_0.mha'];
    % have to break the DVF into 3 "images" before applying transform
    [info,dvf] = MhaRead(names_dvfs{1,5});
    dvf_sum = 0*dvf;      %sum of DVF elements so far
    infok = info;
    infok = rmfield(infok,'ElementNumberOfChannels');
    for jj = 1:10;
        [info,dvf] = MhaRead(names_dvfs{1,jj});
        dvf_sum = dvf_sum+dvf;
        for kk = 1:3;
            %write one component of DVF sum to an "image"
            MhaWrite(infok,dvf_sum(:,:,:,kk),[dir_out,'\circ\cons\dvf_',num2str(kk),'.mha']);
            %apply phase j DVF to sum
            [fail, msg_transformix] = system([name_transformix,' ',...
                ' -tp ',names_elastixpams{1,1},...
                ' -in ',dir_out,'\circ\cons\dvf_',num2str(kk),'.mha ',...
                ' -out ',dir_out,'\circ\cons ',...
                ]);
            % load sum so far (so it is now in the correct spot
            [infok,dvfk] = MhaRead([dir_out,'\circ\cons\result.mha']);
            dvf_sum(:,:,:,kk) =dvfk;
        end
        if jj == 5;
            dvf_half = dvf_sum;
        end
    end
    % writing the "uncorrected" result
    MhaWrite(info,dvf_sum,[dir_out,'\circ\cons\circ_sum.mha']);
    % creating the mask (no point matching spurious regions
    % computing the correction
    v_sum = dvf_sum(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),:);
    circ_scale = [v_sum(:);circ_lam]\[zeros(length(v_sum(:)),1);circ_lam];
    MhaWrite(info,circ_scale*dvf_sum,[dir_out,'\circ\cons\circ_sum_scaled.mha']);
    %peak consistency
    MhaWrite(info,dvf_half,[dir_out,'\circ\cons\circ_half.mha']);
    [info,dvf_inex] = MhaRead([dir_out,'\inex\dvf\dvf_0_to_5.mha']);
    v_half = dvf_half(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),:);
    v_inex = dvf_inex(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),:);
    circ_scale_peak = [v_half(:);v_sum(:)]\...
        [v_inex(:);zeros(length(v_sum(:)),1)];
    MhaWrite(info,circ_scale_peak*dvf_sum,...
        [dir_out,'\circ\cons\circ_sum_scaled_peak.mha']);
    MhaWrite(info,circ_scale_peak*dvf_half,...
        [dir_out,'\circ\cons\circ_inex_scaled_peak.mha']);
    % additive residual method
    dvf_res = -0.1*dvf_sum;
    for jj = 1:10;
        [info,dvf] = MhaRead(names_dvfs{1,jj});
        dvf = dvf+dvf_res;
        name_dvf = names_dvfs{1,jj};
        name_dvf(end-4) = [];
        name_dvf = [name_dvf,'_rescor.mha'];
        MhaWrite(info,dvf,name_dvf);
    end
    save([dir_out,'\circ\cons\parameters.mat'],'circ_lam','circ_scale','circ_scale_peak');
end
%% Compare DVFs. We will treat the inhale-exhale DVF as ground truth,
% and compare against the circular DVF scaled peak to peak. No particular
% reason to consider inhale-exhale as more accurate, just for demonstration
%NOTE: We defined an arbitrary section of the lung at inhale as the
%region/structure of interest as this patient did not have particularly
%mobile tumours.
if do_compare
    mkdir([dir_out,'\analysis']);
    name_gt = [dir_out,'\inex\dvf\dvf_0_to_5.mha'];     %"ground truth"
    name_comp = [dir_out,'\circ\cons\circ_sum_scaled_peak.mha'];    %"comparison"
    %
    [info_dvf,dvf_gt] = MhaRead(name_gt);
    [info_dvf,dvf_comp] = MhaRead(name_comp);
    dvf_res = dvf_gt-dvf_comp;      %residual/ error dvf
    % build the mask (usually only worth doing analysis inside mask)
    mask = 0*dvf_res(:,:,:,1);
    mask(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3)) = 1;
    mask=logical(mask);
    %error as lengths
    vol_res = (dvf_res(:,:,:,1).^2+dvf_res(:,:,:,2).^2+dvf_res(:,:,:,3).^2).^0.5;
    %vol_res = vol_res.*mask;        %masking
    info_vol = rmfield(info_dvf,'ElementNumberOfChannels');
    MhaWrite(info_vol,vol_res,[dir_out,'\analysis\residual.mha']);
    % plotting
    vol_res_roi = vol_res(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3));
    dvf_res_roi_x = dvf_res(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),1);
    dvf_res_roi_y = dvf_res(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),2);
    dvf_res_roi_z = dvf_res(mask_mid(1)-mask_rad(1):mask_mid(1)+mask_rad(1),...
        mask_mid(2)-mask_rad(2):mask_mid(2)+mask_rad(2),...
        mask_mid(3)-mask_rad(3):mask_mid(3)+mask_rad(3),3);
    %%
    figure; set(gcf,'color','w');
    edges = [-15:1:30];
    h1 = histogram(dvf_res_roi_x(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0 0.4470 0.7410]);
    hold on
    h2 = histogram(dvf_res_roi_y(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0.8500 0.3250 0.0980]);
    h3 = histogram(dvf_res_roi_z(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0.9290 0.6940 0.1250]);
    h4 = histogram(vol_res_roi(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0 0 0],'FaceColor',[0 0 0]);
    legend('x residual','y residual','z residual','distance residual');
    axis tight
    % volume difference
    for jj = 1:3;
        temp = squeeze(dvf_gt(:,:,:,jj)).*mask;
        com_gt(jj) = sum(temp(:));
    end
    com_gt = com_gt./sum(mask(:));
    for jj = 1:3;
        temp = squeeze(dvf_comp(:,:,:,jj)).*mask;
        com_comp(jj) = sum(temp(:));
    end
    com_comp = com_comp./sum(mask(:));
    com_diff = com_gt-com_comp;
    %surface difference
    % requires a hull (edge of mask where mask is e.g. structure of
    % interest)
    k = ones([2,2,2]);
    hull = convn(mask,k,'same');
    hull(hull>=8) = 0;
    hull(hull>0) = 1;
    hull = logical(hull);
    % getting the location of each point on the hull
    x = info_dvf.PixelDimensions(1)*[1:info_dvf.Dimensions(1)];
    y = info_dvf.PixelDimensions(2)*[1:info_dvf.Dimensions(2)];
    z = info_dvf.PixelDimensions(3)*[1:info_dvf.Dimensions(3)];
    [X,Y,Z] = meshgrid(x,y,z);
    pc_hull(:,1) = X(hull);
    pc_hull(:,2) = Y(hull);
    pc_hull(:,3) = Z(hull);
    for jj = 1:3;
        dvf_hull = logical(zeros(size(dvf_gt)));
        dvf_hull(:,:,:,jj) = hull;
        pc_gt(:,jj) = pc_hull(:,jj)+dvf_gt(dvf_hull);
        pc_comp(:,jj) = pc_hull(:,jj)+dvf_comp(dvf_hull);
    end
    %we now have a point cloud pc of where the surface ended up
    % compute distances between point clouds
    % note: This potentially takes a long time. Consider skipping a
    % fraction of points (pc_skip>1)
    pc_n = length(pc_gt(:,1));
    pc_skip = 10;
    q = 1;
    %tic
    for jj = 1:pc_skip:pc_n
        %difference between every point
        pc_distj = pc_gt-ones(pc_n,1)*pc_comp(q,:);
        %distance between every point
        pc_distj = (pc_distj(:,1).^2+pc_distj(:,2).^2+pc_distj(:,3).^2).^0.5;   
        pc_mindist(q) = min(pc_distj);
        q = q+1;
    end
    %toc
    figure; set(gcf,'color','w');
    plot(sort(pc_mindist,'descend'));
    title('Distribution of Minimum Surface Distances');
    xlabel('Surface Distance Point Cloud (Largest to Smallest)');
    ylabel('Minimum Surface Distance (mm)');
end
%%
if do_compress;
    [dvfraw_info,dvfraw_vol] = MhaRead([dir_out,'\circ\cons\circ_sum_scaled_peak.mha']);
    rank = 5000;
    coeffs = zeros(rank,3);
    inds = zeros(rank,3);
    for jj = 1:3;
        %extract component
        P = dvfraw_vol(:,:,:,jj);
        %perform DCT
        Q = dct(P,[],1);
        R = dct(Q,[],2);
        S = dct(R,[],3);
        %sort to remove small components
        X = S(:);
        [~,ind] = sort(abs(X),'descend');
        % coeffs = 1;
        % while norm(X(ind(1:coeffs)))/norm(X) < 0.98
        %    coeffs = coeffs + 1;
        % end
        % cS(abs(S) < abs(X(ind(coeffs)))) = 0;
        cS = S;
        cS(ind(rank+1:end))=0;
        inds(:,jj) = ind(1:rank);
        coeffs(:,jj) = X(1:rank);
        cR = idct(cS,[],3);
        cQ = idct(cR,[],2);
        cP = idct(cQ,[],1);
        dvfcomp(:,:,:,jj) = cP;
    end
    %MhaWrite(dvfraw_info,dvfcomp,[dir_out,'\exref\dvf\dvf_1_to_3_compress.mha'])
    inds = single(inds);
    save([dir_out,'\analysis\compressed.mat'],'inds','coeffs');
    dvfloss = dvfraw_vol - dvfcomp;
    dvfloss_length = (squeeze(dvfloss(:,:,:,1).^2)+...
        squeeze(dvfloss(:,:,:,2).^2)+squeeze(dvfloss(:,:,:,3).^2)).^0.5;
    submm = round(100*(sum(dvfloss_length(:)<1)/length(dvfloss_length(:))));
    dvfgt_length = (squeeze(dvfraw_vol(:,:,:,1).^2)+...
        squeeze(dvfraw_vol(:,:,:,2).^2)+squeeze(dvfraw_vol(:,:,:,3).^2)).^0.5;
    submm = round(100*(sum(dvfloss_length(:)<1)/length(dvfloss_length(:))));
    %%
    figure; set(gcf,'color','w');
    edges = [0:1:20];
    h1 = histogram(dvfgt_length(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0 0.4470 0.7410]);
    hold on
    edges = [0:0.07:3];
    h2 = histogram(dvfloss_length(:),edges,'FaceAlpha',0.3,'LineWidth',2,...
        'EdgeAlpha',0.7,'EdgeColor',[0.8500 0.3250 0.0980]);
    legend('Original DVF Lengths','Compressed Representation Error Lengths',...
        'fontsize',20);
    axis tight
    title([num2str(submm),' percent of vectors sub 1mm accurate']);
    xlabel('millimetres','fontsize',20);ylabel('counts','fontsize',20);
    %%
    p = squeeze(P(:,70,:));
    figure; set(gcf,'color','w');
    imagesc(rot90(p,3));
    colormap(jet);
    c = colorbar;
    c.Label.String = 'DVF z Component (mm)';
    c.Location = 'southoutside';
    axis equal tight;
    %xticklabels({});yticklabels({});

    %%
    x = 0*p;
    xx = x;
    q = 1;
    for jj = 1:100
        for kk = 1:jj
            x(jj,kk) = X(ind(q));
            q = q+1;
            xx(kk,jj+1) = X(ind(q));
            q = q+1;
        end
    end
    x = x+xx;
    x = x(1:20,1:20);
%x = squeeze(S(1:100,2,1:100));
%x = squeeze(S(:,70,:));
    figure; set(gcf,'color','w');
    %imagesc(abs(x)); 
    %caxis([0 100]);
    %axis equal tight;
    surf(abs(x));
    view(130,35);
    caxis([0 1500]);
    colormap(jet);
    c = colorbar;
    c.Location = 'southoutside';
    c.Label.String = 'DCT Components';
end


