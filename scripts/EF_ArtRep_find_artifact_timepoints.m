function [zscoreA_cell, delta_cell] = AG1_ArtRep_find_artifact_timepoints(subpar)

%% runs functions from Paul Mazaika's ArtRepair5 toolbox to determine artifacts in the fMRI timeseries data
%% prepared by Vincent Bell

% -----------------------
% Initialize, begin loop
% -----------------------

% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = AG1Params(subpar);
end

%par.Tasks = 'RS';

if ~exist(par.artrepdir)
    mkdir (par.artrepdir)
end

cd(par.artrepdir);

%nRuns = par.numscans;
scans_to_include = par.scans_to_include;

pfig = [];


      Percent_thresh = 1.5; 

      %z_thresh = zthresh %2;  % Currently not used for default.

      mv_thresh = 0.5;  % try 0.3 for subjects with intervals with low noise
 
      
%%K= current run
for k = 1:length(par.scans_to_include)
    x = scans_to_include(k);
    

     
    num_sess = 1;%d.nsess(subj); % spm_input('How many sessions?',1,'n',1,1);
    global_type_flag = 2; %spm_input('Which global mean to use?', 1, 'm', ...
    realignfile = 1;   % Default case is like artdetect4.

    M = []; P = [];
    for i = 1:num_sess
        P{i} = par.ascanfiles{x}; %d.imagefilenames{subj}{i} %spm_select(Inf,'image',['Select data images for session'  num2str(i) ':']);
        if realignfile == 1

            scan_dir = fullfile(par.funcdir, ['scan' prepend(num2str(x))]);
            mvmt_file = dir([scan_dir, '/rp_ascan*']);
            Mh = load(fullfile(scan_dir, mvmt_file.name));  

            M{i} = Mh;
            
        end
    end
    repair1_flag = 1;%spm_input('Always repair 1st scan of each session?', '+1', 'y/n', [1 0], 1);
    GoRepair = 0;       % GUI repair


if global_type_flag==3
    maskY = spm_read_vols(spm_vol(mask));

    maskcount = sum(sum(sum(maskY)));  %  Number of voxels in mask.
    voxelcount = prod(size(maskY));    %  Number of voxels in 3D volume.
end
if global_type_flag == 4   %  Automask option
    disp('Generated mask image is written to file ArtifactMask.img.')
    Pnames = P{1};
    Automask = art_automask(Pnames(1,:),-1,1);
    maskcount = sum(sum(sum(Automask)));  %  Number of voxels in mask.
    voxelcount = prod(size(Automask));    %  Number of voxels in 3D volume.
end
spm_input('!DeleteInputObj');

P = char(P);

    mv_data = [];
    for i = 1:length(M)
        mv_data = vertcat(mv_data,M{i});
    end


fprintf('%-4s: ','Mapping files...')                                  
VY     = spm_vol(P);
fprintf('%3s\n','...done')                                          


temp = any(diff(cat(1,VY.dim),1,1),1);
if length(temp) == 4       % SPM2 case
    if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0])
        error('images do not all have the same dimensions')
    end
elseif length(temp) == 3   % SPM5 nifti case
    if ~isempty(find(diff(cat(1,VY.dim)) ~= 0 ))   
	    error('images do not all have the same dimensions (SPM5)')
    end
end

nscans = size(P,1);

g      = zeros(nscans,1);

fprintf('%-4s: %3s','Calculating globals...',' ')
if global_type_flag==1  % regular mean
    for i  = 1:nscans  
	    g(i) = spm_global(VY(i));
    end
elseif global_type_flag==2  % every voxel
    for i = 1:nscans
        g(i) = mean(mean(mean(spm_read_vols(VY(i)))));
    end
elseif global_type_flag == 3 % user masked mean
     Y = spm_read_vols(VY(1));
    %voxelcount = prod(size(Y));
    %vinv = inv(VY(1).mat);
    %[dummy, idx_to_mask] = intersect(XYZmm', maskXYZmm', 'rows');
    %maskcount = length(idx_to_mask
    for i = 1:nscans
        Y = spm_read_vols(VY(i)); 
        Y = Y.*maskY;
        %Y(idx_to_mask) = [];   
        g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        %g(i) = mean(Y(idx_to_mask));
    end
else   %  global_type_flag == 4  %  auto mask
    for i = 1:nscans
        Y = spm_read_vols(VY(i));
        Y = Y.*Automask;
        if realignfile == 0
            output = art_centroid(Y);
            centroiddata(i,1:3) = output(2:4);
            g(i) = output(1)*voxelcount/maskcount;
        else     % realignfile == 1
            g(i) = mean(mean(mean(Y)))*voxelcount/maskcount;
        end
    end
    % If computing approximate translation alignment on the fly...
    %   centroid was computed in voxels
    %   voxel size is VY(1).mat(1,1), (2,2), (3,3).
    %   calculate distance from mean as our realignment estimate
    %   set rotation parameters to zero.
    if realignfile == 0    % change to error values instead of means.
        centroidmean = mean(centroiddata,1);
        for i = 1:nscans
            mv0data(i,:) = - centroiddata(i,:) + centroidmean;
        end
        % THIS MAY FLIP L-R  (x translation)
        mv_data(1:nscans,1) = mv0data(1:nscans,1)*VY(1).mat(1,1);
        mv_data(1:nscans,2) = mv0data(1:nscans,2)*VY(1).mat(2,2);
        mv_data(1:nscans,3) = mv0data(1:nscans,3)*VY(1).mat(3,3);
        mv_data(1:nscans,4:6) = 0;
    end
end

% Convert rotation movement to degrees
mv_data(:,4:6)= mv_data(:,4:6)*180/pi; 
    
    
fprintf('%s%3s\n','...done\n')
if global_type_flag==3
    fprintf('\n%g voxels were in the user mask.\n', maskcount)
end
if global_type_flag==4
    fprintf('\n%g voxels were in the auto generated mask.\n', maskcount)
end

% ------------------------
% Compute default out indices by z-score, or by Percent-level is std is small.
% ------------------------ 
%  Consider values > Percent_thresh as outliers (instead of z_thresh*gsigma) if std is small.
    gsigma = std(g);
    gmean = mean(g);
    pctmap = 100*gsigma/gmean;
    mincount = Percent_thresh*gmean/100;
    %z_thresh = max( z_thresh, mincount/gsigma );
    z_thresh = mincount/gsigma;        % Default value is PercentThresh.
    z_thresh = 0.1*round(z_thresh*10); % Round to nearest 0.1 Z-score value
    zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
    glout_idx = (find(abs(zscoreA) > z_thresh))';

% ------------------------
% Compute default out indices from rapid movement
% ------------------------ 
%   % Rotation measure assumes voxel is 65 mm from origin of rotation.
    if realignfile == 1 | realignfile == 0
        delta = zeros(nscans,1);  % Mean square displacement in two scans
        for i = 2:nscans
            delta(i,1) = (mv_data(i-1,1) - mv_data(i,1))^2 +...
                    (mv_data(i-1,2) - mv_data(i,2))^2 +...
                    (mv_data(i-1,3) - mv_data(i,3))^2 +...
                    1.28*(mv_data(i-1,4) - mv_data(i,4))^2 +...
                    1.28*(mv_data(i-1,5) - mv_data(i,5))^2 +...
                    1.28*(mv_data(i-1,6) - mv_data(i,6))^2;
            delta(i,1) = sqrt(delta(i,1));
        end
    end
    
     % Also name the scans before the big motions (v2.2 fix)
    deltaw = zeros(nscans,1);
    for i = 1:nscans-1
        deltaw(i) = max(delta(i), delta(i+1));
    end
    delta(1:nscans-1,1) = deltaw(1:nscans-1,1);
    mvout_idx = find(delta > mv_thresh)';
    
    % Total repair list
    out_idx = unique([mvout_idx glout_idx]);
    if repair1_flag == 1
        out_idx = unique([ 1 out_idx]);
    end
    % Initial deweight list before margins
    outdw_idx = out_idx; 
    % Initial clip list without removing large displacements
    clipout_idx = out_idx;
     

% % -----------------------
% % Draw initial figure
% % -----------------------
% 
% 
% figure('Units', 'normalized', 'Position', [0.2 0.2 0.6 0.7]);
% rng = max(g) - min(g);   % was range(g);
% pfig = gcf;
% 
% subplot(5,1,1);
% plot(g);
% %xlabel(['artifact index list [' int2str(out_idx') ']'], 'FontSize', 8, 'Color','r');
% %ylabel(['Range = ' num2str(rng)], 'FontSize', 8);
% ylabel('Global Avg. Signal');
% xlabel('Red vertical lines are to depair. Green vertical lines are to deweight.');
% title('ArtifactRepair GUI to repair outliers and identify scans to deweight');
% %if ( global_type_flag == 1 ) title('Global Mean - Regular SPM'); end
% %if ( global_type_flag == 2 ) title('Global Mean - Every Voxel'); end
% %if ( global_type_flag == 3 ) title('Global Mean - User Defined Mask'); end
% %if ( global_type_flag == 4 ) title('Global Mean - Generated ArtifactMask'); end
% 
% % Add vertical exclusion lines to the global intensity plot
% axes_lim = get(gca, 'YLim');
% axes_height = [axes_lim(1) axes_lim(2)];
% for i = 1:length(outdw_idx)   % Scans to be Deweighted
%     line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'g');
% end
% if GoRepair == 2
%     for i = 1:length(outdw_idx)   % Scans to be Deweighted
%         line((outdw_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
%     end
% end
% subplot(5,1,2);
% %thresh_axes = gca;
% %set(gca, 'Tag', 'threshaxes');
% zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
% plot(abs(zscoreA));
% ylabel('Std away from mean');
% xlabel('Scan Number  -  horizontal axis for all plots');
% 
% thresh_x = 1:nscans;
% thresh_y = z_thresh*ones(1,nscans);
% line(thresh_x, thresh_y, 'Color', 'r');
% 
% %  Mark global intensity outlier images with vertical lines
% axes_lim = get(gca, 'YLim');
% axes_height = [axes_lim(1) axes_lim(2)];
% for i = 1:length(glout_idx)
%     line((glout_idx(i)*ones(1, 2)), axes_height, 'Color', 'r');
% end
% 
% if realignfile == 1
% 	subplot(5,1,3);
%     xa = [ 1:nscans];
% 	plot(xa,mv_data(:,1),'b-',xa,mv_data(:,2),'g-',xa,mv_data(:,3),'r-',...
%        xa,mv_data(:,4),'r--',xa,mv_data(:,5),'b--',xa,mv_data(:,6),'g--');
%     %plot(,'--');
% 	ylabel('ReAlignment');
% 	xlabel('Translation (mm) solid lines, Rotation (deg) dashed lines');
% 	legend('x mvmt', 'y mvmt', 'z mvmt','pitch','roll','yaw',0);
% 	h = gca;
% 	set(h,'Ygrid','on');
% elseif realignfile == 0
%     subplot(5,1,3);
% 	plot(mv0data(:,1:3));
% 	ylabel('Alignment (voxels)');
% 	xlabel('Scans. VERY APPROXIMATE EARLY-LOOK translation in voxels.');
% 	legend('x mvmt', 'y mvmt', 'z mvmt',0);
% 	h = gca;
% 	set(h,'Ygrid','on');
% end 
% 
% subplot(5,1,4);   % Rapid movement plot
% plot(delta);
% ylabel('Motion (mm/TR)');
% xlabel('Scan to scan movement (~mm). Rotation assumes 65 mm from origin');
% y_lim = get(gca, 'YLim');
% legend('Fast motion',0);
% h = gca;
% set(h,'Ygrid','on');
% 
% thresh_x = 1:nscans;
% thresh_y = mv_thresh*ones(1,nscans);
% line(thresh_x, thresh_y, 'Color', 'r');
%    
% % Mark all movement outliers with vertical lines
% subplot(5,1,4)
% axes_lim = get(gca, 'YLim');
% axes_height = [axes_lim(1) axes_lim(2)];
% for i = 1:length(mvout_idx)
%     line((mvout_idx(i)*ones(1,2)), axes_height, 'Color', 'r');
% end
% 
% %keyboard;
% h_rangetext = uicontrol(gcf, 'Units', 'characters', 'Position', [10 10 18 2],...
%         'String', 'StdDev of data is: ', 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_rangenum = uicontrol(gcf, 'Units', 'characters', 'Position', [29 10 10 2], ...
%         'String', num2str(gsigma), 'Style', 'text', ...
%         'HorizontalAlignment', 'left',...
%         'Tag', 'rangenum',...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_threshtext = uicontrol(gcf, 'Units', 'characters', 'Position', [25 8 16 2],...
%         'String', 'Current threshold (std devs):', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_threshnum = uicontrol(gcf, 'Units', 'characters', 'Position', [44 8 10 2],...
%         'String', num2str(z_thresh), 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'threshnum');
% h_threshmvtext = uicontrol(gcf, 'Units', 'characters', 'Position', [106 8 18 2],...
%         'String', 'Motion threshold  (mm / TR):', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_threshnummv = uicontrol(gcf, 'Units', 'characters', 'Position', [126 8 10 2],...
%         'String', num2str(mv_thresh), 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'threshnummv');
% h_threshtextpct = uicontrol(gcf, 'Units', 'characters', 'Position', [66 8 16 2],...
%         'String', 'Current threshold (% of mean):', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_threshnumpct = uicontrol(gcf, 'Units', 'characters', 'Position', [86 8 10 2],...
%         'String', num2str(z_thresh*pctmap), 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'threshnumpct');
% h_deweightlist = uicontrol(gcf, 'Units', 'characters', 'Position', [150 6 1 1 ],...
%         'String', int2str(outdw_idx), 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'deweightlist');
% h_clipmvmtlist = uicontrol(gcf, 'Units', 'characters', 'Position', [152 6 1 1 ],...
%         'String', int2str(clipout_idx), 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'clipmvmtlist');
% h_indextext = uicontrol(gcf, 'Units', 'characters', 'Position', [10 3 15 2],...
%         'String', 'Outlier indices: ', 'Style', 'text', ...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8], ...
%         'ForegroundColor', 'r');
% h_indexedit = uicontrol(gcf, 'Units', 'characters', 'Position', [25 3.25 40 2],...
%         'String', int2str(out_idx), 'Style', 'edit', ...
%         'HorizontalAlignment', 'left', ...
%         'Callback', 'art_outlieredit',...
%         'BackgroundColor', [0.8 0.8 0.8],...
%         'Tag', 'indexedit');
% h_indexinst = uicontrol(gcf, 'Units', 'characters', 'Position', [66 3 40 2],...
%         'String', '[Hit return to update after editing]', 'Style', 'text',...
%         'HorizontalAlignment', 'left', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_deweighttext = uicontrol(gcf, 'Units', 'characters', 'Position', [115 1 21 2],...
%         'String', 'Click to add margins for deweighting', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% h_clipmvmttext = uicontrol(gcf, 'Units', 'characters', 'Position', [94 1 21 2],...
%         'String', 'Mark > 3mm movments for repair', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% if realignfile == 1
%    h_repairtext = uicontrol(gcf, 'Units', 'characters', 'Position', [137 1 17 2],...
%         'String', 'Writes repaired volumes', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% else  % realignfile == 0
%    h_repairtext = uicontrol(gcf, 'Units', 'characters', 'Position', [137 1 17 2],...
%         'String', 'WARNING!! DATA NOT REALIGNED', 'Style', 'text', ...
%         'HorizontalAlignment', 'center', ...
%         'BackgroundColor', [0.8 0.8 0.8]);
% end
% h_addmargin = uicontrol(gcf, 'Units', 'characters', 'Position', [120 3.25 10 2],...
%         'String', 'Margin', 'Style', 'pushbutton', ...
%         'Tooltipstring', 'Adds margins to deweight in estimation', ...
%         'BackgroundColor', [ 0.7 0.9 0.7], 'ForegroundColor','k',...
%         'Callback', 'art_addmargin');
% h_clipmvmt = uicontrol(gcf, 'Units', 'characters', 'Position', [100 3.25 10 2],...
%         'String', 'Clip', 'Style', 'pushbutton', ...
%         'Tooltipstring', 'Marks displacements > 3 mm for repair', ...
%         'BackgroundColor', [ 0.7 0.7 0.8 ], 'ForegroundColor','k',...
%         'Callback', 'art_clipmvmt');
% h_repair = uicontrol(gcf, 'Units', 'characters', 'Position', [140 3.25 10 2],...
%         'String', 'REPAIR', 'Style', 'pushbutton', ...
%         'Tooltipstring', 'Writes repaired images', ...
%         'BackgroundColor', [ 0.9 0.7 0.7], 'ForegroundColor','r',...
%         'Callback', 'art_repair');
% h_up = uicontrol(gcf, 'Units', 'characters', 'Position', [10 8 10 2],...
%         'String', 'Up', 'Style', 'pushbutton', ...
%         'TooltipString', 'Raise threshold for outliers', ...
%         'Callback', 'art_threshup');
% h_down = uicontrol(gcf, 'Units', 'characters', 'Position', [10 6 10 2],...
%         'String', 'Down', 'Style', 'pushbutton', ...
%         'TooltipString', 'Lower threshold for outliers', ...
%         'Callback', 'art_threshdown');
% 
% %guidata(gcf, g);
% guidata(gcf,[g delta mv_data]);
% setappdata(h_repair,'data',P);
% setappdata(h_addmargin,'data2',repair1_flag);
% setappdata(h_clipmvmt,'data3',repair1_flag);
% %setappdata(h_repair,'data3',GoRepair);
% 
% % For GUI or RepairAlone, add deweighting margins to top plot
% % Don't apply margins when motion adjustment was used.(v2.2)
% if GoRepair == 0 | GoRepair == 1
%     art_clipmvmt;
%     art_addmargin;
% end
% 
% if GoRepair == 1 | GoRepair == 2
%     ImageFile = P(1,:);
%     [subjectpath, sess] = fileparts(ImageFile);
%     [uppath1, subname ] = fileparts(subjectpath);
%     [uppath2, up1name ] = fileparts(uppath1);
%     [uppath3, up2name ] = fileparts(uppath2);
%     [uppath4, up3name ] = fileparts(uppath3);
%     titlename = [ up3name, ' / ', up2name,' / ', up1name,' / ', subname ];
%     gcf; subplot(5,1,1);
%     title(titlename);
%     figname = ['artglobal', subname, '.jpg'];
%     try
%         saveas(gcf,figname);
%     catch
%         disp('Warning from art_global: Could not print figure');
%     end
%     art_repair(P);

zscoreA_cell{k}=zscoreA;
delta_cell{k}=delta;
allSigArt{k} = (abs(zscoreA_cell{k}) < par.art.sigThresh);
allMotArt{k} = (delta_cell{k} < par.art.motThresh);

allArt{k} = (allSigArt{k} .* allMotArt{k});
end

sum(find(vertcat(allSigArt{:}))==0)
sum(find(vertcat(allMotArt{:}))==0)

filename=['art_global_modified_' par.substr '.mat'];
save (filename, 'zscoreA_cell', 'delta_cell', 'allSigArt', 'allMotArt', 'allArt');





