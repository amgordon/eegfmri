function filenames = AG_mvpa_imgs_list_maker (S)



yEnc = load(fullfile(univar_dir.Enc, 'SPM'));

yRet = load(fullfile(univar_dir.Ret, 'SPM'));

yLoc = load(fullfile(univar_dir.Loc, 'SPM'));

yRS = load(fullfile(univar_dir.RS, 'SPM'));

v1 = cellstr(yEnc.SPM.xY.P);
v2 = cellstr(yRet.SPM.xY.P);
v3 = cellstr(yLoc.SPM.xY.P);
v4 = cellstr(yRS.SPM.xY.P);

filenames = vertcat(v1, v2);

%filenames = v4;