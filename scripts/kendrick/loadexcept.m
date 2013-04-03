function loadexcept(file,vars)

% function loadexcept(file,vars)
%
% <file> is a string referring to a .mat file
% <vars> is a variable name or a cell vector of variable names to NOT load
%
% load <file> into the base workspace, excluding variables named by <vars>.
%
% example:
% x = 1; y = 2;
% save('temp.mat','x','y');
% clear x y;
% loadexcept('temp.mat','x')
% ~exist('x','var')
% exist('y','var')

% input
if ~iscell(vars)
  vars = {vars};
end

% figure out variable names
varlist = whos('-file',file);
varlist = cat(2,{varlist.name});

% exclude the ones we don't want
ok = cellfun(@(x) ~ismember(x,vars),varlist);
varlist = varlist(ok);

% load in the data
data = load(file,varlist{:});

% assign to caller's workspace
for p=1:length(varlist)
  assignin('base',varlist{p},data.(varlist{p}));
end
