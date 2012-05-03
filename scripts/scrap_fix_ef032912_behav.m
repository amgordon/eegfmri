y= load('Acc1_retrieve_15_29Mar12out(1).mat');
y2 = load('AG3ONStudy29Mar11_9.mat');

d = dir('Acc*retrieve_9*.mat');

for i=1:length(d)
    z{i} =  load(d(i).name);
end

cd ../

idx.nonfix = y.theData.oldNew~=0;
itemsNonFix = y.theData.item(idx.nonfix);
oldItems = y2.respSelData.item;

y.theData.oldNew(idx.nonfix) = 2 - ismember(itemsNonFix, oldItems);

save('Acc1_retrieve_9_29Mar12out(0).mat', '-struct', 'y')

%%


for i=1:length(z)
    idx.repeatedTestStim = ismember(z{i}.theData.item, itemsNonFix)';
    z{i}.theData.REPEATED_JUNK = idx.repeatedTestStim;
    q = z{i};
    save(sprintf('Acc1_retrieve_9_29Mar12out(%g).mat',i), '-struct', 'q')
end