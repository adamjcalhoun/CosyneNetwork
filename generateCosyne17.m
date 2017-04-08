% analyzeCosyne2017

files{1} = 'Data/COSYNE04_names.mat';
files{2} = 'Data/COSYNE05_names.mat';
files{3} = 'Data/COSYNE06_names.mat';
files{4} = 'Data/COSYNE07_names.mat';
files{5} = 'Data/COSYNE08_names.mat';
files{6} = 'Data/COSYNE09_names.mat';
files{7} = 'Data/COSYNE10_names.mat';
files{8} = 'Data/COSYNE11_names.mat';
files{9} = 'Data/COSYNE12_names.mat';
files{10} = 'Data/COSYNE13_names.mat';
files{11} = 'Data/COSYNE14_names.mat';
files{12} = 'Data/COSYNE15_names.mat';
files{13} =  'Data/COSYNE16_names.mat';
files{14} =  'Data/COSYNE17_names.mat';

cosyneData.lastnames = [];
cosyneData.firstnames = [];
cosyneData.posterID = [];
for i=1:length(files)
    load(files{i});

    cosyneData.lastnames = [cosyneData.lastnames strtrim(lastnames)];
    cosyneData.firstnames = [cosyneData.firstnames strtrim(firstnames)];
    cosyneData.posterID = [cosyneData.posterID posterID];
end

save 2017/cosyneData.mat cosyneData

%% manual corrections
% do I need to change Larry, too?

load 2017/cosyneData.mat

for i=1:length(cosyneData.lastnames)
    if length(cosyneData.lastnames{i}) > 5
        if strcmp(cosyneData.lastnames{i}(end-3:end),'ding') && strcmp(cosyneData.firstnames{i}(1),'K') && strcmp(cosyneData.lastnames{i}(1),'K')
            cosyneData.lastnames{i} = 'Kording';
        elseif strcmp(cosyneData.lastnames{i},'Abbott') && strcmp(cosyneData.firstnames{i},'Laurence')
%             cosyneData.firstnames{i} = 'Aurence';
        end
    end
end

save 2017/cosyneData.mat cosyneData

%% now assign IDs to each author and poster
load 2017/cosyneData.mat

% - get unique author IDs
names = cosyneData.lastnames;
for i=1:length(names)
    names{i} = [cosyneData.firstnames{i}(1) names{i}];
end

cosyneData.authorID = zeros(1,length(cosyneData.lastnames));
[uniquenames,~,nameInd] = unique(names);
for i=1:length(uniquenames)
    newnames = find(nameInd == i);
    cosyneData.authorID(newnames) = i;
end

save 2017/cosyneData.mat cosyneData

%% - then assign each poster a unique integer ID
load 2017/cosyneData.mat
posters = cell2mat(cosyneData.posterID);
[uniquePosters,posterInd] = unique(posters);
cosyneData.posterHash = [uniquePosters',(1:length(uniquePosters))'];

save 2017/cosyneData.mat cosyneData

%% - then go through each poster and make a poster x author matrix with a 1
% for each assignment
load 2017/cosyneData.mat

posterMatrix = zeros(size(cosyneData.posterHash,1),length(unique(cosyneData.authorID)));
for j=1:length(cosyneData.posterID)
    for k=1:length(cosyneData.posterID{j})
        posterMatrix(cosyneData.posterHash(cosyneData.posterHash(:,1) == cosyneData.posterID{j}(k),2),cosyneData.authorID(j)) = 1;
    end
end

adjMatrix = zeros(length(unique(cosyneData.authorID)));
for i=1:length(unique(cosyneData.authorID))
    posters = find(posterMatrix(:,i) ~= 0);
    for j=1:length(posters)
        adjMatrix(i,posterMatrix(posters(j),:) ~= 0) = adjMatrix(i,posterMatrix(posters(j),:) ~= 0) + 1;
        adjMatrix(i,i) = 0;
    end
end

save 2017/cosyneFinalData.mat cosyneData adjMatrix posterMatrix

%% Or we could make the adjacency matrix only for 2017!
load 2017/cosyneData.mat

posterMatrix = zeros(size(cosyneData.posterHash,1),length(unique(cosyneData.authorID)));
for j=1:length(cosyneData.posterID)
    for k=1:length(cosyneData.posterID{j})
%         if round(rem(cosyneData.posterID{j}(k),1)*10000) > 2009
        if round(rem(cosyneData.posterID{j}(k),1)*10000) == 2017
            posterMatrix(cosyneData.posterHash(cosyneData.posterHash(:,1) == cosyneData.posterID{j}(k),2),cosyneData.authorID(j)) = 1;
        end
    end
end

adjMatrix = zeros(length(unique(cosyneData.authorID)));
for i=1:length(unique(cosyneData.authorID))
    posters = find(posterMatrix(:,i) ~= 0);
    for j=1:length(posters)
        adjMatrix(i,posterMatrix(posters(j),:) ~= 0) = adjMatrix(i,posterMatrix(posters(j),:) ~= 0) + 1;
        adjMatrix(i,i) = 0;
    end
end

save 2017/cosyne2017FinalData.mat cosyneData adjMatrix posterMatrix

%% Find the people with the most posters/etc
% load 2017/cosyne2017FinalData.mat
load 2017/cosyneFinalData.mat

numPosters = sum(posterMatrix);
[sorted,index] = sort(numPosters,2,'descend');
for i=1:30
    x = find([cosyneData.authorID]==index(i));
    display(['(' num2str(sorted(i)) '): ' cosyneData.firstnames{x(1)}(1) '. ' cosyneData.lastnames{x(1)} ' [' num2str(cosyneData.authorID(x(1))) ']'])
    display([cosyneData.firstnames{x(1)}(1) '. ' cosyneData.lastnames{x(1)}])
    display([num2str(sorted(i))])
    data(i) = sorted(i);
    data2{i} = [cosyneData.firstnames{x(1)}(1) '. ' cosyneData.lastnames{x(1)}];
end
data = data';

%% Write out cosyne adjacency matrix for Gephi
% 2017 data only
load 2017/cosyne2017FinalData.mat
fid = fopen('2017/cosyneFinalData.csv','w');

[uniqueID,ind] = unique([cosyneData.authorID]);
% for i=1:length(unique(cosyneData.authorID))
for i=1:length(uniqueID)
    posters = find(posterMatrix(:,i) ~= 0);
    author = [cosyneData.firstnames{ind(i)}(1) '.' cosyneData.lastnames{ind(i)}];
    for j=1:length(posters)
        otherAuthors = find(posterMatrix(posters(j),:) ~= 0);
        for k=1:length(otherAuthors)
            if (i ~= otherAuthors(k))
                othername = [cosyneData.firstnames{ind(otherAuthors(k))}(1) '.' cosyneData.lastnames{ind(otherAuthors(k))}];
                fprintf(fid,'%s;%s\n',author,othername);
            end
        end
    end
end

fclose(fid);
