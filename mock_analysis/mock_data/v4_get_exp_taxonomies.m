% this script puts together expected taxonomies from isolated v4 amplicons
% doing it in matlab b/c R is too slow

% I figured out how to subsample, randomly subsampling from the unique ASVs
% of each database. 

clear; close all; 

v4p = readtable('pr2_v4amps.csv');
v4p = table2array(v4p(:,[2,7]));

v4s = readtable('silva_v4amps.csv');
v4s = table2array(v4s(:,[2,7]));

pr2db = fastaread('~/Documents/R/pr2_version_4.12.0_18S_dada2.fasta');
pr2db = struct2cell(pr2db)';

silvadb = fastaread('~/Documents/R/silva_nr_v138_train_set.fa');
silvadb = struct2cell(silvadb)';
xx = count(v4s,'|');
x2 = count(silvadb(:,1), ';');
for i = 1:5
    idx = find(x2 == i);
    addme = repmat('NA;',1,6-i);
    silvadb(idx,1) = cellfun(@(x) cat(2,x,addme),silvadb(idx,1),'UniformOutput',false);
end

v4p = unique(v4p(:,2));
v4s = unique(v4s(:,2));

nn = 2500;
% randomly select 2500 sequences from both silva and pr2 amplicons
is = randperm(length(v4s),nn);
ip = randperm(length(v4p),nn);

v4s = v4s(is);
v4p = v4p(ip);

us = unique([v4s;v4p],'stable');
pexptax = cell(length(us),9);
sexptax = cell(length(us),max(x2)+1);
pexptax(:,1) = us;
sexptax(:,1) = us;
pr2flag = cell(nn,1); % for sequences without exact matches in pr2
silflag = cell(nn,1); % for sequences without exact matches in pr2
for i = 1:length(us)
    disp(['asv ',num2str(i), ' of ',num2str(length(us))]);
    s = us{i};
    % pr2 check:
    inpr2 = contains(pr2db(:,2), s);
    idx = find(inpr2 == 1);
    if length(idx) == 1
        et = split(pr2db(idx,1),';')';
        pexptax(i,2:end) = et(cellfun(@isempty,et) == 0);
    elseif length(idx) > 1
        et = split(pr2db(idx,1),';');
        et = et(:,sum(cellfun(@isempty,et)) < size(et,1));
        j = 1;
        while length(unique(et(:,j))) == 1
            j = j+1;
            if j > 8
                break
            end
        end
        j = j-1; % ^that will end when j is 9, or when j is 1 col beyond where the names agree
        et = table2cell(unique(cell2table(et(:,1:j)),'rows','stable'));
        pexptax(i,2:length(et)+1) = et;
    elseif isempty(idx) 
        pr2flag{i} = s;
    end
    
    % silva check:
    insil = contains(silvadb(:,2), s);
    idx = find(insil == 1);
    if length(idx) == 1
        et = split(silvadb(idx,1),';')';
        et = et(cellfun(@isempty,et) == 0);
        sexptax(i,2:1+length(et)) = et;
    elseif length(idx) > 1
        et = split(silvadb(idx,1),';');
        et = et(:,sum(cellfun(@isempty,et)) < size(et,1));
        j = 1;
        while length(unique(et(:,j))) == 1
            j = j+1;
            if j > size(et,2)
                break
            end
        end
        j = j-1; % ^that will end when j is 7, or when j is 1 col beyond where the names agree
        et = table2cell(unique(cell2table(et(:,1:j)),'rows','stable'));
        sexptax(i,2:length(et)+1) = et;
    elseif isempty(idx) 
        silflag{i} = s;
    end
end

pexptax(cellfun(@isempty, pexptax)) = {NaN};
sexptax(cellfun(@isempty, sexptax)) = {NaN};
pexptax = cell2table(pexptax);
sexptax = cell2table(sexptax);
writetable(pexptax, 'pr2_v4amp_exptax.csv');
writetable(sexptax, 'silva_v4amp_exptax.csv');

pr2flag = pr2flag(~cellfun(@isempty, pr2flag));
silflag = silflag(~cellfun(@isempty, silflag));
writetable(cell2table(pr2flag), 'v4_asvs_no_exact_pr2.csv');
writetable(cell2table(silflag), 'v4_asvs_no_exact_silva.csv');

