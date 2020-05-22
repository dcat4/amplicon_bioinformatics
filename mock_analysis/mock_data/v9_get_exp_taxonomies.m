% this script puts together expected taxonomies from isolated v9 amplicons
% doing it in matlab b/c R is too slow

% seems to work, takes ~4 hr to run locally on 

clear; close all; 

v9p = readtable('pr2_v9amps_bothPrimers.csv');
v9p = table2array(v9p(:,[2,7]));

v9p = unique(v9p(:,2));
pr2db = fastaread('~/Documents/R/pr2_version_4.12.0_18S_dada2.fasta');
pr2db = struct2cell(pr2db)';

% v9p = v9p(1:100); % for testing the code
exptax = cell(length(v9p),9);
exptax(:,1) = v9p;
for i = 1:length(v9p)
    disp(['asv ',num2str(i), ' of ',num2str(length(v9p))]);
    s = v9p{i};
    inpr2 = contains(pr2db(:,2), s);
    idx = find(inpr2 == 1);
    if length(idx) == 1
        et = split(pr2db(idx,1),';')';
        exptax(i,2:end) = et(cellfun(@isempty,et) == 0);
    elseif length(idx) > 1
        et = split(pr2db(idx,1),';');
        et = et(:,sum(cellfun(@isempty,et)) < size(et,1));
        bloop = et(:,1);
        j = 1;
        while length(unique(et(:,j))) == 1
            j = j+1;
            if j > 8
                break
            end
        end
        j = j-1; % ^that will end when j is 9, or when j is 1 col beyond where the names agree
        et = table2cell(unique(cell2table(et(:,1:j)),'rows','stable'));
        exptax(i,2:length(et)+1) = et;
    end
end

exptax(cellfun(@isempty, exptax)) = {NaN};
exptax = cell2table(exptax);
writetable(exptax, 'pr2_v9amp_exptax_bothPrimers.csv');


