% [UNDER CONSTRUCTION]

% this script puts together expected taxonomies from isolated v4 amplicons
% doing it in matlab b/c R is too slow

% not sure how to subsample -- prob remove big stuff (silva will be hard to
% do that since names aren't mapped...), then make a set of intersecting
% amplicons from both? and create expected taxonomies from that?

clear; close all; 

v4p = readtable('pr2_v4amps.csv');
v4p = table2array(v4p(:,[2,7]));

v4s = readtable('silva_v4amps.csv');
v4s = table2array(v4s(:,[2,7]));

pr2db = fastaread('~/Documents/R/pr2_version_4.12.0_18S_dada2.fasta');
pr2db = struct2cell(pr2db)';

silvadb = fastaread('~/Documents/R/silva_nr_v138_train_set.fa');

v4p = unique(v4p(:,2));
v4s = unique(v4s(:,2));

v4p = v4p(1:100); % for testing the code
exptax = cell(length(v4p),9);
exptax(:,1) = v4p;
for i = 1:length(v4p)
    disp(['asv ',num2str(i), ' of ',num2str(length(v4p))]);
    s = v4p{i};
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
% writetable(exptax, 'pr2_v4amp_exptax.csv');
