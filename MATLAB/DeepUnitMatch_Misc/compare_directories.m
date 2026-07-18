% Compare contents of two directories and flag any differences
dir1 = 'C:\Users\EnnyB\Dropbox\DeepUnitMatch_Data\ALL_DATA_nomatchtables\';
dir2 = '\\znas.cortexlab.net\Lab\Share\Suyash\ALL_DATA\';
outputFile = fullfile(fileparts(mfilename('fullpath')), 'directory_comparison.xlsx');

fprintf('Comparing:\n  Dir1: %s\n  Dir2: %s\n\n', dir1, dir2);

entries1 = getEntries(dir1);
entries2 = getEntries(dir2);

only_in_1 = setdiff(entries1, entries2);
only_in_2 = setdiff(entries2, entries1);

if isempty(only_in_1) && isempty(only_in_2)
    fprintf('Directories are identical.\n');
else
    if ~isempty(only_in_1)
        fprintf('--- Only in Dir1 (%d items) ---\n', numel(only_in_1));
        for i = 1:numel(only_in_1)
            fprintf('  %s\n', only_in_1{i});
        end
        fprintf('\n');
    end
    if ~isempty(only_in_2)
        fprintf('--- Only in Dir2 (%d items) ---\n', numel(only_in_2));
        for i = 1:numel(only_in_2)
            fprintf('  %s\n', only_in_2{i});
        end
        fprintf('\n');
    end
end

% --- Detect potential renames: Dir1 uses long _PyKS_output names, Dir2 uses shorter names ---
% Pattern: strip _PyKS_output suffix, then take everything from __ onward
% e.g. _2019-09-30_1__2019-09-30_run1_g1_imec0_PyKS_output -> _2019-09-30_run1_g1_imec0
potentialRenames = cell(numel(only_in_1), 5);
nMatches = 0;
for i = 1:numel(only_in_1)
    p1 = only_in_1{i};
    parts = strsplit(p1, '\');
    leaf = parts{end};

    derived = regexprep(leaf, '_pyks_output$', '', 'ignorecase');
    dbl = strfind(derived, '__');
    if ~isempty(dbl)
        derived = derived(dbl(1)+1:end);  % keep from the second underscore onward
    end

    if isempty(derived) || strcmp(derived, leaf)
        continue;
    end

    candidate = strjoin([parts(1:end-1), {derived}], '\');
    if ~ismember(candidate, only_in_2)
        continue;
    end

    prefix1 = [p1 '\'];
    prefix2 = [candidate '\'];
    n1 = sum(strncmp(entries1, prefix1, length(prefix1)));
    n2 = sum(strncmp(entries2, prefix2, length(prefix2)));
    nMatches = nMatches + 1;
    potentialRenames(nMatches, :) = {p1, candidate, n1, n2, n1 == n2};
end
potentialRenames = potentialRenames(1:nMatches, :);

% Build exclusion masks: remove matched folders AND all their children
if nMatches > 0
    matchedFolders1 = potentialRenames(:, 1);
    matchedFolders2 = potentialRenames(:, 2);

    exclude1 = false(size(only_in_1));
    for j = 1:numel(matchedFolders1)
        mf = matchedFolders1{j};
        exclude1 = exclude1 | strcmp(only_in_1, mf) | strncmp(only_in_1, [mf '\'], length(mf)+1);
    end

    exclude2 = false(size(only_in_2));
    for j = 1:numel(matchedFolders2)
        mf = matchedFolders2{j};
        exclude2 = exclude2 | strcmp(only_in_2, mf) | strncmp(only_in_2, [mf '\'], length(mf)+1);
    end

    truly_missing_1 = only_in_1(~exclude1);
    truly_missing_2 = only_in_2(~exclude2);
else
    truly_missing_1 = only_in_1;
    truly_missing_2 = only_in_2;
end

% --- Save to Excel: one sheet per nesting level (renamed items excluded) ---
getLevel = @(p) numel(strsplit(p, {'\\', '/'}));

if ~isempty(only_in_1) || ~isempty(only_in_2)
    if exist(outputFile, 'file')
        delete(outputFile);
    end

    allPaths = [truly_missing_1(:); truly_missing_2(:)];
    if ~isempty(allPaths)
        maxLevel = max(cellfun(getLevel, allPaths));
        wroteAny = false;
        for lv = 1:maxLevel
            in1 = truly_missing_1(cellfun(getLevel, truly_missing_1) == lv);
            in2 = truly_missing_2(cellfun(getLevel, truly_missing_2) == lv);
            if isempty(in1) && isempty(in2)
                continue;
            end
            nRows = max(numel(in1), numel(in2));
            in1(end+1:nRows) = {''};
            in2(end+1:nRows) = {''};
            header = {['Only in Dir1 (level ' num2str(lv) ')'], ['Only in Dir2 (level ' num2str(lv) ')']};
            writecell([header; in1(:), in2(:)], outputFile, 'Sheet', ['Level_' num2str(lv)]);
            wroteAny = true;
        end
        if wroteAny
            fprintf('Truly missing items saved by level.\n');
        end
    end

    if ~isempty(potentialRenames)
        header = {'Dir1 path (long name)', 'Dir2 path (short name)', 'Items in Dir1 subtree', 'Items in Dir2 subtree', 'Same count?'};
        writecell([header; potentialRenames], outputFile, 'Sheet', 'PotentialRenames');
        fprintf('Found %d potential renames. See PotentialRenames sheet.\n', nMatches);
    end

    fprintf('Results saved to: %s\n', outputFile);
end

function entries = getEntries(rootDir)
    entries = {};
    pending = {''};  % relative paths of dirs still to scan
    while ~isempty(pending)
        relDir = pending{1};
        pending(1) = [];
        absDir = fullfile(rootDir, relDir);
        items = dir(absDir);
        items = items(~ismember({items.name}, {'.', '..'}));
        for i = 1:numel(items)
            if isempty(relDir)
                relPath = items(i).name;
            else
                relPath = fullfile(relDir, items(i).name);
            end
            entries{end+1} = relPath; %#ok<AGROW>
            if items(i).isdir
                pending{end+1} = relPath; %#ok<AGROW>
            end
        end
    end
    entries = sort(entries(:));
end
