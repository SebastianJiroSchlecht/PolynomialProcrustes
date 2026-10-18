function audit_figures
% Run each paper plotting script in a fresh workspace, retaining all failures.
addpath(fullfile(pwd,'PEVDToolbox','General'));
addpath(fullfile(pwd,'PEVDToolbox','Decompositions'));
set(groot,'defaultFigureVisible','off');
scripts = {'Generate_Figs12','Generate_Fig3','Generate_Figs4to7',...
    'Generate_Figs910','Generate_Fig12','Generate_Figs13to16',...
    'Generate_Fig17','Generate_Fig18','Generate_Figs1920',...
    'ProcrustesRandomMatrixEnsembleEval'};
diary('reproduction-log.txt');
cleanup = onCleanup(@() diary('off'));
failures = {};
for k = 1:numel(scripts)
    fprintf('\nRUN %s\n',scripts{k});
    try
        run_one(scripts{k});
        fprintf('PASS %s\n',scripts{k});
    catch exception
        fprintf('FAIL %s\n%s\n',scripts{k},getReport(exception));
        failures{end+1} = scripts{k}; %#ok<AGROW>
    end
end
assert(isempty(failures),'Figure reproduction failed: %s',strjoin(failures,', '));
end

function run_one(scriptName)
run(scriptName);
close all;
end
