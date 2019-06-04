format long;

% Load the model if not already done so
% if ~exist('model', 'var')
    model=mphload('/Users/chriswilen/Documents/grad-school/research/mcdermott/Tunable Dissipation/comsol/ringmon.mph');
% end

ri_list =     [40,   50, 60,  70,  80,  90,  100];
ro_min_list = [40.5, 52, 67,  85, 100, 150, 200];
ro_max_list = [43,   60, 85, 130, 350, 450, 650];
n = 10;

ro_list = ones(length(ri_list), n);
c = zeros(length(ri_list),n);

% import custom progressbar package and initialize other variables
progressbar('Overall Progress','Runs Failed');
runs_complete = 0;
failures = 0;


tic;

for i = 1:length(ri_list)
    ro_list(i,:) = linspace(ro_min_list(i), ro_max_list(i),n);
    for j = 1:n
        ri = ri_list(i);
        ro = ro_list(i,j);
        
        model.param.set('r_outer', num2str(ro));
        model.param.set('r_inner', num2str(ri));
        model.mesh('mesh1').run();

        % need to use try-catch here in case the geometry fails to mesh
        try
            model.study('std1').run();
            model.result().numerical("gev1").run();
            c(i,j) = mphglobal(model,'-es.C12','dataset','dset1');
        catch

            % in the case of an error, the mutual for that run will be 0
            warning('Run failed');
            failures = failures + 1;
            progressbar([], failures/(length(ri_list)*n))
        end

        figure(8)
        set(gca,'FontSize',24)
        plot(ro_list', 1e15*c','LineWidth',2)
        xlabel('Outer Diameter [um]');
        ylabel('C [fF]');

        % increment the progressbar
        runs_complete = runs_complete + 1;
        progressbar(runs_complete/(length(ri_list)*n), []) 
    end
end

% end the timer and finalize the progressbar
toc
progressbar(1, []);

