format long;

% Load the model if not already done so
if ~exist('model', 'var')
    model=mphload('/Users/chriswilen/Documents/grad-school/research/mcdermott/Tunable Dissipation/comsol/cfm.mph');
end

x_list = -246:8:0;
y_list = -99:8:-11;

% import custom progressbar package and initialize other variables
progressbar('Overall Progress','Runs Failed');
runs_complete = 0;
failures = 0;

q1 = zeros(length(x_list),length(y_list));
q2 = zeros(length(x_list),length(y_list));
q3 = zeros(length(x_list),length(y_list));

tic;

for i = 1:length(x_list)
    for j = 1:length(y_list)
        x = x_list(i);
        y = y_list(j);
        
        % do not evaluate on the cap pad
        if ((x + 0.4 > -90) && (y + 0.4 > -63) && (y - 0.4 < -23) )
            continue
        end
        
        model.param.set('xx', num2str(x));
        model.param.set('yy', num2str(y));
        model.mesh('mesh1').run();

        % need to use try-catch here in case the geometry fails to mesh
        try
            model.study('std1').run();
            model.result().numerical("gev1").run();
            q1(i,j) = mphglobal(model,'es.Q0_1','dataset','dset1');
            q2(i,j) = mphglobal(model,'es.Q0_2','dataset','dset1');
            q3(i,j) = mphglobal(model,'es.Q0_3','dataset','dset1');
        catch

            % in the case of an error, the mutual for that run will be 0
            warning('Run failed');
            failures = failures + 1;
            progressbar([], failures/(length(x_list)*length(y_list)))
        end

%         figure(8)
%         set(gca,'FontSize',24)
%         plot(washers,Qs,'LineWidth',2)
%         xlabel('Number of washers');
%         ylabel('Q');

        % increment the progressbar
        runs_complete = runs_complete + 1;
        progressbar(runs_complete/(length(x_list)*length(y_list)), []) 
    end
end

% end the timer and finalize the progressbar
toc
progressbar(1, []);

