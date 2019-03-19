format long;

% Load the model if not already done so
if ~exist('model', 'var')
    model=mphload('../../../Tunable Dissipation/comsol/3DCavity.mph');
end

% desired chip-to-chip gaps and alignment offsets to simulate
washers = 0:2:10; % in mm
% [short, long]
pin_l = [2.03 4.56];
barrel_l = [4.57 5.64];
socket_depth = 9;
t_washer = 0.71;

% import custom progressbar package and initialize other variables
progressbar('Overall Progress','Runs Failed');
J=length(washers);
K=length(pin_l);
runs_complete = 0;
failures = 0;

% note: if a simulation fails, the corresponding entry in this matrix will
%       be left at zero (should be obviously different from successes)
Qs=zeros(J,K);
tic;

% loop over the different port lengths
for k=1:K
    % loop over lengths from pin to cavity
    for j=1:J
        model.param.set('port_l', [num2str(socket_depth - barrel_l(k) + j*t_washer) '[mm]']);
        model.param.set('pin_in_l', [num2str(pin_l(k)) '[mm]']);
        model.param.set('pin_out_l', [num2str(0.5) '[mm]']); % num2str(port_l(k)-pin_to_cavity_d(j)-1)
        model.mesh('mesh1').run()

        % need to use try-catch here in case the geometry fails to mesh
        try
            model.study('std1').run();
            q = mphglobal(model,'real(emw.Qfactor)','dataset','dset3');
            Qs(j,k) = q(1);
        catch

            % in the case of an error, the mutual for that run will be 0
            warning('Run failed');
            failures = failures + 1;
            progressbar([], failures/(J*K))
        end

        figure(8)
        set(gca,'FontSize',24)
        plot(washers,Qs,'LineWidth',2)
        xlabel('Number of washers');
        ylabel('Q');

        % increment the progressbar
        runs_complete = runs_complete + 1;
        progressbar(runs_complete/(J*K), []) 
    end
end

% end the timer and finalize the progressbar
toc
progressbar(1, []);

