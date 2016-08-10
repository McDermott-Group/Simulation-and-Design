function [indLength] = calcMutual(w1Width,s1Width,dWidth,s2Width,w2Width)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code contained herein calculates the mutual inductance between 2    %
% superconducting CPW traces separated by a a strip of ground plane of    %
% width d. The inputs of the function as for the center trace width of    %
% CPW 1 (w1Width), the gap width of CPW 1 (s1Width), the separation       %
% between center traces (dWidth), the gap width of CPW 2 (s2Width) and    %
% finally the width of the CPW 2 center trace. All these values are       %
% expected to be in microns.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%_____________      _________     _______      ________      __________   %
%      ground |<-->|<--W1-->|<-->|<--d-->|<-->|<--W2-->|<-->| ground      %
%             | S1 |        | S1 |       | S2 |        | S2 |             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4*pi*1e-7;            % perm of Free Space
w1Width = w1Width*1000;     % convert center trace 1 into nm
w2Width = w2Width*1000;     % convert center trace 2 int nm
s1Width = s1Width*1000;     % Trace 1 Gap width into nm
s2Width = s2Width*1000;     % Trace 2 Gap width into nm
dWidth = dWidth*1000;       % convert CPW's spacing into nm
mult = 10;                  % Ground Plane multiplier
gWidth = mult*w1Width;      % Ground plane is mult larger than center trace
points = 300;               % num points in center traces
r0 = 1e6;                   % return current loop (arbitrarily large)
lambda = 88;                % penetration depth (nm, 100 nm thick Nb)
w1Points = points;                              % Number of points in trace 1
w2Points = round(w2Width/w1Width)*points;       % Number of points in trace 2
groundNumPoints = mult*points;                  % Number of points in groundplane
dNumPoints = round(dWidth/w1Width)*points;      % Number of points in metal separation


% Set up the 1st ground plane
ground1x = linspace(-dWidth/2 - s1Width - w1Width - s1Width - gWidth,...
                    -dWidth/2 - s1Width - w1Width - s1Width,groundNumPoints);

% Set up Trace 1
trace1x = linspace(-dWidth/2 - s1Width - w1Width,-dWidth/2 - s1Width,w1Points);

% Set up separation
separationX = linspace(-dWidth/2,dWidth/2,dNumPoints);

% Set up trace 2
trace2x = linspace(dWidth/2 + s2Width,dWidth/2 + s2Width + w2Width,w2Points);

% Set up 2nd ground plane 
ground2x = linspace(dWidth/2 + s2Width + w2Width + s2Width,...
                    dWidth/2 + s2Width + w2Width + s2Width + gWidth,groundNumPoints);

% Combine all vectors into master vector
xVec = [ground1x trace1x separationX trace2x ground2x];

% 0 Thickness approximation
yVec = 0;

% Dummy vector for KRON multiplication
xDum = ones(1,length(xVec));

% Dummy for KRON multiplication
yDum = 1;

% Distance Matrix
rMat = kron(xVec,yDum);

% Distance Vector
rVec = reshape(rMat',length(xVec)*length(yVec),1);

% Memory management
clear xVec yVec rMat xDum yDum

% Dummy mutual vector for KRON
mdum=ones(1,length(rVec)); 

% Distances between vector sites
r_ij=abs(kron(mdum,rVec)-kron(conj(rVec'),mdum'));

% Distance step in ground plane
dxG = gWidth/groundNumPoints;

% Distance step in trace separation
dxD = dWidth/dNumPoints;

% Distance step in Trace 1
dxW1 = w1Width/w1Points;

% Distance step in Trace 2
dxW2 = w2Width/w2Points;

% Self inductance of ground plane 1
selfGround1 = dxG/2*ones(1,length(ground1x));

% Seld inductance of Trace 1
selfW1 = dxW1/2*ones(1,length(trace1x));

% Self inductance of separation
selfD = dxD/2*ones(1,length(separationX));

% Self inductance of Trace 1
selfW2 = dxW2/2*ones(1,length(trace2x));

% Self inductance of ground plane 2
selfGround2 = dxG/2*ones(1,length(ground2x));

% Master self inductance vector
self = [selfGround1 selfW1 selfD selfW2 selfGround2];

% Add self inductance of diagonal of site matrix
r_ij = r_ij + diag(self);

% Induce return current radius
m_exp = r0./r_ij;

% Geometric inductance contribution
mGeo = mu0/2/pi*log(m_exp);

% Kinetic inductance contribution
lKinVec = mu0*lambda^2/dxW1^2*ones(1,length(rVec));

% put kinetic inductancde on diagonal
lKin = diag(lKinVec);

% Total mutual inductance matrix
m = mGeo + lKin;

% memory management
clear mGeo lKinVec lKin r_ij self mdum m_exp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code is how we actually compute the mutual between        %
% traces 1 and 2. What happens is that we are looking to solve the Eqn    %
% Phi  = I * M. Specifically, we want the element M_(1,2) which satisfies %
% Phi_2 = I_1 * M_(1,2). In order to do this, we set all the ground plane %
% phi values to -1 and set trace 1's phi value to +1. We then adjust the  %
% phi value for trace 2 until the current I_2 is zero which enforces that %
% all the couling is between 1 --> 2 with no 2 --> 1 contribution.        %
% We them take the value for phi2 and divide it by the sum of currents in %
% trace one and get the mutual element we desire.                         %
% 
% A numerical root finding technique is employed where the value of phi2  %
% is incremented between -1 to 1 and the current in the trace 2 is calc'd %
% The algorithm checks for parity between sucessive elements in the       %
% calc'd currents vector of trace 2 to find when it has crossed  0. The   %
% phi value is then reset to its previous value and the incremental change% 
% in phi is reduced by a factor of 10. The routine breaks when the value  %
% of I in trace 2 is below 1E-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phi value for ground plane 1
phiG1 = -1*ones(length(ground1x),1);

% phi value for trace 1
phiW1 = 1*ones(length(trace1x),1);

% phi value for separation
phiD = -1*ones(length(separationX),1);

% phi value for ground plane 2
phiG2 = -1*ones(length(ground2x),1);

% i is the incremented value for phi of trace 2
i = -1;

% j is a counter that simply keeps
j = 1;

% incremental change in value for i (phi value for trace 2).
di = .05;

while i < 1 % Max value of i should not exceed 1
   
   % Set current value of phi 2 vector 
   phiW2 = i*ones(length(trace2x),1);
   
   % update master phi vector
   phi = [phiG1' phiW1' phiD' phiW2' phiG2']';
   
   % right matrix divide master phi vector by mutual inductance matrix
   current = m\phi;
  
   % Sum up current in trace 2
   iSum(j) = sum(current(groundNumPoints+points+dNumPoints+1:groundNumPoints+points+dNumPoints+points));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Portion of while loop that checks parity of succesive trace 2 currents  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if j>1
        % Most current value for current in trace 2
        curValue = iSum(j);
        
        % previous value for current 
        pastValue = iSum(j-1);
        
        % multiply last two values of trace 2 current vector
        sign = curValue*pastValue;
        
        % if I_i * I_(i-1) < 0, current has gone through 0 
        if sign < 0
            
            % Step back in phi value
            i = i-di;
            
            % reduce step in phi value by factor of 10
            di = di/10;
            
            % Step back trace 2 current counter
            j = j-1;
         
        end
        
        % Check to see if current in trace 2 is zero-ish
        if  abs(curValue) < 1e-6
           iCurrent(j) = i;
           iCurrent(j+1) = i;
           % Break the loop is current in trace 2 is zero
           break
        end
      
  end
  
  % update phi 2 value vector
  iCurrent(j) = i;
  
  % update j
  j = j+1;
  
  % update phi 2 value
  i = i + di; 
  
  clear phiW2;
  
end

% Set final value of trace 2 phi
phiW2I = iCurrent(end);

% Set final trace 2 phi vector
phiW2 = phiW2I*ones(length(trace2x),1);

% set final phi master vector
phi = [phiG1' phiW1' phiD' phiW2' phiG2']';

% calculate current
current = m\phi;

% Sum up currents in trace 1
i1Sum = sum(current(groundNumPoints+1:groundNumPoints+points));

% M_(1,2) = phi_2 / I_1
indLength = (phiW2I+1)/i1Sum;


end