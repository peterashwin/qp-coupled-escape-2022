%% stochastic simulation of BIST model UNI - Heun method
% comuting realisations and monitoring the order/time of escape
% n-nodes
%
% Jennifer Creaser 
% July improved code and now do kmax realisations at once
% the coupling function should be correct now.
% Jan 2022 modified to run diff beta, alpha and nu values through bash
% Aug 2022 saves whole sequence of escape times

% clear
% close all
% %parsed in with the bash:
% nu = 0.01
% alpha = 0.05
% beta = 0.18

n = 3; % number of nodes in network

kmax = 5000;    % how many to compute at once
FIG = 0;            % set to 0 to not plot if kmax>1

k = kmax*n;     % how many rows to start

% set parameters for type of network (all, chain, ring, dis)
coup = 'chain'; % chain is uni
paras = set_paras(n, coup, kmax);

%%%% parameter values parsed in by the bash %%%%%
paras.nu = nu;
paras.beta = beta;
paras.alpha = alpha;

% node function - bistable model
nodeFunc = @(x,paras) -(x - 1).*(x.^2 - paras.nu);

% coupling function - diffusive
coupFunc = @(x,beta,B)  beta*(paras.Aext)*x - B.*x;

% random numbers
rng(0);
rmax = 10000;

init_eq = -sqrt(paras.nu);


tic
%% compute realisations

% find mean escape time for each beta value
beta = paras.beta;

B = beta.* sum( paras.Aext, 2 );     % coupling strength and structure fixed for each beta

% set up noise vals - this takes time
myrands1 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
myrands2 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
r=1;

y = ones(k,1).*(init_eq);                 % initial value
t = 0;                          % initial time

tau = cell(k,1);  tau(:,:) = {zeros(1)};  % time array
X = cell(k,1); X(:,:) = {zeros(1)}; % state array

if FIG==1; figure; hold on; end

while 1                 % no maximum time

    if r>rmax
        myrands1 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
        myrands2 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
        r=1; % reset r
        %fprintf('more rands! \n')
    end

    % evaluate slope of deterministic bit left side of interval
    fLeft = nodeFunc(y,paras) + coupFunc(y,beta,B);
    % prediction euler's step
    yBar = y + paras.h.*fLeft + myrands1(:,r);
    % correction
    fRight = nodeFunc(yBar,paras) + coupFunc(yBar,beta,B);
    % save current values for interpolation
    yOrig = y;
    % next values
    y = y + paras.h .* (fLeft+fRight) ./ 2 +  myrands2(:,r);

    % escape test
    clear escNode escZeros escTest
    escNode = find(y > paras.thresh1)';     % get index of escaped nodes
    escZeros = cellfun(@(x) x(end)==0,X);   % find all state 0 (unescaped)
    escTest = find(escZeros(escNode) == 1 );% check just escaped

    if ~isempty(escTest)
        for Nextesc = 1:length(escTest)
            X{escNode(escTest(Nextesc))}(end+1) = 1; % add 1
            % linear interpolation for escape times
            timeInterp = (paras.thresh1 - yOrig(escNode(escTest(Nextesc))))./(y(escNode(escTest(Nextesc)))-yOrig(escNode(escTest(Nextesc))));
            % add time
            tau{escNode(escTest(Nextesc))}(end+1)= t + timeInterp.*paras.h;
        end
    end


    % backwards escape test
    clear escNode escZeros escTest
    escNode = find(y < paras.thresh2)';                 % if one or more nodes are NOT past the threshold
    escZeros = cellfun(@(x) x(end)==1,X);   % find all state 1 (escaped)
    escTest = find(escZeros(escNode) == 1 ); % check just crossed

    if ~isempty(escTest)
        for Nextesc = 1:length(escTest)
            X{escNode(escTest(Nextesc))}(end+1) = 0; % add 0 to mark unescaped
            % linear interpolation for escape times
            timeInterp = (paras.thresh1 - yOrig(escNode(escTest(Nextesc))))./(y(escNode(escTest(Nextesc)))-yOrig(escNode(escTest(Nextesc))));
            % add time
            tau{escNode(escTest(Nextesc))}(end+1)= t + timeInterp.*paras.h;
        end
    end

    t=t+paras.h;                          % next time step
    r = r+1;

    % check if all escaped
    escZeros = cellfun(@(x) x(end)==1,X);   % find all state 1 (escaped)
    if ~any(escZeros==0)% all nodes have escaped
        break;
    end
end

Taunam=sprintf(['bist_returns_' num2str(n) coup '_kmax' num2str(kmax) '_nu' strrep(num2str(paras.nu),'.','pt') '_alpha' strrep(num2str(alpha),'.','pt') '_beta' strrep(num2str(beta),'.','pt')]);
save(Taunam,'tau')


%fprintf(['\n Beta ' num2str(beta) ' done \n'])



