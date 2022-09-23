%% stochastic simulation of BIST model UNI - Heun method
% comuting realisations and monitoring the order/time of escape
% n-nodes
%
% Jennifer Creaser May 2021
% July improved code and now do kmax realisations at once
% the coupling function should be correct now.
% Jan 2022 modified to run diff beta, alpha and nu values through bash

%clear 
%close all

n = 2; % number of nodes in network

kmax = 2000; %2000;    % how many to compute at once
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

rmax = 100000;
init_eq = -sqrt(paras.nu);



%% compute realisations

for j = 1:length(paras.beta)  %1:length(paras.beta)                    % find mean escape time for each beta value
    beta = paras.beta(j);
    
    tic,
    
    B = beta.* sum( paras.Aext, 2 );     % coupling strength and structure fixed for each beta
    
    %escTimes = zeros(kmax,n);        % reshape in post
    %escOrder = zeros(kmax,n);
    
    tic
    %for k=1:paras.kmax
    
    % set up noise vals - this takes time
    myrands1 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
    myrands2 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
    r=1;
    
    y = ones(k,1).*(init_eq);                 % initial value
    t = 0;                          % initial time
    %esc_count = 1;                  % counter for escapes
    
    tau = zeros(k,1);    % time placeholder 
    tau1 = zeros(k,1); % time of first escape
    X = zeros(k,1);              % initial state
    Z = zeros(k,1);              % initial state
    

    
    if FIG==1; figure; hold on; end
    
    while 1                 % no maximum time
    %while t < 10000
    
        if r>rmax
            myrands1 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
            myrands2 = paras.alpha .* sqrt(paras.h) .* randn(k,rmax);
            r=1; % reset r
            fprintf('more rands! \n')
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
        escNode = find(y > paras.thresh1)';                 % if one or more nodes are past the threshold
        escTest = find(X(escNode) == 0 );                   % get index on which one(s) have not escaped (sort of)
        if ~isempty(escTest)                                % check node escape
            X(escNode(escTest)) = 1;                        % flag as escaped
            Z(escNode(escTest)) = Z(escNode(escTest))+1;    % add up how many threshold crossings
            % linear interpolation for escape times
            timeInterp = (paras.thresh1 - (yOrig(escNode(escTest))))./((y(escNode(escTest)))-(yOrig(escNode(escTest))));
            tau(escNode(escTest))= t + timeInterp.*paras.h; % updates each time the threshold is crossed
            if tau1(escNode(escTest)) == 0; % record first crossing time only
                tau1(escNode(escTest)) = t + timeInterp.*paras.h;
            end
        end
        
        % backwards escape test
        escNode = find((y)<paras.thresh2)';                 % if one or more nodes are NOT past the threshold
        escTest = find(X(escNode) == 1 );                   % get index on which one(s) have escaped
        if ~isempty(escTest)                                % check node escape
            X(escNode(escTest)) = 0;                        % remove escape flag
            Z(escNode(escTest)) = Z(escNode(escTest))+1;    % add one to treshold counter
        end        
        
        t=t+paras.h;                          % next time step
        r = r+1;
        
        if ~any(X==0)% all nodes have escaped 
            break;
        end
    end
    
    Taunam=sprintf(['bist_Ztimes_' num2str(n) coup '_kmax' num2str(kmax) '_beta' strrep(num2str(beta),'.','pt') '_nu' strrep(num2str(paras.nu),'.','pt') '_alpha' strrep(num2str(alpha),'.','pt') '.dat']);
    fileID = fopen(Taunam,'w'); % creates file
    fprintf(fileID,'%15.15f %15.15f %15.15f\n',[Z tau1 tau]'); % each row becomes a column
    fclose(fileID);
    
    toc
    fprintf(['\n Beta ' num2str(beta) ' done \n'])
end


