function params = set_paras(n,coup, kmax)
% set parameters for the simulation of a network of hopf normal form
% coup = coupling and can be
%   all: all to all connected
%   chain: unidirectional chain
%   dis: fully disconnected

% parameters in the system
%params.alpha = 0.03;                   % noise amplitude based on prev calcs of mean esc in bist 2 node
%params.nu = 0.1;                 

% parameters for the simulations
params.h = 1e-3;                       % time step 0.001
params.kmax = kmax;                    % max number of realisations for each beta
params.thresh1 = 0.55;              %
params.thresh2 = 0;                   % escape threshold

% for the coupling
%params.beta =[0:0.05:0.5 0.6];%[0 0.05 0.1  0.15:0.01:0.25 0.3 0.4 0.5]; % 0.01; 0.1; 1.0];    % coupling strength logspace(-3,2, 11)

if isequal(coup,'all')
    A = ones(n,n) - eye(n,n);
elseif isequal(coup,'chain') % chain coupling 1<-2<-3<-4<-5<-...<-n
    A = zeros(n,n);
    for k=1:n-1 % 
        A(mod(k,n)+1,k) = 1;
    end
elseif isequal(coup,'ring') % ring coupling 1->2->3->4->5->...->n->1
    A = zeros(n,n);
    for k=1:n
        A(k,mod(k,n)+1) = 1;
    end
elseif isequal(coup,'dis')
    A = zeros(n,n);
end

% expand to k realisations all at once
I = eye(kmax);
Aext = kron(I,A); % repeat on the diag
Aext = sparse(Aext);

params.A = A;
params.Aext = Aext;


