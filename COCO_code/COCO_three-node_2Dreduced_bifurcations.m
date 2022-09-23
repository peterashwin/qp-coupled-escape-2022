%% bifurcation diagram in COCO
% prototypical bistable model
% 2 nodes

clear

%% Create empty coco problem and add Phi zero problem and initial guess for u0


EQ = [-0.1 0.1 1];

for i=1:3
    for j=1:3
        x0 = [EQ(i), EQ(j)];  % x1 = x2 = sqrt(nu)
        pnames = {'nu','beta'};
        p0 = [0.01, 0];      % nu = 0.01, beta = 0
        
        prob = coco_prob();
        prob = coco_set(prob, 'ode', 'vectorized', false);
        ode_fcns = {@bist};
        ode_args = {ode_fcns{:}, x0, pnames, p0};
        cont_args = {1, 'beta', [0 2]};
        
        coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:});
        
        %
        figure(4); hold on
        thm = struct('special', {{'SN','BP'}});
        coco_plot_bd(thm, 'test2d', 'beta', 'x')
        
        grid on
    end
end
%% Cont in Nu
EQ = [-0.1 0.1 1];

for i=1:3
    for j=1:3
        x0 = [EQ(i), EQ(j)];  % x1 = x2 = sqrt(nu)
        pnames = {'nu','beta'};
        p0 = [0.01, 0];      % nu = 0.01, beta = 0
        
        prob = coco_prob();
        prob = coco_set(prob, 'ode', 'vectorized', false);
        ode_fcns = {@bist};
        ode_args = {ode_fcns{:}, x0, pnames, p0};
        cont_args = {1, 'nu', [0 2]};
        
        coco(prob, 'nu2d', @ode_isol2ep, ode_args{:}, cont_args{:});
        
        %
        figure(5); hold on
        thm = struct('special', {{'SN','BP'}});
        coco_plot_bd(thm, 'nu2d', 'nu', 'x')
        
        grid on
    end
end

%%  for multiple nu
NUS = [0.001, 0.01, 0.1 0.2 0.5 0.9];
figure
t = tiledlayout(1,length(NUS));
t.TileSpacing = 'compact';
t.Padding = 'compact';

for nn = 1:length(NUS)
    nexttile; hold on;
    nu1 = NUS(nn)
EQ = [-sqrt(nu1) sqrt(nu1) 1];

for i=1:3
    for j=1:3
        x0 = [EQ(i), EQ(j)];  % x1 = x2 = sqrt(nu)
        pnames = {'nu','beta'};
        p0 = [nu1, 0];      % nu = 0.01, beta = 0
        
        prob = coco_prob();
        prob = coco_set(prob, 'ode', 'vectorized', false);
        ode_fcns = {@bist};
        ode_args = {ode_fcns{:}, x0, pnames, p0};
        cont_args = {1, 'beta', [0 0.5]};
        
        coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:});
        
        %

        thm = struct('special', {{'SN','BP'}});
        coco_plot_bd(thm, 'test2d', 'beta', 'x')
        
        grid on
    end
end
end



%% follow SN and BP points

% note to follow the BP might need to break symmetry slightly so it becomes two SNs

%% functions

function f = bist(x,p)

nu   = p(1,:);
beta = p(2,:);

x1 = x(1,:);
x2 = x(2,:);

f = [-(x1-1.0).*(x1.^2-nu) + beta.*(x2-x1);
    -(x2-1.0).*(x2.^2-nu) + beta.*(sqrt(nu)-x2)];

end
