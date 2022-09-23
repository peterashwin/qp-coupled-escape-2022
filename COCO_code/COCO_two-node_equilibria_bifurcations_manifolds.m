%% using coco to compute 1D manifolds as BVPs - bistable_2node Uni
clear
close all
clc

%% Values to compute

BetaVals = [0.005 0.01:0.01:0.4];% [0.001 0.04 0.1 0.25 0.4];%  % the beta value you want
PLOTFIG = 0; % if 1 then plots and saves the phase plane plot of eqs (and mans if computed)
PLOTBIF = 1;
SAVEEQS = 0;
SAVESTAB = 0; % UnixQ_stableEQS.dat <- all stable EQ points
COMPMANS = 0; % if 1 computes the manifolds, otherwise just eqs

%% set up

EQ = [-0.1 0.1 1];
delta = 1e-4; % starting dist from eq

StbMk = {'square','^','o'};

% 2 node uni bistable model
func = @(x,p) [-(x(1,:)-1.0).*(x(1,:).^2-p(1,:)) + p(2,:).*(x(2,:)-x(1,:)); % vectorised and autonomous
               -(x(2,:)-1.0).*(x(2,:).^2-p(1,:))                          ]; % sqrt(p(1,:)) is the Q input from

bistJac = @(x,p) [ -3.*(x(1,:).^2) + p(1,:) + 2.*x(1,:) - p(2,:),   p(2,:);
                   0,                                               -3.*x(2,:).^2 + p(1,:) + 2.*x(2,:)];

%FigH = figure('DefaultAxesPosition', [0.02, 0.1, 0.95, 0.9]); hold on
%t = tiledlayout(2,13,'TileSpacing','none');
%% Loop over B's

for BB = 1:length(BetaVals)
    beta_stop = BetaVals(BB);
    EQS = []; % blank EQ values
    clear  mans stEQs   % blank mans values
    mm=1;
    if PLOTFIG; figure(2); clf; end

    for i=1:3
        for j=1:3 % for each equilibria

            pnames = {'nu','beta'};
            p0 = [0.01, 0];      % nu = 0.01, beta = 0 starting params
            x0 = [EQ(i), EQ(j)];

            %continue the equilibirium to suitable beta value
            prob = coco_prob();
            prob = coco_set(prob, 'ode', 'vectorized', false);
            ode_args = {func, x0, pnames, p0};
            cont_args = {1, 'beta', [0 beta_stop]}; % beta stopping condition in here
            bd = coco(prob, 'test2d', @ode_isol2ep, ode_args{:}, cont_args{:});

            if PLOTBIF
                figure(1); hold on % bifurcation diagram part
                thm = struct('special', {{'SN','BP'}});
                coco_plot_bd(thm, 'test2d', 'beta', 'x')
            end

            % extract eq values and compute stability:
            bs = cell2mat(bd(2:end,9));                     % list of beta values
            idx = find(abs(bs - beta_stop) <1e-8);          % could be multiple?
            if ~isempty(idx)
                for k = 1%:length(idx)
                    eq =  bd{idx(k)+1,19}; % eq value

                    p1 = [0.01; beta_stop];

                    J = bistJac(eq,p1);
                    [vs,evs,stab] = checkStab(J);

                    EQS(end+1,:) = [eq' stab]; % save eq and stab valu

                    if PLOTFIG
                        figure(2); hold on % plot of phase space
                        plot(eq(1),eq(2),StbMk{stab},'color','k','markerfacecolor','k','linewidth',2,'markersize',10); % plot the equilibria
                    end

                    if COMPMANS
                        for vcol = 1:2
                            if evs(vcol,vcol)<0 % stable integrate chenge which end you follow

                                for sgn = [-1 1]

                                    x0 = (eq + delta.*sgn.*vs(:,vcol))';  % must be row vector
                                    t0 = 0;                          % one t value for each row

                                    prob = coco_prob();
                                    prob = ode_isol2coll(prob, '', func, t0, x0, pnames, p1);
                                    data = coco_get_func_data(prob, 'coll', 'data');
                                    maps = data.coll_seg.maps;  % this gets all the relevent bits of the solution

                                    prob = coco_add_pars(prob, 'pars', [maps.x0_idx; maps.x1_idx; maps.T_idx], ...
                                        {'y1s' 'y2s' 'y1e' 'y2e' 'T'});
                                    prob = coco_set(prob, 'cont', 'ItMX', 1000, 'NPR', 100, 'NAdapt', 2); % need even smaller NAdapt here
                                    cont_args = {1, {'y1s' 'T' 'y2s' 'coll.err' 'coll.err_TF'}, {[-0.3 1.3] [0 300]}}; % continuation in T
                                    bd2 = coco(prob, 'coll1', [], cont_args{:});

                                    man = [cell2mat(bd2(2:end,13)), cell2mat(bd2(2:end,14)) ];
                                    mans{mm,2} = 1;
                                    mans{mm,1} = man;
                                    mm=mm+1;

                                    if PLOTFIG
                                        plot(cell2mat(bd2(2:end,13)),cell2mat(bd2(2:end,14)),'color',[0 0 0.5], 'linewidth',2)
                                    end
                                end


                            else
                                for sgn = [-1 1] % each sig of eq

                                    x0 = (eq + delta.*sgn.*vs(:,vcol))';  % must be row vector
                                    t0 = 0;                          % one t value for each row

                                    prob = coco_prob();
                                    prob = ode_isol2coll(prob, '', func, t0, x0, pnames, p1);
                                    data = coco_get_func_data(prob, 'coll', 'data');
                                    maps = data.coll_seg.maps;  % this gets all the relevent bits of the solution

                                    prob = coco_add_pars(prob, 'pars', [maps.x0_idx; maps.x1_idx; maps.T_idx], ...
                                        {'y1s' 'y2s' 'y1e' 'y2e' 'T'});
                                    prob = coco_set(prob, 'cont', 'ItMX', 1000, 'NPR', 100, 'NAdapt', 5);
                                    cont_args = {1, {'y1e' 'T' 'y2e' 'coll.err' 'coll.err_TF'}, {[-0.3 1.3] [0 300]}}; % continuation in T {[0 1] [0 100]}
                                    bd2 = coco(prob, 'coll1', [], cont_args{:});

                                    man = [cell2mat(bd2(2:end,15)), cell2mat(bd2(2:end,16)) ];
                                    mans{mm,2} = 1;
                                    mans{mm,1} = man;
                                    mm=mm+1;

                                    if PLOTFIG
                                        plot(cell2mat(bd2(2:end,15)),cell2mat(bd2(2:end,16)),'color',[0.5 0 0], 'linewidth',2)
                                    end

                                end
                            end
                        end
                    end

                end
            end
        end
    end

    fprintf(['\n' strrep(num2str(beta_stop),'.','pt') ' done \n'])

    if PLOTFIG
        figure(2); box on
        axis([-0.2 1.2 -0.2 1.2 ])
        xlabel('x1'); ylabel('x2','rotation',0)
        set(gca,'fontsize',16,'linewidth',1)
        title(['beta ' num2str(beta_stop)])
        hgsave(gcf,['Uni2D_beta' strrep(num2str(beta_stop),'.','pt') '.fig']);
    end

    if COMPMANS; save(['MansUni2D_beta' strrep(num2str(beta_stop),'.','pt') '.mat'],'mans'); end

    if SAVEEQS; save(['EQSUni2D_beta' strrep(num2str(beta_stop),'.','pt') '.mat'],'EQS'); end

    if SAVESTAB
        % extract stable points for QP computation
        %xvals = EQS(EQS(:,3)==3,1:2);
        %xvals = xvals(xvals(:,1) - xvals(:,2)>1e-5,:); % only the AQ state if it exists
        
        % for all stable
        xvals = reshape(EQS(EQS(:,3)==3,1:2)',1,[]);
        
        stEQs= [beta_stop xvals ];

        fileID = fopen('Uni2DstableEQS_all.dat','a'); %'w' is for write(over) 'a' is for append
        fprintf(fileID,[repmat('%12.8f ',[1,length(stEQs)]) '\n'],stEQs);
        fclose(fileID);
    end
    
    
    
    
end







%%
function [vs,evs,stab] = checkStab(J)    % can compare to the 'eigs' col in bd
[vs,evs] = eig(J);    % can compare to the 'eigs' col in bd
id1 = sign(evs(1,1)); id2 = sign(evs(2,2));
if id1>0 && id2>0; stab = 1;
elseif (id1>0 && id2<0) || (id1<0 && id2>0); stab = 2;
else; stab = 3;
end

end






