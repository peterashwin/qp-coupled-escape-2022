%% to plot the mean escape times of the BIST model with two nodes
% code to plot the second run of realisations
% Jen Creaser July 2021
% 
close all
clc
clear

%% load the data
files=dir('bist_Ztimes_2chain_kmax2000_beta*'); % the file names 

filenames={files.name}; % lists the filenames in an array
nof=numel(filenames); % number of files

% extract values
n = str2double(filenames{1}(13)); % number of nodes in network
kmax = str2double(filenames{1}(24:27)); %2000;    % how many to compute at once
coup = filenames{1}(14:18); % chain is uni

paras = set_paras(n, coup, kmax);

% extract all beta values
for k=1:nof %for each file.
    file=filenames{k}    
    params = split(file(12:end-4),'_');
    bval(k) = str2double(strrep(params{4}(5:end),'pt','.'));
end
allVals = sort(bval);

%%

for k=1:nof %for each file.
    file=filenames{k}
    
    params = split(file(12:end-4),'_');

    beta = str2double(strrep(params{4}(5:end),'pt','.'));
    b = find(abs(allVals - beta) <1e-5);
    
    dat = load(file);
    Z = reshape(dat(:,1),[n,kmax]);
    tau1 = reshape(dat(:,2),[n,kmax]);
    tau = reshape(dat(:,3),[n,kmax]);
    
    % times for 'first' escape times
    [Times, Order]=  sort(tau1);   Times = Times';
    Escape1(b,1) = mean(Times(:,1));                % 1|0
    Escape1(b,2) = mean(Times(:,2) - Times(:,1));   % 2|1
    
    % Percentage of realisations that will go back after first escape
    P1(b,1) = (length(find(Order(1,find(Z(2,:)>1))==1))/kmax)*100; % first escape is node 1
    P1(b,2) = (length(find(Order(1,find(Z(2,:)>1))==2))/kmax)*100; % first escape is node 2
    
    % times for 'last' escape times
    [Times, Order]=  sort(tau);   Times = Times';
    Escape(b,1) = mean(Times(:,1));                % 1|0
    Escape(b,2) = mean(Times(:,2) - Times(:,1));   % 2|1
    
    Counter(b,1) = length(find(Z(1,:)>1))/kmax; % proportion of realisations that visted node 1 and returned to 0 before full escape
    Counter(b,2) = length(find(Z(2,:)>1))/kmax; % proportion of realisations that visted node 2 and returned to 0 before full escape
    Cmax(b,1) = mean(Z(1,find(Z(2,:)>1))); % maximum number of times realisation went back and forth between node 1 and 0 before escaping
    Cmax(b,2) = mean(Z(2,find(Z(2,:)>1))); % maximum number of times realisation went back and forth between node 2 and 0 before escaping
    
    % probability of final escape being node 1 or node 2
    P(b,1) = length(find(Order(1,:)==1))/kmax; % first escape is node 1
    P(b,2) = length(find(Order(1,:)==2))/kmax; % first escape is node 2

    % probability for only those realisations that returned
    %P(b,1) = length(find(Order(1,find(Z(2,:)>1))==1))/kmax; % first escape is node 1
    %P(b,2) = length(find(Order(1,find(Z(2,:)>1))==2))/kmax; % first escape is node 2
   
end



%%

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

cols = [0.5 0.5 0.5; 0.8 0.8 0.8];
bifcols = lines(8);
close all

FigH = figure('DefaultAxesPosition', [0.1 0.25 0.88 0.7]);
hold on 
plot(allVals, Escape1(:,1),'-','color',[0.75 0.75 0.75],'linewidth',3)
plot(allVals, Escape1(:,2),'--','color',[0.75 0.75 0.75],'linewidth',3)
plot(allVals, Escape(:,1),'-k','linewidth',3)
plot(allVals, Escape(:,2),'--k','linewidth',3)
plot([0.01 0.01],[0 120],'-','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[0 120],'-','color',bifcols(5,:),'linewidth',3)
plot([0.3025 0.3025],[0 120],'-','color',bifcols(5,:),'linewidth',3)
plot([0.18 0.18],[0 120],'-','color',bifcols(4,:),'linewidth',3)
legend({'first-first','second','first-last','second'},'location','east','interpreter','latex');
xlabel('$\beta$','interpreter','latex')
ylabel('Mean Escape Time','interpreter','latex')
box on
set(gca,'linewidth',2)
set(gca,'position',[0.15 0.23 0.8 0.72])
xlim([0 0.4])
hgsave(gcf,'Bist_Ztimes_2chain_kmax2000_nu0pt01_times.fig');
s=hgexport('readstyle','18x10x20'); %read the style
hgexport(gcf,['Bist_Ztimes_2chain_kmax2000_nu0pt01_times.eps'],s);

% FigH = figure('DefaultAxesPosition', [0.08, 0.08, 0.90, 0.85]);
% b= bar(allVals,Cmax); 
% for k = 1:2
%    b(k).FaceColor = cols(k,:);
% end
% hold on 
% plot([0.01 0.01],[0 4],'-','color',bifcols(4,:),'linewidth',1.5)
% plot([0.2025 0.2025],[0 4],'-','color',bifcols(4,:),'linewidth',1.5)
% plot([0.3025 0.3025],[0 4],'-','color',bifcols(4,:),'linewidth',1.5)
% plot([0.18 0.18],[0 4],'-','color',bifcols(2,:),'linewidth',1.5)
% ylabel('Max number of thresh crossings in each dir')
% xlabel('\beta')
% legend({'node 1','node 2'})
% xlim([0 0.4])
% hgsave(gcf,'Bist_Ztimes_2chain_kmax2000_nu0pt01_maxreturns.fig');

FigH = figure('DefaultAxesPosition',[0.1 0.25 0.88 0.7]); 
%subplot(r,c,2);
b= bar(allVals,P1(:,2));
b.FaceColor = [0.75 0.75 0.75];
b.LineWidth = 1;
ylabel('$\%$ of returns','interpreter','latex')
hold on 
plot([0.01 0.01],[0 1.2],'-','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[0 1.2],'-','color',bifcols(5,:),'linewidth',3)
plot([0.3025 0.3025],[0 1.2],'-','color',bifcols(5,:),'linewidth',3)
plot([0.18 0.18],[0 1.2],'-','color',bifcols(4,:),'linewidth',3)
%ylabel('%','interpreter','latex') % given that their first-first escape is node 2
xlabel('$\beta$','interpreter','latex')
xlim([0 0.4]); set(gca,'ygrid','on'); ylim([0 1])
set(gca,'linewidth',2)
set(gca,'position',[0.15 0.23 0.8 0.72])
hgsave(gcf,'Bist_Ztimes_2chain_kmax2000_nu0pt01_perc.fig');
s=hgexport('readstyle','18x10x20'); %read the style
hgexport(gcf,['Bist_Ztimes_2chain_kmax2000_nu0pt01_perc.eps'],s);


FigH = figure('DefaultAxesPosition', [0.1 0.25 0.88 0.7]); 
b= bar(allVals,P(:,[2,1]));
for k = 1:2
   b(k).FaceColor = cols(k,:);
   b(k).LineWidth = 1;
end
ylabel('P(dir. final escape)','interpreter','latex')
hold on 
plot([0.01 0.01],[0 1],'-','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[0 1],'-','color',bifcols(5,:),'linewidth',3)
plot([0.3025 0.3025],[0 1],'-','color',bifcols(5,:),'linewidth',3)
plot([0.18 0.18],[0 1],'-','color',bifcols(4,:),'linewidth',3)
%ylabel('P'); 
set(gca,'ygrid','on');
legend({'$x_1$','$x_2$'},'interpreter','latex','location','northwest')
xlabel('$\beta$','interpreter','latex')
xlim([-0.01 0.41])
set(gca,'linewidth',2)
set(gca,'position',[0.15 0.23 0.8 0.72])
hgsave(gcf,'Bist_Ztimes_2chain_kmax2000_nu0pt01_prob.fig');
s=hgexport('readstyle','18x10x20'); %read the style
hgexport(gcf,['Bist_Ztimes_2chain_kmax2000_nu0pt01_prob.eps'],s);





