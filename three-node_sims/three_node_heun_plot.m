%% to plot the mean escape times of the BIST model with two nodes
% code to plot the second run of realisations
% Jen Creaser July 2021
%
close all
clc
clear

%% load the data

files=dir('bist_returns_3chain_kmax10000*'); % the file names

filenames={files.name}; % lists the filenames in an array
nof=numel(filenames); % number of files

% extract values
setUp = split(filenames{1},'_');
n = str2double(setUp{3}(1)); % number of nodes in network
coup = setUp{3}(2:end); % chain is uni
kmax = str2double(setUp{4}(5:end)); %2000;    % how many to compute at once
paras = set_paras(n, coup, kmax);
paras.nu = str2double(strrep(setUp{5}(3:end),'pt','.'));
paras.alpha = str2double(strrep(setUp{6}(6:end),'pt','.'));

% extract all beta values
for k=1:nof %for each file.
    fileParts = split(filenames{k}(1:end-4),'_');
    betaValues(k) = str2double(strrep(fileParts{7}(5:end),'pt','.'));
end
allVals = sort(betaValues);

% set up seqs with perms
allSeqs = num2cell(perms([1 2 3])',[1,3])';
allSeqs = cellfun(@(x) x',allSeqs,'UniformOutput',false);

%%
for k=1:nof %for each file.
    clear Times Order tau
    fileParts = split(filenames{k}(1:end-4),'_');
    beta= str2double(strrep(fileParts{7}(5:end),'pt','.'));
    b = find(abs(allVals - beta) <1e-5);

    load(filenames{k}); % tau

    tau = reshape(tau,[n,kmax]);
    for col  = 1:kmax
        realisTimes = []; realisOrder = [];
        for node = 1:n
            times = tau{node,col};
            order = ones(size(times)).*(4-node);
            order(1:2:end) = order(1:2:end)*-1; % make odds -ve marks returns
            realisTimes = [realisTimes times];
            realisOrder = [realisOrder order];
        end
        [tvals,j] = sort(realisTimes);
        Times{col,1} = realisTimes(j(4:end));% ignoring first 3 zeros
        Order{col,1} = realisOrder(j(4:end));% ignoring first 3 zeros
    end

    % percentage with returns
    IndReturns = cellfun(@(x) length(x)>3,Times);
    PercReturns(b) = (sum(IndReturns)/length(IndReturns))*100;

    % escape sequences
    IndLengths = unique(cellfun(@(x) length(x),Order));
    % for each length identify unique seqs and count how many
    for l = 1:length(IndLengths)
        lengthL = cellfun(@(x) length(x)==IndLengths(l),Order);
        SeqsL = table2cell(unique(cell2table(Order(lengthL))));
        for stest = 1:length(SeqsL)
            % compare each seq to allSeqs list and save new ones
            seqtest = sum(cellfun(@(x) isequal(x,SeqsL{stest}),allSeqs));
            if seqtest==0
                allSeqs{end+1} = SeqsL{stest};
            end
        end
    end

    for s = 1:length(allSeqs)
        SeqCounts(s,b) = sum(cellfun(@(x) isequal(x,allSeqs{s}),Order));
        SeqProbs(s,b) = sum(cellfun(@(x) isequal(x,allSeqs{s}),Order))/kmax;
    end

    % direction of first escape
    for n=1:3; FirstDir(n,b) = sum(cellfun(@(x) x(1)==n,Order))/kmax; end

end



%%

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

cols = bone(9);
bifcols = lines(9);
seqcols = bone(length(allSeqs)+2);
close all

%% Plot perc returns
FigH = figure('DefaultAxesPosition',[0.15 0.23 0.8 0.72]);
b = bar(allVals,PercReturns);
b.FaceColor = cols(5,:);
b.LineWidth = 1;
hold on
plot([0.0101 0.0101],[0 3.5],':','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[0 3.5],':','color',bifcols(5,:),'linewidth',3)
plot([0.1528 0.1528],[0 3.5],'--','color',[0.5 0.5 0.5],'linewidth',3) % Gate height
xlabel('$\beta$','interpreter','latex')
ylabel('$\%$ of returns','interpreter','latex')
%xlim([0 0.4]); 
set(gca,'ygrid','on'); %ylim([0 1])
set(gca,'linewidth',2); box on;
nam = 'Bist_returns_3chain_perc';
hgsave(gcf,[nam '.fig']);
s=hgexport('readstyle','18x10x20'); %read the style 18cm x 10cm with 20 font
hgexport(gcf,[nam '.eps'],s);

%% plot dir of first escape

FigH = figure('DefaultAxesPosition', [0.15 0.23 0.8 0.72]);
b= bar(allVals,FirstDir,'lines','none');
for k = 1:3
    b(k).FaceColor = cols(k+3,:);
    b(k).LineWidth = 1;
end
ylabel('P(dir. first escape)','interpreter','latex')
hold on
plot([0.0101 0.0101],[0 1.2],':','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[0 1.2],':','color',bifcols(5,:),'linewidth',3)
plot([0.1528 0.1528],[0 1.2],'--','color',[0.5 0.5 0.5],'linewidth',3) % Gate height
set(gca,'ygrid','on');
legend({'$x_1$','$x_2$','$x_3$'},'interpreter','latex','location','northwest')
xlabel('$\beta$','interpreter','latex')
xlim([-0.01 0.26]);ylim([0 1])
set(gca,'linewidth',2);box on
nam = 'Bist_returns_3chain_dir';
hgsave(gcf,[nam '.fig']);
s=hgexport('readstyle','18x10x20'); %read the style
hgexport(gcf,[nam '.eps'],s);

%% plot sequence probabilites

FigH = figure('DefaultAxesPosition', [0.15 0.23 0.8 0.72]);hold on
for n = 1:size(SeqCounts,1)
    plot(allVals,SeqCounts(n,:)./kmax,'o-','color',seqcols(n,:),'linewidth',2); 
end
plot([0.0101 0.0101],[.00005 30000],':','color',bifcols(5,:),'linewidth',3)
plot([0.2025 0.2025],[.00005 3e4],':','color',bifcols(5,:),'linewidth',3)
plot([0.1528 0.1528],[.00005 3e4],'--','color',[0.5 0.5 0.5],'linewidth',3) % Gate height

ylabel('Counts','interpreter','latex')
hold on
set(gca,'ygrid','on');
legend(cellfun(@(x) num2str(x),allSeqs,'UniformOutput',false),'interpreter','latex','location','northwest')
xlabel('$\beta$','interpreter','latex')
xlim([-0.01 0.26])
ylim([.00005 1])
set(gca,'YScale','log')
set(gca,'YMinorGrid','off')
set(gca,'linewidth',2); box on;
nam = 'Bist_returns_3chain_seqProb';
hgsave(gcf,[nam '.fig']);
s=hgexport('readstyle','18x10x20'); %read the style
hgexport(gcf,[nam '.eps'],s);



