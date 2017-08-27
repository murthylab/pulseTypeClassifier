cc()
load('data/NM91_pulses.mat', 'pulseShapesNorm')
pulseShapesNorm = double(pulseShapesNorm);
load('data/NM91_pulseLabels.mat')
%% classification of P1/2
% project
templates = normalize(grpstats(pulseShapesNorm, pulseLabels));

projectedPulses = pulseShapesNorm*templates';
% xval QDA
for cvRun = 1:100
   c = cvpartition(length(pulseLabels), 'HoldOut',0.2);
   training = c.training;
   test = c.test;
   
   trainingShuffled = training( randperm(length(training)) );
   testShuffled     = test( randperm(length(test)) );
   qda = fitcdiscr(projectedPulses(training,:), pulseLabels(training), 'DiscrimType','quadratic');
   predLabels = qda.predict(projectedPulses(test,:));
   Ctmp = confusionmat(pulseLabels(test), predLabels);
   pcorr(cvRun) = sum(diag(Ctmp))./sum(Ctmp(:));
   C(:,:,cvRun) = Ctmp./sum(Ctmp(:));
   
   % test classifier with shuffled labels
   predLabels = qda.predict(projectedPulses(test ,:));
   Ctmp = confusionmat(pulseLabels(testShuffled), predLabels);
   pcorrShuffle(cvRun) = sum(diag(Ctmp))./sum(Ctmp(:));
   CShuffle(:,:,cvRun) = Ctmp./sum(Ctmp(:));
end
save('pulseTypeClassifier','qda', 'templates')
%% crossvalidate on other melanogaster strains
try % this will fail on most machines 
   strainLabels = {'NM91','CarM03', 'CM07', 'CSTul', 'N30', 'TZ58', 'ZW109', 'ZH23'};
   pcorrStrain = nan(length(strainLabels), 1);
   for strain = 1:length(strainLabels)
      load(['../../pulseShapes.results/res/' strainLabels{strain} '_pulses.mat'], 'pulseShapesNorm');  % load normalized pulses
      load(['../../pulseShapes.results/res/' strainLabels{strain} '_clustering.mat'], 'Gw');           % load original pulse labels
      strainPulseLabels = Gw==1;    % convert pulse labels such that P1=1, P2=2
      Gc = classifyPulses(double(pulseShapesNorm)/1000);   % classify pulses
      pcorrStrain(strain) = mean(Gc==strainPulseLabels);   % get fraction of correct classifications
      save('data/classifierXval', 'pcorrStrain')
   end
catch ME
   disp(ME.getReport())
   warning('Could not cross-validate (probably missing data from other strains). Loading results instead.')
   load('data/classifierXval', 'pcorrStrain')
end

%%
% QDA decision boundary
K = qda.Coeffs(1,2).Const;
L = qda.Coeffs(1,2).Linear;
Q = qda.Coeffs(1,2).Quadratic;
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + [x1, x2]*Q*[x1, x2]';

clf
subplot(1,4,1:2)
hl = gscatter(projectedPulses(:,1), projectedPulses(:,2), pulseLabels,lines(2));
hold on
h2 = ezplot(f,[-1 1 -0.53 1]);
set(h2, 'LineWidth', 2, 'Color','k');
axis('tight')
xlabel('projected onto P1 template')
ylabel('projected onto P2 template')
axis('tight', 'square')
legend([hl; h2], {'P1', 'P2', 'qda'}, 'Location', 'SouthWest')
set(gca, 'XTick', -.5:.5:1, 'YTick', -1:.5:1)
hline(0)
vline(0)

subplot(143)
hs = plotSpread([pcorr' pcorrShuffle']);
set(hs{1}, 'MarkerSize', 18, 'Color','k')
ylabel('p_{correct}')
set(gca, 'YLim',[0 1], 'XTick',1:2,'XTickLabel',{'test orginal', 'test shuffled'},...
   'XTickLabelRotation', 60)
title('classification from prj')
hline(0.5)

subplot(144)
hb1 = bar(1, pcorrStrain(1), 'EdgeColor', 'none', 'FaceColor', [1 .6 .6]);
hold on
hb2 = bar(2:length(pcorrStrain), pcorrStrain(2:end), 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6]);
hline(pcorrStrain(1))
set(gca, 'XTick', 1:length(strainLabels), 'XTickLabel', strainLabels, 'XTickLabelRotation',60, 'YLim', [0.9 1.0])
ylabel('p_{correct}')
title('test on different strains')
set(gcas, 'Color', 'none', 'Box','off','TickDir', 'out')
figexp('pulseTypeClassifierBuild', 1.5, .5)
