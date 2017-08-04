load('data/NM91_pulses.mat', 'pulseShapes');    % load pulses
pulses = single(pulseShapes)/1000;              % saved as int16 - cast to single and rescale
pulsesNorm = normalizePulses(pulses);           % normalize pulses
pulseLabels = classifyPulses(pulsesNorm);       % classify pulses

%% plot results
pulseIdx = uint16(linspace(1, size(pulses,1), 100));  % select subset of pulses for plotting
clf
subplot(311)
plot(pulses(pulseIdx,:)', 'Color', [.4 .4 .4])
title('raw pulses')
subplot(312)
plot(pulsesNorm(pulseIdx,:)', 'Color', [.4 .4 .4])
title('normalized pulses')
subplot(313)
hold on
h2 = plot(pulsesNorm(pulseIdx(pulseLabels(pulseIdx)==1),:)', 'r');
h1 = plot(pulsesNorm(pulseIdx(pulseLabels(pulseIdx)==0),:)', 'k');
title('labeled pulses')
legend([h1(1), h2(1)], {'P_{slow}', 'P_{fast}'})
legend('boxoff')
set(gcas, 'Color', 'none', 'Box','off','TickDir', 'out')
axis(gcas,'tight', 'off')
figexp('demo', 0.6, 1)
