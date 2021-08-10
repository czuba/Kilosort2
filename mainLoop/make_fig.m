function make_fig(W, U, mu, nsp)

% xy limits
ampLim = [0 1.1*gather(max(mu(:)))];
nspLim = [0, 1.1*gather(max(nsp(:)))];

subplot(2,2,1)
imagesc(W(:,:,1))
title('Temporal Components')
xlabel('Unit number');
ylabel('Time (samples)');

subplot(2,2,2)
imagesc(U(:,:,1))
title('Spatial Components')
xlabel('Unit number');
ylabel('Channel number');

subplot(2,2,3)
plot(mu)
ylim(ampLim)
title('Unit Amplitudes')
xlabel('Unit number');
ylabel('Amplitude (arb. units)');

subplot(2,2,4)
semilogx(1+nsp, mu, '.')
ylim(ampLim)
xlim(nspLim)
title('Amplitude vs. Spike Count')
xlabel('Spike Count');
ylabel('Amplitude (arb. units)');
drawnow