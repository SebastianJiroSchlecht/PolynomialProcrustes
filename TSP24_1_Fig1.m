% TSP24_1_Fig1.m
%
% Generate Fig. 1 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., Mar. 2025
%
% Example for an analytic SVD with a negative singular value

clear all; close all;
set_figure_style;

% singular values
S = zeros(2,2,3);
S(1,1,:) = [0 1 0];
S(2,2,:) = [.5 .5 .5];

%------------------------------------------------------------------------------
%  Figure: singular values, Fourier domain
%------------------------------------------------------------------------------
FS = 12; % Font size for labels

% Define signals
s1 = zeros(1, 1024); s1([end,1,2]) = S(1,1,:);
s2 = zeros(1, 1024); s2([end,1,2]) = S(2,2,:);

% Call the function to plot FFT results
plotFFT(s1, s2, 'Fig1.eps');


% Compute the matrices
U = zeros(2,2,2);
U(:,:,1) = [1 1; 0 0]/sqrt(2);
U(:,:,2) = [0 0; 1 -1]/sqrt(2);
V = zeros(2,2,1);
V(:,:,1) = [1 -1; 1 1]/sqrt(2);
A = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

A2 = zeros(2,2,4);
A2(:,:,1) = [1 1; 0 0];
A2(:,:,2) = [3 -1; -1 -1];
A2(:,:,3) = [1 1;   1  -3];
A2(:,:,4) = [0 0;  -1 -1];
A2 = A2/4;
V3 = zeros(2,2,3); V3(:,:,2) = V(:,:,1);
Ahat1 = PolyMatConv(U,ParaHerm(V3));
D = zeros(2,2,1); D(:,:,1) = diag([1,1]);
Ahat2 = PolyMatConv(U,PolyMatConv(D,ParaHerm(V3)));


error1 = sum(sum(sum(abs(A-Ahat1).^2)))
error2 = sum(sum(sum(abs(A-Ahat2).^2)))

%------------------------------------------------------------------------------
%  Figure: Procrustes solution, time domain
%------------------------------------------------------------------------------

figure(2);
plotSubplot(1, A(1,1,:), Ahat1(1,1,:), '$a_{1,m}[n], \; q_{1,m}[n]$', '', 1);
plotSubplot(2, A(1,2,:), Ahat1(1,2,:), '', '', 2);
plotSubplot(3, A(2,1,:), Ahat1(2,1,:), '$a_{2,m}[n], \; q_{2,m}[n]$', 'time index $n$', 1);
plotSubplot(4, A(2,2,:), Ahat1(2,2,:), '', 'time index $n$', 2);

set(gcf, 'OuterPosition', [230 250 300 400]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
% print -depsc Fig2a.eps

figure(3);
plotSubplot(1, A(1,1,:), Ahat2(1,1,:), '$a_{1,m}[n], \; q_{\ast,1,m}[n]$', '', 1);
plotSubplot(2, A(1,2,:), Ahat2(1,2,:), '', '', 2);
plotSubplot(3, A(2,1,:), Ahat2(2,1,:), '$a_{2,m}[n], \; q_{\ast,2,m}[n]$', 'time index $n$', 1);
plotSubplot(4, A(2,2,:), Ahat2(2,2,:), '', 'time index $n$', 2);

set(gcf, 'OuterPosition', [230 250 300 400]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
% print -depsc Fig2b.eps




function plotFFT(s1, s2, figName)
    figure;
    
    % Compute FFT
    freq_axis = (0:1023) / 1024;
    fft_s1 = real(fft(s1, 1024));
    fft_s2 = real(fft(s2, 1024));

    % Plot the FFT results
    plot(freq_axis, fft_s1, 'b-'); hold on;
    plot(freq_axis, fft_s2, 'r--');

    % Set axis limits and grid
    axis([0 1 -0.5 1.5]); grid on;

    % Y-axis label
    ylabel('$\sigma_m(\mathrm{e}^{\mathrm{j}\Omega})$');

    % Configure tick labels
    set(gca, 'TickLabelInterpreter', 'latex', ...
        'XTick', (0:1/6:1), 'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'}, ...
        'YTick', (-1.5:0.5:1.5), 'YTickLabel', {'$-1.5$', '$-1$', '$-.5$', '$0$', '$.5$', '$1$', '$1.5$'});

    % Add legend
    legendFontSize = get(0, 'DefaultAxesFontSize') - 2;
    legend({'$m=1$', '$m=2$'},  'fontsize', legendFontSize, 'location', 'West');

    % X-axis label
    xlabel('normalised angular frequency $\Omega$');

    % Adjust figure size and insets
    set(gcf, 'OuterPosition', [230 250 570 250]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
    
end



function plotSubplot(pos, A, Ahat, ylabel_text, xlabel_text, m)
    % This function plots a single subplot with stem and plot

    subplot(2,2,pos); % Position in the subplot grid
    stem(0:3, squeeze(A), 'b'); hold on; % Plot stem plot
    plot(0:3, squeeze(Ahat), 'r*'); % Overlay red stars
    
    % Set axis limits and grid
    axis([-.5 3.5 -1.1 1.1]); grid on;
    
    % Add text annotation for m value
    text(0, 0.75, ['$m=' num2str(m) '$']);
    
    % Add ylabel if provided
    if ~isempty(ylabel_text)
        ylabel(ylabel_text);
    end
    
    % Add xlabel if provided
    if ~isempty(xlabel_text)
        xlabel(xlabel_text);
    end
end