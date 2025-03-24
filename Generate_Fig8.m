% Generate Fig. 8 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., Mar. 2025

clear all; close all;
set_figure_style;

% Check for saved data
if exist('Fig8.mat', 'file') ~= 2
    [Af, Qfb, Qf, f, Seg1, Seg2, Seg3, Nfft] = initializeMatrices();
    save('Fig8.mat', 'Af', 'Qfb', 'Qf', 'f', 'Seg1', 'Seg2', 'Seg3', 'Nfft');
else
    load('Fig8.mat');
end

% Define parameters
ri = 2; ci = 2; 
Seg1 = (1:342); Seg2 = (343:683); Seg3 = (684:1024);
f = (0:1023) / 1024;

% Plot legend placeholders
plotLegend();

% plot target curves
plotTargetCurves(Af, f, ri, ci)

% Plot main curves
plotMainCurves(Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci);

% Axis and labels
axis([0 1 -1.1 1.1]);
set(gca, 'XTick', (0:1/6:1), 'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'});

% Legend
legend({'$A_{22}(\mathrm{e}^{\mathrm{j}\Omega})$', ...
       '$\hat{Q}_{22}^{(N \rightarrow \infty)}(\Omega)$', ...
       '$\hat{Q}_{\ast,22}(\mathrm{e}^{\mathrm{j}\Omega})$'}, ...
       'Interpreter', 'latex', 'Location', 'SouthEast');

% Labels
xlabel('normalised angular frequency $\Omega$');
ylabel('$\Re\{\cdot\}$, $\Im\{\cdot\}$');


% Plot insets
position = [0.175 0.65 0.18 0.25];
b2 = [11.5/36 12.5/36 -0.8 1];
plotInset(position, b2, Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci)
set(gca,'XTick',(11.5:.5:12.5)/36,'XTickLabel',{'$\frac{23\pi}{36}$',...
     '$\frac{2\pi}{3}$','$\frac{25\pi}{36}$'});

position = [0.71 0.65 0.18 0.25];
b1 = [23.5/36 24.5/36 -1 .2];
plotInset(position, b1, Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci)
set(gca,'XTick',(23.5:.5:24.5)/36,'XTickLabel',{'$\frac{47\pi}{36}$',...
     '$\frac{4\pi}{3}$','$\frac{49\pi}{36}$'});

% Adjust figure settings
set(gcf, 'OuterPosition', [230 250 570 350]);
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save figure
print('-depsc', 'Figures/Fig8.eps');



function [Af, Qfb, Qf, f, Seg1, Seg2, Seg3, Nfft] = initializeMatrices()
    % Matrix definition
    S = zeros(2,2,3);
    S(1,1,:) = [0 1 0];
    S(2,2,:) = [0.5 0.5 0.5];

    U = zeros(2,2,2);
    U(:,:,1) = [1 1; 0 0]/sqrt(2);
    U(:,:,2) = [0 0; 1 -1]/sqrt(2);
    
    V = zeros(2,2,1);
    V(:,:,1) = [1 -1; 1 1]/sqrt(2);

    A3 = PolyMatConv(U, PolyMatConv(S, ParaHerm(V)));

    A = zeros(2,2,4);
    A(1,1,:) = [-1 1 -1 0]/4;
    A(1,2,:) = [1 3 1 0]/4;
    A(2,1,:) = [0 1 3 1]/4;
    A(2,2,:) = [0 -1 1 -1]/4;

    % Compute DFT
    Af = fft(A,1024,3);

    % Binwise solution
    Qfb = zeros(2,2,1024);
    for i = 1:1024
        [u,~,v] = svd(Af(:,:,i));
        Qfb(:,:,i) = u*v';
    end    

    % Analytic solution
    N = 16; % Order of allpass filter
    A2 = zeros(2,2,size(A3,3) + N + 1);
    A2(:,:,N+2:end) = A3;

    B = zeros(2,2,1); 
    B(:,:,1) = eye(2);
    Qhat = PUProcrustes(A2, B, 2^12, 0, N);

    ProcrustesMetrics(A3, S, U, V, Qhat);

    Nfft = size(Qhat,3);
    Qhat2 = zeros(2,2,4*Nfft);
    Qhat2(:,:,1:Nfft) = Qhat;
    Nfft = 4*Nfft;

    Qf = fft(circshift(Qhat2, [0 0 -N-1]), Nfft, 3);

    % Frequency axis segmentation
    f = (0:1023)/1024;
    Seg1 = (1:342);
    Seg2 = (343:683);
    Seg3 = (684:1024);
end

function plotLegend()
    % Plot invisible objects for legend styling
    h = plot([-1 -1], [-1 -1], '-', 'linewidth', 3); 
    set(h(1), 'color', [1 1 1] * 0.75);
    hold on;
    plot([-1 -1], [-1 -1], 'k--');
    plot([-1 -1], [-1 -1], 'k-');
end

function plotTargetCurves(Af, f, ri, ci)
    % Define colors
    blue_shade = [.8 0.8 1];
    red_shade = [1 0.6 0.6];

    % Plot real part
    h = plot(f, real(squeeze(Af(ri, ci, :))), '-', 'linewidth', 3); 
    set(h(1), 'color', blue_shade);

    % Plot imaginary part
    h = plot(f, imag(squeeze(Af(ri, ci, :))), '-', 'linewidth', 3); 
    set(h(1), 'color', red_shade);
end

function plotMainCurves(Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci)
    grid on;

    % Plot real part
    plot((0:Nfft-1)/Nfft, real(squeeze(Qf(ri, ci, :))), 'b-');
    plotSegmentedQfb(f, real(Qfb), ri, ci, Seg1, Seg2, Seg3, 'b--');

    % Plot imaginary part
    plot((0:Nfft-1)/Nfft, imag(squeeze(Qf(ri, ci, :))), 'r-');
    plotSegmentedQfb(f, imag(Qfb), ri, ci, Seg1, Seg2, Seg3, 'r--');

    % markers
    plot(1/3, squeeze(real(Qfb(ri, ci, Seg1(end)))), 'bo', 'MarkerFaceColor', 'r','LineWidth', 1.5); 
    plot(1/3, squeeze(real(Qfb(ri, ci, Seg2(1)))), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 1.5); 
    plot(2/3, squeeze(real(Qfb(ri, ci, Seg2(end)))), 'bo', 'LineWidth', 1.5);
    plot(2/3, squeeze(real(Qfb(ri, ci, Seg3(1)))), 'ro', 'MarkerFaceColor', 'b', 'LineWidth', 1.5); 
    plot(1/3, squeeze(imag(Qfb(ri, ci, Seg2(1)))), 'ro', 'LineWidth', 1.5); 
    plot(2/3, squeeze(imag(Qfb(ri, ci, Seg2(end)))), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    
    % plot boxes
    b1 = [23.5/36 24.5/36 -1 .2];
    plotBoxes(b1);

    b2 = [11.5/36 12.5/36 -0.8 1];
    plotBoxes(b2);
end

function plotBoxes(b1)
    boxx1 = [b1(1) b1(2) b1(2) b1(1) b1(1)]; boxy1 = [b1(4) b1(4) b1(3) b1(3) b1(4)];
    plot(boxx1, boxy1, 'k');
end

function plotSegmentedQfb(f, Qfb, ri, ci, Seg1, Seg2, Seg3, lineStyle)
    % Function to plot real or imaginary parts of Qfb across different segments
    plot(f(Seg1), squeeze(Qfb(ri, ci, Seg1)), lineStyle);
    plot(f(Seg2), squeeze(Qfb(ri, ci, Seg2)), lineStyle);
    plot(f(Seg3), squeeze(Qfb(ri, ci, Seg3)), lineStyle);
end

function plotInset(position, b, Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci)
    % Create inset axis at the specified position
    axes('position', position);
    hold on;

    plotMainCurves(Qf, Qfb, f, Seg1, Seg2, Seg3, Nfft, ri, ci);
   
    axis(b)
end


