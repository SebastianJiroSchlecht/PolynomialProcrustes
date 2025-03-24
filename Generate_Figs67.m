% Generate Figs. 6 and 7 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., Mar. 2025

close all; clear all;

Nfft = 2^10;
W = (0:Nfft-1)/Nfft;
Sigma = .5 + .5*cos(2*pi*W) + sin(4*pi*W);

% approximate zero crossings:
Wo = [0.2826 0.5 0.809924 0.907476]*2*pi;

%------------------------------------------------------------------------------
%  generate and plot switching functions
%------------------------------------------------------------------------------
%  *** higher order design (Markus Lang)
NN = 2.^(1:4);
H2 = zeros(Nfft,length(NN));
Angle2 = zeros(Nfft,length(NN));
for n = 1:length(NN),
   [b,a] = SingularValueAllpassSwitchCompact(Wo,NN(n));
   h = impz(b,a,Nfft);
   h2 = circshift(h,-NN(n)+2);
   Angle2(:,n) = unwrap(angle(fft(h2,Nfft)));
   H2(:,n) = cos(unwrap(angle(fft(h2,Nfft))));
end;


%------------------------------------------------------------------------------
%  Figure 6a
%------------------------------------------------------------------------------
FS=12;

Angle0 = zeros(size(W));
for i =1:4,
    Angle0 = Angle0 - (W>Wo(i)/2/pi);
end;

plotPhaseResponse(W, Angle0, Angle2, 'Figures/Fig6a.eps');
 
%------------------------------------------------------------------------------
%  Figure 6b
%------------------------------------------------------------------------------

plotFrequencyResponse(W, Sigma, H2, 'Figures/Fig6b.eps')


%------------------------------------------------------------------------------
%  Figure 6c
%------------------------------------------------------------------------------
Ssum = sum(abs(Sigma));
Int = zeros(1,size(H2,2));
for n = 1:size(H2,2),
   H2(:,n) = H2(:,n).*Sigma.';
   Int(n) = (Ssum-sum(H2(:,n)))/Ssum;
end;

plotMagnitudeResponse(W, Sigma, H2, Wo, 'Figures/Fig6c.eps');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Figure 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------
%  find zero crossings
%---------------------------------------------------
W0 = [pi/2    2*pi/3;
      7/8*pi  9/8*pi;
      4*pi/3  5*pi/3;
      5*pi/3  2*pi];
Theta = zeros(4,1);
for i = 1:4,
   Wlower = W0(i,1); Wupper = W0(i,2);
   flower = .5 + .5*cos(Wlower) + sin(2*Wlower); fupper =   .5 + .5*cos(Wupper) + sin(2*Wupper);
   fmid = 1;
   while (abs(fmid)>1e-12);
      Wmid = (abs(flower)*Wlower + abs(fupper)*Wupper)/(abs(flower) + abs(fupper));
      fmid =   .5 + .5*cos(Wmid) + sin(2*Wmid);
      if sign(fmid) == sign(flower),
         flower = fmid; Wlower = Wmid;
      else
         fupper = fmid; Wupper = Wmid;
      end;
   end;         
   Theta(i) = Wmid;
end;

%----------------------------------------------------
%  calculate allpasses
%----------------------------------------------------
Nfft = 2^14;
NN = 2.^(1:8);
H2 = zeros(Nfft,length(NN));
Angle2 = zeros(Nfft,length(NN));
for n = 1:length(NN),
    disp(sprintf('allpass design for order %d',NN(n)));
   [b,a] = SingularValueAllpassSwitchCompact(Theta,NN(n));
   h = impz(b,a,Nfft);
   h2 = circshift(h,-NN(n)+2);
   Angle2(:,n) = unwrap(angle(fft(h2,Nfft)));
   H2(:,n) = cos(unwrap(angle(fft(h2,Nfft))));
%   plot(w,cos(unwrap(angle(fft(h2,1024)))));
end;

%----------------------------------------------------
%  numerical integration
%----------------------------------------------------
W = (0:(Nfft-1))'/Nfft*2*pi;
f = .5 + .5*cos(W) + sin(2*W);
F = zeros(length(NN),1);
for n = 1:length(NN),
   F(n) = sum( (abs(H2(:,n).*f)-1).^2 )/Nfft;  
end;
Flimit = sum( (abs(f)-1).^2 )/Nfft;

plotLogErrorScaling(NN, F, Flimit, 'Figures/Fig7.eps');  

  

%-------------------------------------------------------------------------
%    cost function
%-------------------------------------------------------------------------
function costval = CostFunction2(X)
rho=X(1); phi=X(2); theta=X(3); j = sqrt(-1);     % 1st allpass
b1 = [rho*exp(j*phi) -1]*exp(j*theta);
a1 = [1 -rho*exp(-j*phi)]; 
rho=X(4); phi=X(5); theta=X(6);                   % 2nd allpass
b2 = [rho*exp(j*phi) -1]*exp(j*theta);
a2 = [1 -rho*exp(-j*phi)]; 
b = conv(b1,b2); a = conv(a1,a2);
h = impz(b,a,512);
f = zeros(1023,1);
f(512:1023) = h;
f = f + conj(flipud(f));
f2 = zeros(1024,1);
f2(1:512) = f(512:1023);
f2(514:1024) = f(1:511);
F = fft(f2,1024);

W = (0:1023)'/1024;
S = .5 + .5*cos(2*pi*W) + sin(4*pi*W);
costval = 1/sum(real(F.*S));
end



function plotPhaseResponse(W, Angle0, Angle2, figName)
    % Function to plot phase response curves
    figure(1); clf;
    
    % Plot ideal phase response
    h = plot(W, Angle0, '-', 'LineWidth', 3); 
    set(h(1), 'Color', [1 1 1] * 0.75);
    hold on;
    
    % Plot phase responses for different J values
    plot(W, Angle2(:,1) / pi, 'b-'); % J=2
    plot(W, Angle2(:,2) / pi, 'r--'); % J=4
    h = plot(W, Angle2(:,3) / pi, '-.'); 
    set(h(1), 'Color', [0 0.5 0]); % J=8 (green)
    plot(W, Angle2(:,4) / pi, 'k:'); % J=16

    % Axis settings
    axis([0 1 -4.2 0.2]);
    grid on;

    % Labels
    text(0.01, -0.5, '(a)');
    ylabel('Phase Response $\Phi(\Omega)$');

    % Tick Labels
    set(gca, ...
        'XTick', (0:1/6:1), ...
        'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'}, ...
        'YTick', -4:1:0, ...
        'YTickLabel', {'$-4\pi$', '$-3\pi$', '$-2\pi$', '$-\pi$', '$0$'});

    % Legend
    legendFontSize = 10;
    legend({'Ideal', '$J=2$', '$J=4$', '$J=8$', '$J=16$'}, ...
           'FontSize', legendFontSize, 'Location', 'SouthWest');

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 250]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end

function plotFrequencyResponse(W, Sigma, H2, figName)
    % Function to plot frequency response curves

    figure(2); clf;
    
    % Plot ideal frequency response
    h = plot(W, sign(Sigma), '-', 'LineWidth', 3); 
    set(h(1), 'Color', [1 1 1] * 0.75);
    hold on;
    
    % Plot frequency response curves for different J values
    plot(W, H2(:,1), 'b-');  % J=2
    plot(W, H2(:,2), 'r--'); % J=4
    h = plot(W, H2(:,3), '-.'); 
    set(h(1), 'Color', [0 0.5 0]); % J=8 (green)
    plot(W, H2(:,4), 'k:');  % J=16

    % Axis settings
    axis([0 1 -1.1 1.1]);
    grid on;

    % Labels
    text(0.01, 0.8, '(b)');
    ylabel('frequency response $f(\mathrm{e}^{\mathrm{j}\Omega})$');

    % Tick Labels
    set(gca, ...
        'XTick', (0:1/6:1), ...
        'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'}, ...
        'YTick', -1:0.5:1, ...
        'YTickLabel', {'$-1$', '$-.5$', '$0$', '$.5$', '$1$'});

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 250]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end


function plotMagnitudeResponse(W, Sigma, H2, Wo, figName)
    % Function to plot magnitude response f(e^jΩ)σ(e^jΩ)

    figure(3); clf;

    % Plot ideal magnitude response
    h = plot(W, abs(Sigma), '-', 'LineWidth', 3);
    set(h(1), 'Color', [1 1 1] * 0.75);
    hold on;

    % Plot magnitude response curves for different J values
    plot(W, H2(:,1), 'b-');  % J=2
    plot(W, H2(:,2), 'r--'); % J=4
    h = plot(W, H2(:,3), '-.');
    set(h(1), 'Color', [0 0.5 0]); % J=8 (green)
    plot(W, H2(:,4), 'k:');  % J=16

    % Plot vertical dashed lines at Wo frequencies
    for n = 1:4
        plot(Wo(n) / (2 * pi) * [1 1], [-0.5 1.5], 'k--');
    end

    % Annotate frequency markers
    for n = 1:4
        text(Wo(n) / (2 * pi) - 0.01, 1.75, ['$\Omega_' num2str(n) '$']);
    end

    % Axis settings
    axis([0 1 -0.5 2]);
    grid on;

    % Labels
    text(0.01, -0.25, '(c)');
    ylabel('$f(\mathrm{e}^{\mathrm{j}\Omega}) \sigma(\mathrm{e}^{\mathrm{j}\Omega})$');

    % Tick Labels
    set(gca, ...
        'XTick', (0:1/6:1), ...
        'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'}, ...
        'YTick', -1:0.5:2, ...
        'YTickLabel', {'$-1$', '$-.5$', '$0$', '$.5$', '$1$', '$1.5$', '$2$'});

    % X-axis label
    xlabel('normalized angular frequency $\Omega$');

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 260]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end


function plotLogErrorScaling(NN, F, Flimit, figName)
    % Function to plot log-error scaling with power law reference
    
    figure; clf;
    
    % Dummy plots for legend
    plot(-1, -1, 'b-'); 
    hold on;
    h = plot(-1, -1, '-', 'LineWidth', 3);
    set(h(1), 'Color', [1 1 1] * 0.75);
    
    % Actual plots
    h = plot(log2(NN), -1 - 2 * log2(NN), '-', 'LineWidth', 3); 
    set(h(1), 'Color', [1 1 1] * 0.75);
    plot((1:8), log2(F - Flimit), 'b-');

    % Axis settings
    axis([1 8 -18 -2]);
    grid on;

    % Labels
    xlabel('$\log_{2} J$');
    ylabel('$\log_{2} (e_{\mathrm{LS}}-\epsilon_{\mathrm{LS}})$');

    % Legend
    legend({'$e_{\mathrm{LS}}-\epsilon_{\mathrm{LS}}$', 'power law $1/J^2$'}, ...
           'FontSize', 10, 'Location', 'SouthWest');

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 240]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end