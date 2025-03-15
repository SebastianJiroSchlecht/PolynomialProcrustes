% TSP24_1_Figs2to4.m
%
%
% Generate Figs. 2, 3, and 4 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., Mar. 2025

clear all; close all;

FS = 12;

%------------------------------------------------------------------------------
%  systems and parameters
%------------------------------------------------------------------------------
S = zeros(2,2,3);
S(1,1,:) = [0 1 0];
S(2,2,:) = [.5 .5 .5];

s1 = zeros(1,1024); s1(1:2) = [1 0]; s1(end) = 0;
s2 = zeros(1,1024); s2(1:2) = [.5 .5]; s2(end) = .5;
U = zeros(2,2,2);
U(:,:,1) = [1 1; 0 0]/sqrt(2);
U(:,:,2) = [0 0; 1 -1]/sqrt(2);
V = zeros(2,2,1);
V(:,:,1) = [1 1; -1 1]/sqrt(2);
A = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

A2 = zeros(2,2,4);
A2(:,:,1) = [-1 1; 0 0];
A2(:,:,2) = [1 3;  1 -1];
A2(:,:,3) = [-1 1; 3  1];
A2(:,:,4) = [0 0;  1 -1];
A2 = A2/4;
V3 = zeros(2,2,3); V3(:,:,2) = V(:,:,1);
Ahat1 = PolyMatConv(U,ParaHerm(V3));
D = zeros(2,2,1); D(:,:,1) = diag([1,-1]);
Ahat2 = PolyMatConv(U,PolyMatConv(D,ParaHerm(V3)));


plotFourierApproximation(s1, s2);


%------------------------------------------------------------------------------
%  Procrustes attempt
%------------------------------------------------------------------------------
figure(2);
Nfft = 2^8;
IndexN = (-Nfft/2:(Nfft/2-1));
Vprime = zeros(2,2,Nfft);
Vprime(:,1,Nfft/2+1) = [1; -1]/sqrt(2);
% Fourier series of square wave
Vprime(:,2,1:Nfft/2) = [1; 1] * sqrt(2)./(IndexN(1:Nfft/2)*pi).*sin(IndexN(1:Nfft/2)*pi*2/3);
Vprime(:,2,Nfft/2+2:end) = [1; 1] * sqrt(2)./(IndexN(Nfft/2+2:end)*pi).*sin(IndexN(Nfft/2+2:end)*2*pi/3);
Vprime(:,2,Nfft/2+1) = [1; 1]*(1/3)/sqrt(2);

 
% attempted Procrustes solution
Ahat3 = PolyMatConv(U,ParaHerm(Vprime));
Ahat3 = real(Ahat3(:,:,117:142));
Anew = zeros(2,2,26);
Anew(:,:,11:14) = A;

error3 = sum(sum(sum(abs(Anew-Ahat3).^2)));

plotProcrustes(Anew, Ahat3);


%------------------------------------------------------------------------------
%  error metrics
%------------------------------------------------------------------------------
PowTwos = (3:18);
ErrorNorm = zeros(1,length(PowTwos));
PUError = zeros(1,length(PowTwos));
for i = 1:length(PowTwos),
  Nfft = 2^PowTwos(i);
  IndexN = (-Nfft/2:(Nfft/2-1));
  Vprime = zeros(2,2,Nfft);
  Vprime(:,1,Nfft/2+1) = [1; -1]/sqrt(2);
  % Fourier series of square wave
  Vprime(:,2,1:Nfft/2) = [1; 1] * sqrt(2)./(IndexN(1:Nfft/2)*pi).*sin(IndexN(1:Nfft/2)*pi*2/3);
  Vprime(:,2,Nfft/2+2:end) = [1; 1] * sqrt(2)./(IndexN(Nfft/2+2:end)*pi).*sin(IndexN(Nfft/2+2:end)*2*pi/3);
  Vprime(:,2,Nfft/2+1) = [1; 1]*(1/3)/sqrt(2);

  % solution in the time domain
  Ahat3 = PolyMatConv(U,ParaHerm(Vprime)); 
  % mismatch w.r.t. A
  Error = Ahat3;
  t = (Nfft/2-1:Nfft/2+2);
  Error(:,:,t) = Error(:,:,t) - A;
%  ErrorNorm(i) = sum(sum(sum(abs(Error))));
  for n = 1:size(Error,3),
     %
     %   changed |.|_F to |.|^2_F, SW, 5/1/24
     % 
     ErrorNorm(i) = ErrorNorm(i) + norm(Error(:,:,n),'fro')^2;
  end;   
  % mismatch w.r.t. paraunitarity
  Qf = fft(Ahat3,2*Nfft,3);
  for k = 1:2*Nfft,
     PUError(i) = PUError(i) + sum(sum(abs(squeeze(Qf(:,:,k))*squeeze(Qf(:,:,k))'-eye(2)).^2));
  end;   
  PUError(i) = PUError(i)/(2*Nfft);
end;

plotErrorMetrics(PowTwos, ErrorNorm, PUError);




function plotErrorMetrics(PowTwos, ErrorNorm, PUError)
    % Function to plot Least Squares (LS) mismatch and Paraunitarity (PU) error
    
    FS = get(0, 'DefaultAxesFontSize') + 1; % Increase font size for better readability

    % Create figure
    fig = figure(3); clf;
    
    % Define colors for the two y-axes
    left_color = [1 0 0];  % Red for Least Squares Error
    right_color = [0 0 1]; % Blue for Paraunitarity Error
    
    set(fig, 'defaultAxesColorOrder', [left_color; right_color]);

    % Hold for multiple plots
    hold on; grid on;

    % Legend placeholders
    plot([-1 -1], [1 2], 'b-'); % Blue line for PU Error
    plot([-1 -1], [1 2], 'r--'); % Red dashed line for LS Error

    % Limiting value for LS error
    h = plot([2 18], [1 1] * (17/12 - 2 * sqrt(3) / pi), '-', 'LineWidth', 2);
    set(h, 'Color', [0.7 0.7 0.7]); 

    % Plot Least Squares Error (LS Error)
    plot(PowTwos, ErrorNorm, 'r--');
    
    % Configure left Y-axis (LS Error)
    axis([2 18 0.2 0.325]);
    set(gca, 'YTick', (0.2:0.025:0.325), ...
             'YTickLabel', {'$0.200$', '$0.225$', '$0.250$', '$0.275$', '$0.300$', '$0.325$'});

    % Left Y-axis Label (Least Squares Mismatch)
    ylabel('least squares mismatch, $e_{\mathrm{LS}}$', ...
          'FontSize', FS);
    
    % Switch to right Y-axis
    yyaxis right;
    
    % Plot Paraunitarity Error
    plot(PowTwos, 10 * log10(PUError), 'b-');
    
    % Configure right Y-axis (PU Error)
    axis([2 18 -55 -5]);
    set(gca, 'YTick', (-55:10:-5), ...
             'YTickLabel', {'$-55$', '$-45$', '$-35$', '$-25$', '$-15$', '$-5$'});

    % Right Y-axis Label (Paraunitarity Error)
    ylabel('paraunitarity error, $10\mathrm{log}_{10}\{ e_{\mathrm{PU}} \}$', ...
           'FontSize', FS);

    % X-axis Label
    xlabel('$\mathrm{log}_2 N$', 'FontSize', FS);

    % Configure legend
    legend({'$e_{\mathrm{PU}}$', '$e_{\mathrm{LS}}$', ...
            '$\mathrm{lim}_{N\rightarrow \infty}(e_{\mathrm{LS}})$'}, ...
            'FontSize', FS-2, 'Location', 'East');

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 350]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', 'Fig4.eps');
end


function plotProcrustes(Anew, Ahat3)

% Time range and font size adjustment
t = (-10:15);

% Create Figure
figure(2); clf;

% Plot subplots
plotSubplotProcrustes(1, t, Anew, Ahat3, 1, 1, ...
    '$a_{1,m}[n], \; \hat{q}^{(N)}_{1,m}[n]$', '', 1);

plotSubplotProcrustes(2, t, Anew, Ahat3, 1, 2, ...
    '', '', 2);

plotSubplotProcrustes(3, t, Anew, Ahat3, 2, 1, ...
    '$a_{2,m}[n], \; \hat{q}^{(N)}_{2,m}[n]$', 'time index $n$', 1);

plotSubplotProcrustes(4, t, Anew, Ahat3, 2, 2, ...
    '', 'time index $n$', 2);

% Adjust figure layout
set(gcf, 'OuterPosition', [230 250 570 350]);
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save figure
print('-depsc', 'Fig3.eps');

end

function plotSubplotProcrustes(position, t, Anew, Ahat3, row, col, label_y, label_x, m_value)
    % Function to plot a single subplot for Procrustes results
    
    subplot(2,2,position);
    stem(t, squeeze(Anew(row, col, :)), 'b'); hold on;
    plot(t, squeeze(Ahat3(row, col, :)), 'r*');
    
    legendFontSize = get(0, 'DefaultAxesFontSize') - 2;

    % Axis limits and grid
    axis([-10.5 10.5 -1.1 1.1]); 
    grid on;
    
    % Tick labels
    set(gca, 'XTick', -10:5:10, ...
             'XTickLabel', {'$-10$', '$-5$', '$0$', '$5$', '$10$'});

    % Optional Y-label
    if ~isempty(label_y)
        ylabel(label_y, 'FontSize', legendFontSize);
    end

    % Optional X-label
    if ~isempty(label_x)
        xlabel(label_x, 'FontSize', legendFontSize);
    end

    % Text annotation
    text(-9, 0.75, ['$m=' num2str(m_value) '$'], 'FontSize', legendFontSize);
end



function plotFourierApproximation(s1, s2)
    % Function to plot the Fourier series approximation with non-negative singular values
    figure(1); clf;
    
    % Compute Fourier Transform
    freq_axis = (0:1023) / 1024;
    fft_s1 = real(fft(s1, 1024));
    fft_s2 = abs(real(fft(s2, 1024)));

    % Plot signals
    plot(freq_axis, fft_s1, 'b-');
    hold on;
    plot(freq_axis, fft_s2, 'r--');
    
    % Configure axis limits and grid
    axis([0 1 0 1.5]);
    grid on;

    % Labels
    ylabel('$|\sigma_m(\mathrm{e}^{\mathrm{j}\Omega})|$');
    
    % Configure Tick Labels
    set(gca, 'XTick', (0:1/6:1), ...
        'XTickLabel', {'$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'}, ...
        'YTick', (0:0.5:1.5), ...
        'YTickLabel', {'$0$', '$.5$', '$1$', '$1.5$'});

    % Add legend
    legendFontSize = get(0, 'DefaultAxesFontSize') - 2;
    legend({'$m=1$', '$m=2$'},'FontSize', legendFontSize, 'Location', 'SouthWest');

    % X-axis Label
    xlabel('Normalised angular frequency $\Omega$');

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 200]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', 'Fig2.eps');
end