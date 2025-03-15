% TSP24_1_Figs9to12.m
%
% Generate Figs. 9, 10, 11, and 12 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., Mar. 2025
 
clear all; close all;

M = 3; N = 10; FS = 12;

%------------------------------------------------------------------------------
%  matrix A
%------------------------------------------------------------------------------
% singular values
Nfft = 64;
s = [0  1 0  1 0;
     0 .5 0 .5 0; 
     -1i*.5 .5 0 .5 1i*.5];
st = zeros(3,Nfft); st(:,1:3) = s(:,3:5); st(:,Nfft-1:Nfft) = s(:,1:2);     
sf = fft(st,Nfft,2);
% figure(1); plot((0:Nfft-1)/Nfft,real(sf)','*-'); 

U = PUPolyMatRand(3,10,0,'complex');
V = PUPolyMatRand(3,10,1,'complex');
S = zeros(3,3,5);
for m = 1:3, S(m,m,:) = s(m,:); end;
A2 = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

B = zeros(M,M,1); B(:,:,1) = eye(M);
A = zeros(M,M,size(A2,3)+N+1);
A(:,:,N+2:end) = A2;
if exist('Qhat1024.mat','file'),
   load Qhat2048
else    
   Qhat = PUProcrustes(A,B,2^12,0,15);                           % Procrustes 
end;

%-----------------------------------------------------------------------------
%  Figure 9: singular values
%-----------------------------------------------------------------------------

% Define FFT parameters
Nfft2 = 2048; 
ss = zeros(3, Nfft2); 
ss(:,1:3) = s(:,3:5); 
ss(:,Nfft2-1:Nfft2) = s(:,1:2); 
Sf = real(fft(ss, Nfft2, 2));
Sf = [Sf Sf(:,1)]; % Append first column for wrap-around effect

plotFourierTransform(Sf, 'Fig9.eps');

%-----------------------------------------------------------------------------
%  Figure 10 and 11: polynomial Procrustes solution
%-----------------------------------------------------------------------------

plotMatrixSubplots(real(A2), real(Qhat), 'Fig10.eps');
plotMatrixSubplots(imag(A2), imag(Qhat), 'Fig11.eps');
 
%-----------------------------------------------------------------------------
%  Figure 12: \Pi(z) \Sigma(z)
%-----------------------------------------------------------------------------
P = PolyMatConv(ParaHerm(U),PolyMatConv(Qhat,V));
Ls = size(S,3);
Pi = zeros(2048,3);
for i = 1:3, 
   dummy = squeeze(P(i,i,34:end));
   Pi(1:length(dummy),i) = dummy;  
   Pi(2048-32:end,i) = squeeze(P(i,i,1:33)); 
end;
Pif = fft(Pi,2048,1);
Pif = [Pif; Pif(1,:)];
SigmaPif = real(Pif.*Sf.');

% metric --- diagonalisation
SSS = PolyMatConv(S,P);
diagon = PolyMatNorm(SSS,'OffDiag')/PolyMatNorm(SSS,'OnDiag');
% metric --- positivity
positivity = norm(SigmaPif-abs(Sf.'),'fro').^2/norm(abs(Sf),'fro').^2;
% metrix --- paraunitarity
QQ = PolyMatConv(Qhat,ParaHerm(Qhat));
Lq = size(QQ,3);
QQ(:,:,(Lq+1)/2) = QQ(:,:,(Lq+1)/2) - eye(3);
paraunitarity = PolyMatNorm(QQ);

plotSpectralMagnitude(Nfft2, Sf, SigmaPif, 'Fig12.eps');

     
%-----------------------------------------------------------------------------
%  some numerical evaluations
%-----------------------------------------------------------------------------
Res = ProcrustesMetrics(A,S,U,V,Qhat);   
disp(sprintf('error in paraunitarity:     %2.12g',Res(1)));
disp(sprintf('error in diagonalisation:   %2.12g',Res(2)));
disp(sprintf('error in positivity:        %2.12g',Res(3)));
disp(sprintf('least squares error A-Q:    %2.12g',Res(4)));
[~,dummy,~,~] = PolyMatAlign(A,PolyMatConv(U,ParaHerm(V)));
disp(sprintf('least squares error A-UV^P: %2.12g',dummy));



function plotSpectralMagnitude(Nfft2, Sf, SigmaPif, figName)
    % Function to plot spectral magnitude |σ_m(e^jΩ)| and real part of s_m,m(e^jΩ)
    % including a zoomed-in inset
    
    figure; clf;

    % Create placeholder plots for legend
    plot([-1 -2], [0 0], 'b-'); hold on;
    plot([-1 -2], [0 0], 'r--');
    plot([-1 -2], [0 0], '-.', 'Color', [0 0.5 0]);
    
    plotMagnitudes(Sf, SigmaPif, Nfft2);

    % Axis and grid settings
    axis([0 1 0 2]); 
    grid on;

    % Tick Labels
    set(gca, ...
        'XTick', (0:1/8:1), ...
        'XTickLabel', {'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', ...
                       '$3\pi/2$', '$7\pi/4$', '$2\pi$'}, ...
        'YTick', (0:1:2), ...
        'YTickLabel', {'$0$', '$1$', '$2$'});

    % Labels
    xlabel('normalised angular frequency $\Omega$');
    ylabel('$|\sigma_{m}(\mathrm{e}^{\mathrm{j}\Omega})|$, $\Re\{s_{m,m}(\mathrm{e}^{\mathrm{j}\Omega})\}$');

    % Bounding box for inset
    Box1 = [0.245 0.255 0 0.05];
    plot(Box1([1 1 2 2 1]), Box1([3 4 4 3 3]), 'k-', 'LineWidth', 1);

    % Legend
    legend({'$m=1$', '$m=2$', '$m=3$'}, ...
           'FontSize', 10, 'Location', 'NorthEast');

    % Create inset plot
    axes('Position', [0.29 0.64 0.07 0.25]);
    hold on;
    plotMagnitudes(Sf, SigmaPif, Nfft2);

    % Axis settings for inset
    axis([0.245 0.255 0 0.05]); 
    grid on;

    % Tick Labels for inset
    set(gca, ...
        'XTick', [0.245 0.255], ...
        'XTickLabel', {'$\frac{49\pi}{100}$', '$\frac{51\pi}{100}$'}, ...
        'YTick', [0 0.05], ...
        'YTickLabel', {'$0.00$', '$0.05$'});

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 280]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end

function plotMagnitudes(Sf, SigmaPif, Nfft2)

% Plot absolute values of spectral magnitudes
freq_axis = (0:Nfft2) / Nfft2;
for i = 1:3
    h = plot(freq_axis, abs(Sf(i,:)), '-', 'LineWidth', 3);
    set(h(1), 'Color', [1 1 1] * 0.7);
end

% Plot processed spectral functions
plot(freq_axis, SigmaPif(:,1), 'b-');
plot(freq_axis, SigmaPif(:,2), 'r--');
h = plot(freq_axis, SigmaPif(:,3), '-.');
set(h(1), 'Color', [0 0.5 0]);

end


function plotMatrixSubplots(A2, Qhat, figName)
    % Function to plot 3x3 matrix subplots of A2 vs. Qhat (real or imaginary parts)
    figure; clf;

    % Iterate over matrix indices
    for m = 1:3
        for n = 1:3
            subplot(3, 3, 3*(m-1) + n); 

            % Plot stem plot for A2
            stem(0:24, squeeze((A2(m, n, :))), 'b'); 
            hold on;

            % Plot Qhat values
            plot(0:24, squeeze((Qhat(m, n, 12:36))), 'r*'); 

            % Configure axes
            axis([0 24 -0.4 0.4]); 
            grid on;

            % Set Y-label for first column
            if n == 1
                ylabel(sprintf('$\\ell=%d$', m));
            end  

            % Set Title for first row
            if m == 1
                title(sprintf('$m=%d$', n));
            end

            % Set X-label for last row
            if m == 3
                xlabel('Index $n$');
            end

            % Configure Tick Labels
            set(gca, ...
                'XTick', (0:8:24), ...
                'XTickLabel', {'$0$', '$8$', '$16$', '$24$'}, ...
                'YTick', (-0.4:0.2:0.4), ...
                'YTickLabel', {'$-.4$', '$-.2$', '$0$', '$.2$', '$4$'});
        end
    end

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 400]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end


function plotFourierTransform(Sf, figName)
    % Function to plot the Fourier Transform of signal s with sampling points
    Nfft2 = size(Sf,2) - 1;
    Nfft3 = 16;

    % Create figure
    figure; clf;

    % Plot Fourier Transforms
    freq_axis = (0:Nfft2) / Nfft2;
    plot(freq_axis, Sf(1,:), 'b-'); hold on;
    plot(freq_axis, Sf(2,:), 'r--');
    h = plot(freq_axis, Sf(3,:), '-.');
    set(h(1), 'Color', [0 0.5 0]); % Set color for m=3

    % Plot sampled points
    sample_indices = 1:Nfft2/Nfft3:Nfft2;
    plot(freq_axis(sample_indices), Sf(1, sample_indices), 'b*');
    plot(freq_axis(sample_indices), Sf(2, sample_indices), 'r*');
    h = plot(freq_axis(sample_indices), Sf(3, sample_indices), '*');
    set(h(1), 'Color', [0 0.5 0]);

    % Axis settings
    axis([0 1 -2 2]); 
    grid on;

    % Plot markers
    plot([0 1 2 3 4]/4, [1 0 -1 0 1], 'ko', 'MarkerFaceColor', 'k');
    plot([0.08333449 .5-0.08333449], [1 -1]*1.732036, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'w');

    % Tick Labels
    set(gca, ...
        'XTick', (0:1/8:1), ...
        'XTickLabel', {'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', ...
                       '$3\pi/2$', '$7\pi/4$', '$2\pi$'}, ...
        'YTick', -2:1:2, ...
        'YTickLabel', {'$-2$', '$-1$', '$0$', '$1$', '$2$'});

    % Labels
    xlabel('normalised angular frequency $\Omega$');
    ylabel('$\sigma_{m}(\mathrm{e}^{\mathrm{j}\Omega})$');

    % Legend
    legend({'$m=1$', '$m=2$', '$m=3$'}, ...
           'FontSize', 10, 'Location', 'SouthWest');	

    % Adjust figure layout
    set(gcf, 'OuterPosition', [230 250 570 280]);
    set(gca, 'LooseInset', get(gca, 'TightInset'));

    % Save figure
    print('-depsc', figName);
end
