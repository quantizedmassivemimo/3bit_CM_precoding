% =========================================================================
% -- Simulator for 1-bit Massive MU-MIMO Precoding in VLSI with C3PO
% -------------------------------------------------------------------------
% -- (c) 2018 Christoph Studer, Oscar Castañeda, and Sven Jacobsson
% -- e-mail: studer@cornell.edu, oc66@cornell.edu, and
% -- sven.jacobsson@ericsson.com (version 0.1; February 27, 2018)
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our
% -- journal paper:
% --   Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein,
% --   and Christoph Studer, "VLSI Design of a 3-bit Constant Modulus
% --   Precoder for Massive MU-MIMO", IEEE International Sympsosium on
% --   Circuits and Systems (ISCAS), to appear in 2018
% -- and clearly mention this in your paper
% =========================================================================
function precoder_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runId = 0;       % simulation ID (used to reproduce results)    
    par.U = 16;          % number of single-antenna users
    par.B = 256;         % number of base-station antennas (B>>U)
    par.mod = '16QAM';   % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
    par.trials = 1e3;    % number of Monte-Carlo trials (transmissions)
    par.NTPdB_list = ... % list of normalized transmit power [dB] values
        -10:2:16;        % to be simulated        
    par.precoder = ...   % precoding scheme(s) to be evaluated
        {'ZF','MRT','ZFQ','MRTQ','C2PO','C3PO'};
    par.save = true;     % save results (true,false)
    par.plot = true;     % plot results (true,false)
    
    % *** C2PO and C3PO specific
    %
    % reasonable parameters for C2PO and C3PO with different system confi-
    % gurations
    % please optimize manually for best performance (depends on # of iters)
    %
    % BxU    | mod.  | tau   | delta | rho
    % -------+-------+-------+-------+------
    % 32x16  | BPSK  | 2^-6  | 12.8  | 1.25
    % 64x16  | BPSK  | 2^-7  | 25.6  | 1.25
    % 128x16 | BPSK  | 2^-7  | 25.6  | 1.25
    % 256x16 | BPSK  | 2^-8  | 51.2  | 1.25
    % -------+-------+-------+-------+------
    % 32x16  | QPSK  | 2^-6  | 12.8  | 1.25
    % 64x16  | QPSK  | 2^-7  | 25.6  | 1.25
    % 128x16 | QPSK  | 2^-7  | 25.6  | 1.25
    % 256x16 | QPSK  | 2^-8  | 51.2  | 1.25
    % -------+-------+-------+-------+-------
    % 256x16 | 16QAM | 2^-8  | 51.2  | 1.25
    % -------+-------+-------+-------+-------
    % 256x16 | 64QAM | 2^-8  | 51.2  | 1.25
    
    par.CxPO.tau = 2^(-8); % good for 256x16 with 16-QAM
    par.CxPO.rho = 1.25; % rho = 1/(1-tau*delta) [aka. pushfactor]
    par.CxPO.iterations = 10; % max number of iterations
    
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1};   % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);

% simulation name (used for saving results)
par.simName = ['ERR_',num2str(par.U),'x',num2str(par.B), '_', ...
    par.mod, '_', num2str(par.trials),'Trials'];

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK',
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
end

% compute symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.bps = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% - initialize result arrays (detector x normalized transmit power)
% vector error rate
res.VER = zeros(length(par.precoder),length(par.NTPdB_list));
% symbol error rate
res.SER = zeros(length(par.precoder),length(par.NTPdB_list));
% bit error rate
res.BER = zeros(length(par.precoder),length(par.NTPdB_list));
% error-vector magnitude
res.EVM = zeros(length(par.precoder),length(par.NTPdB_list));
% SINDR
res.SINDR = zeros(length(par.precoder),length(par.NTPdB_list));
% transmit power
res.TxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% receive power
res.RxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% simulation beamforming time
res.TIME = zeros(length(par.precoder),length(par.NTPdB_list));

% compute noise variances to be considered
N0_list = 10.^(-par.NTPdB_list/10);

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.bps,par.trials);

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix and noise vector
    n = sqrt(0.5)*(randn(par.U,1)+1i*randn(par.U,1));
    H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
    
    % algorithm loop
    for d=1:length(par.precoder)
        
        % normalized transmit power loop
        for k=1:length(par.NTPdB_list)
            
            % set noise variance
            N0 = N0_list(k);
            
            % record time used by the beamformer
            starttime = toc;
            
            % beamformers
            switch (par.precoder{d})
                % noise-independent
                case 'ZF',      % ZF beamforming (infinite precision)
                    [x, beta] = ZF(par, s, H);
                case 'ZFQ',     % ZF beamforming (3-bit)
                    [x, beta] = ZF(par, s, H);
                    %Projection to 8PSK
                    xA = abs(real(x)) + 1i*abs(imag(x));
                    x01 = (imag(xA) >= 1/(sqrt(2)-1)*real(xA));
                    x10 = (imag(xA) <= (sqrt(2)-1)*real(xA));
                    x11 = ~x01 & ~x10;
                    xA(x01) = 1i;
                    xA(x10) = 1;
                    xA(x11) = (1+1i)/sqrt(2);
                    x = sign(real(x)).*real(xA) + 1i*sign(imag(x)).*imag(xA);
                    x = x/sqrt(par.B);
                    beta = norm(s,2)^2/(s'*H*x);
                case 'MRT',      % MRT beamforming (infinite precision)
                    [x, beta] = MRT(par, s, H);
                case 'MRTQ',    % MRT beamforming (3-bit)
                    [x, beta] = MRT(par, s, H);
                    %Projection to 8PSK
                    xA = abs(real(x)) + 1i*abs(imag(x));
                    x01 = (imag(xA) >= 1/(sqrt(2)-1)*real(xA));
                    x10 = (imag(xA) <= (sqrt(2)-1)*real(xA));
                    x11 = ~x01 & ~x10;
                    xA(x01) = 1i;
                    xA(x10) = 1;
                    xA(x11) = (1+1i)/sqrt(2);
                    x = sign(real(x)).*real(xA) + 1i*sign(imag(x)).*imag(xA);
                    x = x/sqrt(par.B);                    
                    beta = norm(s,2)^2/(s'*H*x);                
                case 'C2PO',      % C2PO: C1PO with simpler preprocessing
                    [x, beta] = C2PO(par, s, H); 
                case 'C3PO',      % C3PO: C2PO with 3-bits
                    [x, beta] = C3PO(par, s, H);
                otherwise,
                    error('par.precoder not specified')
            end
            
            % record beamforming simulation time
            res.TIME(d,k) = res.TIME(d,k) + (toc-starttime);
            
            % transmit data over noisy channel
            Hx = H*x;
            y = Hx + sqrt(N0)*n;
            
            % extract transmit and receive power
            res.TxPower(d,k) = res.TxPower(d,k) + mean(sum(abs(x).^2));
            res.RxPower(d,k) = res.RxPower(d,k) + mean(sum(abs(Hx).^2))/par.U;
            
            % user terminals can estimate the beamforming factor beta
            shat = beta*y;
            
            % perform user-side detection
            [~,idxhat] = min(abs(shat*ones(1,length(par.symbols)) ...
                -ones(par.U,1)*par.symbols).^2,[],2);
            bithat = par.bits(idxhat,:);
            
            % -- compute error and complexity metrics
            err = (idx~=idxhat);
            res.VER(d,k) = res.VER(d,k) + any(err);
            res.SER(d,k) = res.SER(d,k) + sum(err)/par.U;
            res.BER(d,k) = res.BER(d,k) + ...
                sum(sum(bits(:,:,t)~=bithat))/(par.U*par.bps);
            res.EVM(d,k) = res.EVM(d,k) + 100*norm(shat - s)^2/norm(s)^2;
            res.SINDR(d,k) = res.SINDR(d,k) + norm(s)^2/norm(shat - s)^2;
            
        end % NTP loop
        
    end % algorithm loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',...
            time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.EVM = res.EVM/par.trials;
res.SINDR = res.SINDR/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.TIME = res.TIME/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structures)

if par.save
    save([ par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plots)

if par.plot
    
    % - BER results
    marker_style = {'k-','b:','r--','y-.','g-.','bs--','mv--'};
    figure
    for d=1:length(par.precoder)
        semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    hold off
    grid on
    box on
    xlabel('normalized transmit power [dB]','FontSize',12)
    ylabel('uncoded bit error rate (BER)','FontSize',12);
    if length(par.NTPdB_list) > 1
        axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
    end
    legend(par.precoder,'FontSize',12,'location','southwest')
    set(gca,'FontSize',12);
end

end

%% Zero-forcing beamforming (with infinite precision)
function [x, beta] = ZF(par, s, H)

% normalization constant (average gain)
rho = sqrt((par.B-par.U)/(par.Es*par.U));

% transmitted signal
x = rho*H'/(H*H')*s;

% beamforming factor
beta = 1/rho;

end

%% Maximum ratio transmission (MRT) beamforming (with infinite precision)
function [x, beta, P] = MRT(par, s, H)

% normalization constant
gmrt = 1/sqrt(par.Es*par.U*par.B); % average gain
% gmrt = 1/sqrt(par.Es*trace(H*H')); % instant gain

% precoding matrix
P = gmrt*H';

% transmitted signal
x = P*s;

% scaling factor
beta = sqrt(par.U*par.Es/par.B);

end

%% C2PO: biConvex 1-bit PrecOding with simplified processing
function [x, beta] = C2PO(par,s,H)

% initial guess
x = H'*s;

% preprocessing with approximate inverse
tau = par.CxPO.tau; % step size
Ainvapprox = eye(par.B) - tau*H'*(eye(par.U)-s*s'/norm(s,2)^2)*H ;

% main C1PO algorithm loop
for i=2:par.CxPO.iterations
    x = par.CxPO.rho*(Ainvapprox*x);
    x = min(max(real(x),-1),1) + 1i*min(max(imag(x),-1),1);
end
x = (sign(real(x))+1i*sign(imag(x)))/sqrt(2*par.B);

% scaling factor
beta = norm(s,2)^2/(s'*H*x);

end

%% C3PO
function [x, beta] = C3PO(par,s,H)

% initial guess
x = H'*s;

% preprocessing with approximate inverse
tau = par.CxPO.tau; % step size
Ainvapprox = eye(par.B) - tau*H'*(eye(par.U)-s*s'/norm(s,2)^2)*H ;

% main C1PO algorithm loop
for i=2:par.CxPO.iterations
    x = par.CxPO.rho*(Ainvapprox*x);
    
    %Projection into convex hull of 8PSK
    xA = abs(real(x)) + 1i*abs(imag(x));
    xS = (imag(xA) > (1-sqrt(2))*real(xA)+1) | (imag(xA) > 1/(1-sqrt(2))*(real(xA)-1));
    x01 = xS & (imag(xA) >= 1/(sqrt(2)-1)*real(xA)+1);
    x10 = xS & (imag(xA) <= (sqrt(2)-1)*(real(xA)-1));
    x11 = xS & ((sqrt(2)-1)*(real(xA)+1) <= imag(xA)) & (imag(xA) <= 1/(sqrt(2)-1)*real(xA)-1);
    xP01 = xS & (1/(sqrt(2)-1)*real(xA)-1 <= imag(xA)) & (imag(xA) <= 1/(sqrt(2)-1)*real(xA)+1);
    xP10 = xS & ((sqrt(2)-1)*(real(xA)-1) <= imag(xA)) & (imag(xA) <= (sqrt(2)-1)*(real(xA)+1));
    xA(x01) = 1i;
    xA(x10) = 1;
    xA(x11) = (1+1i)/sqrt(2);
    xA(xP01) = -1/(2*sqrt(2))*(imag(xA(xP01))-1/(sqrt(2)-1)*real(xA(xP01))-1);
    xA(xP01) = xA(xP01) +1i * ((1-sqrt(2))*xA(xP01)+1);
    xA(xP10) = -1/(2*sqrt(2))*(imag(xA(xP10))-(sqrt(2)-1)*real(xA(xP10))+1/(1-sqrt(2)));
    xA(xP10) = xA(xP10) +1i * (1/(1-sqrt(2))*(xA(xP10)-1));
    
    x = sign(real(x)).*real(xA) + 1i*sign(imag(x)).*imag(xA);
end

%Projection into 8PSK
xA = abs(real(x)) + 1i*abs(imag(x));
x01 = (imag(xA) >= 1/(sqrt(2)-1)*real(xA));
x10 = (imag(xA) <= (sqrt(2)-1)*real(xA));
x11 = ~x01 & ~x10;
xA(x01) = 1i;
xA(x10) = 1;
xA(x11) = (1+1i)/sqrt(2);
x = sign(real(x)).*real(xA) + 1i*sign(imag(x)).*imag(xA);
x = x/sqrt(par.B);

% scaling factor
beta = norm(s,2)^2/(s'*H*x);

end
