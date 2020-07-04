%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  PBMMI STRING FDTD ASSIGNMENT 
%%%%%  SINGLE STRING FDTD CODE
%%%%% 
%%%%%  FUNCTION FOR FINITE DIFFERENCE TIME DOMAIN SIMULATION OF 
%%%%%  A STIFF STRING INCLUDING LOSS AND BOUNDARY CONDITIONS
%%%%%
%%%%%  TAKES IN A SET OF PHYSICAL STRING AND SIMULATION
%%%%%  PARAMETERS AND RETURNS A MONOPHONIC AUDIO SIGNAL
%%%%%
%%%%%  USES A MATLAB STRUCTURE FOR SIMULATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = A4_String_FDTD_s1760559_Premkumar(opts,phys_param,sim_param)
   % Print options and parameters
   opts;phys_param;sim_param;

   % Copy over parameters, taking into account some options
   T = phys_param.T;          % tension (N)
   r = phys_param.r;          % string radius (m)
   if opts.add_stiffness
      E = phys_param.E;       % Young's modulus (Pa)
   else
      E = 0;
   end
   rho = phys_param.rho;      % density (�?) (kg/m^3) 

   T60 = phys_param.T60;      % T60 (s)
   L = phys_param.L;          % length (m)

   SR = sim_param.SR;         % sample rate (Hz)

   Tf = sim_param.Tf;         % duration of simulation (s)

   xi = sim_param.xi;         % coordinate of excitation (normalised, 0-1)
   famp = sim_param.famp;     % peak amplitude of excitation (N)
   dur = sim_param.dur;       % duration of excitation (s)
   exc_st = sim_param.exc_st; % start time of excitation (s)
   xo = sim_param.xo;         % coordinate of output (normalised, 0-1)

   %%            Derived parameters         %%
   
   A = pi*r^2;                 % string cross-sectional area
   I = 0.25*pi*r^4;            % string moment of inertia
   c = sqrt(T/(rho*A));        % wave speed
   sig = 6*log(10)/T60;        % loss parameter (σ)
   K = sqrt(E*I/(rho*A));      % stiffness constant (κ)
   k = 1/SR;                   % time step

   
   hmin = sqrt(0.5*(c^2*k^2+sqrt(c^4*k^4+16*K^2*k^2)));      % minimal grid spacing for stability
   N = floor(L/hmin);          % number of grid points to update
   h = L/N;                    % actual grid spacing used

   %% Assert Stability Conditions
   
   assert(h>=hmin)             %for stability
   assert(sig>=0)              %for stability

   lambda = c*k/h;             % Courant number (λ)
   mu = k*K/h^2;               % numerical stiffness constant (μ)

   %% I/O 

   Nf = floor(Tf*SR);        % number of time steps for simulation (based on Tf, SR)

   li = floor(xi*N);         % grid index of excitation (based on xi,N,L)
   lo = floor(xo*N);         % grid index of output (based on xo,N,L)

   assert(is_pinteger(li))
   assert(is_pinteger(lo))

   %% Create Force Signal
   
   f = zeros(Nf,1);                   % input force signal
   durint = round(dur*SR);            % duration of force signal, in samples
   exc_st_int = (floor(exc_st*SR))+1; % start time index for excitation

   assert(is_pinteger(durint))
   assert(is_pinteger(exc_st_int))

   for n=exc_st_int:exc_st_int+durint-1
      % For struck implement full hann window 
      if strcmp(opts.input_type,'struck')
         f(n) = famp*0.5*(1-cos(2*pi*(n/durint)));
      % For plucked implement half hann window
      elseif strcmp(opts.input_type,'plucked')
         f(n) = famp*0.5*(1-cos(pi*(n/durint)));
      end
   end

   %% State Variables
   
   u0 = zeros(N,1);           % state at time index n+1
   u1 = zeros(N,1);           % state at time index n
   u2 = zeros(N,1);           % state at time index n-1

   y = zeros(Nf,1);           % output vector

   %% Start and end l-index for for-loop update 
   % l not an element of {0, 1, lstart, N, N + 1} (and considering MATLAB
   % indexing)
   
   lstart = 3;
   lend = N-2;

   % Check if they are integers
   assert(is_pinteger(lstart)) 
   assert(is_pinteger(lend)) 

   %% Main Loop
   
   tic
   for n=1:Nf
      % interior update
      if opts.useforloop
         for l = lstart:lend
            u0(l) = (2*u1(l) + (sig*k-1)*u2(l) + lambda^2*(u1(l-1)-2*u1(l)+u1(l+1)) - mu^2*(u1(l+2)-4*u1(l+1)+6*u1(l)-4*u1(l-1)+u1(l-2)))/(1+sig*k); 
         end
      else
      % Vectorized
         u0(lstart:lend) = (2*u1(lstart:lend) + (sig*k-1)*u2(lstart:lend) + lambda^2*(u1(lstart-1:lend-1)-2*u1(lstart:lend)+u1(lstart+1:lend+1)) - mu^2*(u1(lstart+2:lend+2)-4*u1(lstart+1:lend+1)+6*u1(lstart:lend)-4*u1(lstart-1:lend-1)+u1(lstart-2:lend-2)))/(1+sig*k);
      end
      
      % Boundary updates
      
      if strcmp(opts.bctype,'clamped')
         u0(1) = (2*u1(1) + (sig*k-1)*u2(1) + lambda^2*(-2*u1(1)+u1(2)) - mu^2*(u1(3)-4*u1(2)+6*u1(1)))/(1+sig*k);
         u0(2) = (2*u1(2) + (sig*k-1)*u2(2) + lambda^2*(u1(1)-2*u1(2)+u1(3)) - mu^2*(u1(4)-4*u1(3)+6*u1(2)-4*u1(1)))/(1+sig*k);
         u0(N-1) = (2*u1(N-1) + (sig*k-1)*u2(N-1) + lambda^2*(u1(N-2)-2*u1(N-1)+u1(N)) - mu^2*(u1(N-3)-4*u1(N-2)+6*u1(N-1)-4*u1(N)))/(1+sig*k);
         u0(N) = (2*u1(N) + (sig*k-1)*u2(N) + lambda^2*(-2*u1(N)+u1(N-1)) - mu^2*(u1(N-2)-4*u1(N-1)+6*u1(N)))/(1+sig*k);
      elseif strcmp(opts.bctype,'simply_supported')
         u0(1) = (2*u1(1) + (sig*k-1)*u2(1) + lambda^2*(-2*u1(1)+u1(2)) - mu^2*(u1(3)-4*u1(2)+5*u1(1)))/(1+sig*k);
         u0(2) = (2*u1(2) + (sig*k-1)*u2(2) + lambda^2*(u1(1)-2*u1(2)+u1(3)) - mu^2*(u1(4)-4*u1(3)+6*u1(2)-4*u1(1)))/(1+sig*k);
         u0(N-1) = (2*u1(N-1) + (sig*k-1)*u2(N-1) + lambda^2*(u1(N-2)-2*u1(N-1)+u1(N)) - mu^2*(u1(N-3)-4*u1(N-2)+6*u1(N-1)-4*u1(N)))/(1+sig*k);
         u0(N) = (2*u1(N) + (sig*k-1)*u2(N) + lambda^2*(-2*u1(N)+u1(N-1)) - mu^2*(u1(N-2)-4*u1(N-1)+5*u1(N)))/(1+sig*k);
      end

      % send in input (with force coefficient)
      u0(li) = u0(li) + (k^2/(rho*A*h))*f(n);

      % read output
      if strcmp(opts.output_type,'displacement')
         y(n) = u0(lo);
      elseif strcmp(opts.output_type,'velocity')
      % ''velocity'' read-out here
         y(n) = (u0(lo) - u1(lo))*SR;
      end

      %% Plotting
      
      if (opts.plot_on)
         % draw current state
         if n==1
            figure;
            h1=plot([1:N]'*h, u0, 'k');
            xlim([0 L]);
            ylim([-0.005 0.005]); 
            xlabel('position (m)');
         else
            set(h1,'ydata',u0);
            drawnow;
         end
         fprintf('n=%d out of %d\n',n,Nf);
      end

      % Shift states to step forward in time
      u2 = u1;
      u1 = u0;
   end

   %read last samples of output
   for n=Nf-4:Nf
      fprintf('y(%d) = %.15g\n',n,y(n));
   end
   toc

   %%%%% plot spectrum
   if (opts.plot_on)
      figure
      yfft = 10*log10(abs(fft(y)));
      plot([0:Nf-1]'/Nf*SR, yfft, 'k')
      xlabel('freq (Hz)')
   end
end

%is positive integer?
function y=is_pinteger(x)
   y=((mod(x,1)==0) && (x>0));
end
