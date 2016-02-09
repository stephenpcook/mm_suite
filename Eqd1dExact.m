function xend=Eqd1dExact(xin,Monitor)
% Function for a 1-d equidistributed mesh, equidistributing exactly.

% Function takes input of a starting mesh and a Monitor structure
% and outputs a mesh which is equidistributed with MMPDE5, ODE15s for
% up to t_max.
%%%%%%%%%
% INPUT %
%%%%%%%%%
%
% xin     - starting mesh
% Monitor - Structure. Currently accepting Monitor.type as 'points' or
%           'function'. Points provides points and M values then splines
%           the points. Function is defined everywhere. In both cases the 
%           function is differentiated with FD at the following steps.
%
%%%%%%%%%%
% OUTPUT %
%%%%%%%%%%
%
% xout - Final mesh
%
%%%%%%%%%%%%%%
% PARAMATERS %
%%%%%%%%%%%%%%
%
% Stephen Cook, 15-11-2012
if nargin==0
    xin=(0:21)./21;
    Monitor=@(x)exp(-10*x);
end % if nargin
if nargout==0
end

mf= getM(Monitor);
N=length(xin);

h= 1./(N-1);      % Computational grid width
xi= h*((1:N)-1)'; % Computational grid

x0 = xin;   % Initial Physical Grid
xp0= zeros(size(x0));

% MAIN FUNCTIONALITY. WHAT A MESS!
if strcmp(Monitor.type,'points')
    % Should I put a check in to see that the start and end points are the
    % same for x0 and monitor.x? Probably, not now.
    x_m = Monitor.x;
    Mx  = Monitor.M;
    x = equidistribute(x0,x_m,Mx);
elseif strcmp(Monitor.type,'function')
    Mx = Monitor.function(x0);
    x=equidistribute(x0,x0,Mx);
end % if monitor.type==points
xend = x;

function xout = equidistribute(x0,x_m,Mx)
    % PUT THIS INTO SUBFUNCTION. FOR ANALYTIC M CAN ITERATE THIS. (Converge?)
    
    % PART 1
    if 0 % OLD CODE
    %  % Integrate M over x_m by the composite Trapezium Rule
    %  Nm = length(x_m);
    %  IntMi = zeros(1,Nm-1); % Integral of each element
    %  cIM = zeros(1,Nm-1); % Cumilitive Integral of M.
    %  for ii = 1:Nm-1 
    %    IntMi(ii) = 1/2*(Mx(ii)+Mx(ii+1))*(x_m(ii+1)-x_m(ii)); % Trapezium 
    %  end % for ii
    %  cIM=cumsum(IntMi); % Cumilitive integral, from 0 to x_m.
    %  IntM = sum(IntMi); % Integral of whole interval. 
    %  % AGAIN NEED TO BE CLEVERERER

    %  % Place the points of x to satisty equidistribution.
    %  xout = zeros(size(x0));
    %  Nx = length(x0);
    %  x_h = 1/(Nx-1);
    %  % if (x(1) == x_m(1)) good, else do clever, end
    %  xout(1) = x0(1);
    %  Mi=1; % Index of cIM (left point), 1 to Nm-1
    %  for ii = 2:N
    %      Target = (ii-1)*x_h*IntM; % Target integral.
    %      while and(Mi<Nm-1,cIM(Mi)<=Target); % Find the right interval.
    %      				   % Mi is the index of interval and the
    %      				   % left side interval
    %          Mi = Mi + 1;
    %      end % while
    %      xL = x_m(Mi);
    %      xR = x_m(Mi+1);
    %      if Mi==1
    %      ML=0;
    %  else
    %      ML = cIM(Mi-1);
    %  end
    %      M_in = IntMi(Mi);
    %      xout(ii) = xL + (Target-ML)*(xR - xL)/(M_in);
    %      s = (Mx(Mi+1) - Mx(Mi))/(xR-xL);
    %      xout(ii) = xL + sqrt(2*(Target - ML)/s + (Mx(Mi)/s)^2) - Mx(Mi)/s;
    %  end % for ii
    end
    X = x_m;
    U = Mx;
    N = length(x0) - 1;
    
    Y = zeros(size(x0));
    
    intUi = 1/2*(U(1:end-1)+U(2:end)).*(diff(X));
    intU = [0;cumsum(intUi(:))];
    theta = intU(end);
    
    Y(1) = X(1);
    jj = 1;
    for ii = 1:N-1
    Target = ii/N * theta;
    while and(intU(jj) < Target, jj<N)
        jj = jj + 1;
    end
    jj = jj - 1;
    
    XL = X(jj);
    UL = U(jj);
    intUL = intU(jj);
    
    XR = X(jj+1);
    UR = U(jj+1);
    
    target_local = Target - intUL;
    m = (UR-UL)/(XR-XL);
    
    if m==0
      Y(ii+1) = XL + target_local/UL;
    else
      Y(ii+1) = XL + (-UL + sqrt(UL^2 + 2*m*target_local))/m;
    end % if m==0
    
    end % for ii
    Y(N+1) = X(N+1);

    xout = Y;
end

function M=getM(Monitor)
    if strcmp(Monitor.type,'function')
        % Monitor function given
        M = Monitor.function;
    end
    if strcmp(Monitor.type,'points')
        % Pointwise
        % Not yet positivity preserving
        Monitor.pp= spline(Monitor.x, Monitor.M); 
        M= @(x) ppval(Monitor.pp,x);
    end
end % function getM

end % Main function mesh
