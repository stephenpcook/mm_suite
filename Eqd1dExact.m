function xend = Eqd1dExact(xin,Monitor)
%EQD1DEXACT Equidistribute to a monitor function structure.
%
% xend = Eqd1dExact(xin, Monitor)
%
% Function takes input of a starting mesh and a Monitor structure
% and outputs a mesh which is equidistributed with MMPDE5, ODE15s for
% up to t_max.
%
% %%%%%%%%%
% % INPUT %
% %%%%%%%%%
%
% xin     - starting mesh
% Monitor - Structure. Currently accepting Monitor.type as 'points' or
%           'function'. Points provides points and M values then splines
%           the points. Function is defined everywhere. In both cases the
%           function is differentiated with FD at the following steps.
%
% %%%%%%%%%%
% % OUTPUT %
% %%%%%%%%%%
%
% xend - Final mesh
%
% %%%%%%%%%%%%%%
% % PARAMATERS %
% %%%%%%%%%%%%%%
%
% Stephen Cook, 15-11-2012
%
% See also:
% EQD1DODE15

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

x0 = xin(:);   % Initial Physical Grid
xp0= zeros(size(x0));

% MAIN FUNCTIONALITY.
if strcmp(Monitor.type,'points')
    % Should I put a check in to see that the start and end points are the
    % same for x0 and monitor.x? Probably, not now.
    x_m = Monitor.x(:);
    Mx  = Monitor.M(:);
    x = equidistribute(x0,x_m,Mx);
elseif strcmp(Monitor.type,'function')
    Mx = Monitor.function(x0);
    x = equidistribute(x0,x0,Mx);
end % if monitor.type==points

xend = reshape(x,size(xin));

function xout = equidistribute(x0,x_m,Mx)
    % Find points xout that equidistribute the lin. interp. of (x_m,Mx)
    %
    %  IN
    % x0  - the output will be the same size as x0 (length N)
    % x_m - the knots of the Monitor function (X)
    % Mx  - the values of the monitor function to be equidistributed to (U)
    %  OUT
    % xout - (Y)


    X = x_m;
    U = Mx;
    N = length(x0) - 1;
    
    Y = zeros(size(x0));
    
    %   $$intUi(i) = \int_{X(i)}^{X(i+1)} u(x) dx ,$$
    %   $$intU(i) = \int_{X(1)}^{X(i)} u(x) dx$$
    % and
    %   $$theta = \int_{X(1)}^{X(N)} u(x) dx$$
    % Assuming U is piecewise linear.
    intUi = 1/2*(U(1:end-1)+U(2:end)).*(diff(X));
    intU = [0;cumsum(intUi(:))];
    theta = intU(end);
    
    Y(1) = X(1);
    jj = 1;
    for ii = 1:N-1
    Target = ii/N * theta;
    % Find the interval in which Y(ii) lies.
    while and(intU(jj) < Target, jj<N)
        jj = jj + 1;
    end % while
    jj = jj - 1;
    
    XL = X(jj);
    UL = U(jj);
    intUL = intU(jj);
    
    XR = X(jj+1);
    UR = U(jj+1);
    
    target_local = Target - intUL;
    m = (UR-UL)/(XR-XL);
    
    % In this interval
    %   $$U(x) = m*(x-XL) + UL .$$
    % We want to find Y such that
    %   $$\int_XL^Y U(x) dx = target_local , $$
    % so we integrate to get
    %   $$ Y = XL + (-UL + sqrt(UL^2 + 2*m*target_local))/m . $$
    % For small m, this gives huge error, so we rearrange as
    Y(ii+1) = XL + 2*target_local/(UL + sqrt(UL^2 + 2*m*target_local));
    
    end % for ii
    Y(N+1) = X(N+1);

    xout = Y;
end

function M = getM(Monitor)
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
