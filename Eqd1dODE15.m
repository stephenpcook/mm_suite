function xend=Eqd1dODE15(xin,Monitor,params)
%EQD1DODE15 Construct  non-uniform 1-d mesh using MMPDE5 and ODE15s.
%
% xend = Eqd1dODE15(xin, Monitor, params)
%
% Function takes input of a starting mesh and a Monitor structure
% and outputs a mesh which is equidistributed with MMPDE5, ODE15s for
% up to t_max.
% %%%%%%%%%
% % INPUT %
% %%%%%%%%%
%
% xin     - starting mesh
% Monitor - Structure. Currently accepting Monitor.type as 'points' or
%           'function'. Points provides points and M values then splines
%           the points. Function is defined everywhere. In both cases the
%           function is differentiated with FD at the following steps.
%            .type - 'points' or 'function'
%            .function - @(x)f(x)
%            .x - vector
%            .M - vector
% params  - Structure.
%           .tau
%           .t_max
%           .t_N
%
% %%%%%%%%%%
% % OUTPUT %
% %%%%%%%%%%
%
% xout - Final mesh
%
% %%%%%%%%%%%%%%
% % PARAMATERS %
% %%%%%%%%%%%%%%
%
% tau   - Relaxation parameter in MMPDE5
% t_max - What time final time to be output is
% t_N   - The number of points to evaluate t at. (Need this?)
%
% Contains the subfunctions xt and getM
%
% See also:
% EQD1DEXACT

if nargin==0
    xin=(0:21)./21;
    Monitor=@(x)exp(-10*x);
elseif nargin==2                      % Parameters
    tau= 1;
    t_max=1;
    t_N= 200;
else
    try
        tau = params.tau;
    catch
        tau=1;
    end
    t_max=params.t_max;
    t_N=params.t_N;
end % if nargin
if nargout==0
end

mf= getM(Monitor);
N=length(xin);

h= 1./(N-1);      % Computational grid width
xi= h*((1:N)-1)'; % Computational grid
T= t_max*(0:(t_N-1))./(t_N-1);
T= [0 t_max];

x0 = xin;   % Initial Physical Grid
xp0= zeros(size(x0));

u0=@(x) sin(2*pi*x)+1/2*sin(pi*x);

f= @(t,x) xt(x);
[tout,xout]= ode15s(f,T,x0,xp0);

xend=xout(end,:);

function xtout= xt(x)
%XT Subfunction for constructing approximation to x_t
% x_t = 1/tau*(x_xi*M(x))_xi;

% Forward/Backwards difference estimate of dx/d(xi).
x_xi= (x(2:end)-x(1:end-1))/h;
x_mid= (x(1:end-1)+x(2:end))*0.5;

M= mf(x);

% Linear interpolation for midpoints of M.
M_mid= (M(1:end-1)+M(2:end))*0.5;
% Or why not just extract the values from M?
M_mid= mf(x_mid);

xtout= 1/(tau*h) * [0;...
                (M_mid(2:end).*x_xi(2:end)...
                - M_mid(1:end-1).*x_xi(1:end-1));...
                0];
end % function xt

function M=getM(monitor)
    if strcmp(monitor.type,'function')
        % Monitor function given
        M = monitor.function;
    end
    if strcmp(monitor.type,'points')
        % Pointwise
        % Not yet positivity preserving
        monitor.pp= spline(monitor.x, monitor.M);
        M= @(x) ppval(monitor.pp,x);
    end
end % function getM

end % Main function mesh
