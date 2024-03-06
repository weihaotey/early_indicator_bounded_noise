function y = saddle_discrete(x,r,epsilon,varargin)
    p = inputParser;
    addParameter(p, 'dist' ,'Normal',@isstring);
    parse(p,varargin{:});
    if p.Results.dist == 'Normal'
        pd = makedist( 'Normal','sigma',0.5*epsilon);
        t = truncate(pd,-epsilon,epsilon);
        y = x + r - x.^2 + random(t,size(x));
    elseif p.Results.dist == 'Uniform'
        y = x + r - x.^2 + unifrnd(-epsilon,epsilon,size(x));
    elseif p.Results.dist == 'Unboundednormal'
        y = x + r - x.^2 + normrnd(0,epsilon/2,size(x));
    elseif p.Results.dist == 'Deterministic'
        y = x + r - x.^2;
    elseif p.Results.dist == '1duni'
        f = @(x) 3*((1-exp(-x))./(1+exp(-x))); % an invertible map
        y = f(x) - r + unifrnd(-epsilon,epsilon,size(x));
    elseif p.Results.dist == '1dnorm'
        f = @(x) 3*((1-exp(-x))./(1+exp(-x)));
        pd = makedist( 'Normal','sigma',0.5*epsilon);
        t = truncate(pd,-epsilon,epsilon);
        y = f(x) - r + random(t,size(x));
    elseif p.Results.dist == '1dnormunbounded'
        f = @(x) 3*((1-exp(-x))./(1+exp(-x)));
        y = f(x) - r + normrnd(0,epsilon/2,size(x));
    elseif p.Results.dist == '1dlinear'
        epsilon = epsilon*(1-r); % adjust the noise so that boundary is fixed
        y = r*x + unifrnd(-epsilon,epsilon,size(x));
    elseif p.Results.dist == '1dlinearnormal'
        epsilon = epsilon*(1-r); % adjust the noise so that boundary is fixed
        pd = makedist( 'Normal','sigma',0.5*epsilon);
        t = truncate(pd,-epsilon,epsilon);
        y = r*x + random(t,size(x));
    elseif p.Results.dist == '1dlinearunbounded'
        epsilon = epsilon*(1-r); % adjust the noise so that boundary is fixed
        y = r*x + normrnd(0,epsilon/2,size(x));
    elseif p.Results.dist == '1dmodifieduni' % modified for early warning signal
        % f = @(x) (3+4*r)*((1-exp(-x-r))./(1+exp(-x-r)))-4*r; % an invertible map
        f = @(x) erf((2+r)*x-r^2);
        y = f(x) + unifrnd(-epsilon,epsilon,size(x));
    end
end