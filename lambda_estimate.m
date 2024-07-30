% s = settings;
% s.matlab.editor.AllowFigureAnimation.TemporaryValue = 1;
tic; close all;
paramchange = 0; % changing paramter slowly across time
timeserieslength = 1e6+1; % total length of the time series to be produced
nsample = 500; % number of experiment per parameter
nsample1 = nsample; % for single trajectory
nsample2 = nsample1;
sigma = 0.1; % bound of noise or 2*std
h = 0.01; % time step for continuous time
n_avg = 0; % the number of iteration to ensure convergence
transient = 100; % initial amount of iteration to get rid of transient
nn = 101; % number of parameters to run
if paramchange==1
    nn = timeserieslength;
end
sample_size = 100; % try different sample size (or repeat experiment for single trajectory)
n = sample_size;

d2 = 0; % if use two dimensional map
plot_on = 0; % whether to plot the histogram
mapType = "1duni"; % type of deterministic map "1duni"
linear_map = 0; % deterministic map is linear
negativestart = 0; % whether to start at negative initial point if less than 0
discrete = 1; % whether the system used is discrete
average = 0; % whether to use lots of points, instead of one traj
rup = 0.3; rdown = -0.5;
rr = linspace(rdown,rup,nn);
if d2==1
    s = ones(2,nn);
else
    s = ones(1,nn);
end
x_all = {}; % to record all x values
x_all_all = {}; % to record different x
limitx = ones(2,nn); % record the limit values of data
rparam = linspace(rdown,rup,timeserieslength); % paramter to run in changing param case

continuous_change = 1; % whether to use initial data continuously from previous parameter
forcenegative = 0; % whether to change initial condition to -1 after 0.32
j = 1; % for different sample size
% for n = sample_size
for nsample = nsample1
    if paramchange == 1
        % non autonomous case where parameter changes slowly 
        x0 = 1; % initial point
        for k = 1:nn
            r = rparam(k);
            if k==1
                x = ones(1,n)*x0; % start all points from x0
                for i = 1:transient % run the system a few round to ensure it is in minimal invariant set
                    x = saddle_discrete(x,r,sigma,"dist",mapType);
                end
            end
            x = saddle_discrete(x,r,sigma,"dist",mapType);
            x_all{k}{j} = x; % record the data
            if mod(k,nsample)==0
            % if k >= nsample
                disp(k)
                xx = ones(1,nsample);
                for ii = 1:nsample
                    xx(ii) = x_all{k-ii+1}{j}(1);
                end
                if plot_on
                    histogram(xx)
                    title(num2str(r))
                    drawnow
                end
                if d2 == 0
                    s(k) = std(xx);
                end
            end
        end
        
    else
        % for simulate data for each parameter
        for k = 1:nn
            r = rr(k); % parameter
            if d2 == 0
                if mapType{1}(1)=='1'
                    x0 = 1;
                else
                    x0 = 1;
                end
                if continuous_change == 0 || k == 1 % restart initial point i   f not continuous change
                    x = ones(1,n)*x0; % start all points from x0
                    for i = 1:transient % run the system a few round to ensure it is in minimal invariant set
                        x = saddle_discrete(x,r,sigma,"dist",mapType);
                    end
                end
                for i = 1:n_avg+nsample
                    x = saddle_discrete(x,r,sigma,"dist",mapType);
                    if i>n_avg % take the next few iterations as different samples
                        x_all{k}{i-n_avg}{j} = x;
                        if any(x<0) && negativestart && mapType{1}(1)=='1' && not(continuous_change)
                            break
                        end
                    end
                end
                if any(x<0) && negativestart && mapType{1}(1)=='1' && not(continuous_change) || (r >= 0.32 && forcenegative)  % use the initial value -3 if any data goes below 0 
                    x0 = -1;
                    x = ones(1,n)*x0; % start all points from x0
                    for i = 1:transient
                        x = saddle_discrete(x,r,sigma,"dist",mapType);
                    end
                    for i = 1:n_avg+nsample
                        x = saddle_discrete(x,r,sigma,'dist',mapType);
                        if i>n_avg
                            x_all{k}{i-n_avg}{j} = x;
                        end
                    end
                end
            else
                x0 = 1; 
                x = ones(n,2)*x0; % start all points from x0
                for i = 1:n_avg+nsample
                    rrrr = unifrnd(0,1,[n,1]);
                    theta = unifrnd(-pi,pi,[n,1]);
                    % randnum = normrnd(0,sigma,n,2);
                    randnum = sigma*rrrr.*[cos(theta),sin(theta)];
                    % x = henon(x,r,0.3,2,0);
                    x = henon(x,r,0.3,2,0)+randnum;
                    if i>n_avg % take the next few iterations as different samples
                        x_all{k}{i-n_avg}{j} = x;
                    end
                end
            end
            if plot_on
                if average == 1
                    histogram(x)
                else
                    xx = ones(1,nsample);
                    for ii = 1:nsample
                        xx(ii) = x_all{k}{ii}{j}(1);
                    end
                    histogram(xx)
                end
                title(num2str(r))
                drawnow
            end
            if d2 == 0
                if average == 0
                    xx = ones(1,nsample);
                    for ii = 1:nsample
                        xx(ii) = x_all{k}{ii}{j}(1);
                    end
                end
                s(k) = std(xx);
            end
            disp(k)
            % if plot_on
            %     if d2 == 0
            %         limitx(:,k) = [min(x),max(x)];
            %         hh = histogram(x,20);
            %     else
            %         s(1,k) = std(x(:,1)); s(2,k) = std(x(:,2));
            %         limitx(:,k) = [min(x(:,1)),max(x(:,1))];
            %         hh = histogram(x(:,1),20);
            %     end
            %     title(num2str(r))
            %     % axis([-1,2,0,1000])
            %     drawnow
            % end
        end
    end
    disp(j)
    j = j + 1;
end
toc;

%% CURRENT LATEST
% find the optimal tail distribution from the stationary distribution
% using least square model
tic;
clear x0 
% close all
writevideo =0;
if writevideo
    vidfile1 = VideoWriter('exponential_loglog','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end
linear_map = 0; % change whether linear or nonlinear
boun_compare = 0; % whether to test different boundary points
sample_size = 10; % change the sample size
edge_bin = 1; % use minmax data as edge bin
plot_on = 0; % whether to plot the figures
plot_on_large = 0; % whether to plot the large step
display_on = 0; % whether to display fval and others
boundary_cheat = 0; % whether to have known boundary
lambda_all_all = {}; % record all lambda
lambda_each = {}; % record lambda for nsample jj
hyperlen = 31; % number of hyperparameters to run
norm_all = zeros(sample_size,hyperlen,3); % record the norm of real lambda and estimated
n_n = nn;
weight_on = 0; weightpower = 1; % whether to introduce weight to the least square objective function
brute = 0; brutequad = 1; brutequad11 = 0; brutequad2 = 0; brutequadnolog=0; linearmethod = 0; extrapolate = 0; % which method to use
unknowns3 = 0; % whether to use 3 unknowns for quad method
extol_maxh = 1; % use 0.01 portion of max h value
minheight = 2; % minimum height of histogram to include if extol_max is false
bounleft = 0; % whether to use x-(bounmin) of used data
fval_tol1 = 3e-1; fval_tol = fval_tol1; % tolerance of function evaluation
left_data = 1; % whether to test left or right data
disp_on = 0; % whether to display result
testboun = 0; % whether to test some value away from boundary
boundata = 0.12; % only take data this value away from edge (if testboun==1)
smooth_on = 0; % whether to smoothen the histrogram data
boun_data = 0; % whether to use the limit from data
uneven = 0; % whether to use uneven bins
if uneven == 1
    k = (atan(linspace(-1.3,1.3,BinNo+1))/atan(1.3)+1)/2; % uneven histogram
end
calc_real_lambda = 1; % whether to calculate real lambda

% initialisation
if boun_compare == 1
    boun_n = 30;
    bounmin_all = ones(boun_n,n_n,sample_size); % record data boundary
    lambda_boun_all = ones(boun_n,n_n,sample_size);
    boun_vec_all = ones(boun_n,n_n);
else
    bounmin_all = ones(1,n_n); % record data boundary
end
lambda_all = ones(1,n_n); % record lambda estimate
real_boun_all = ones(1,n_n); % record real boundary point
real_lambda_all = ones(1,n_n); % record real lambda
target_lambda_all = ones(1,n_n); % record targeted lambda given data boundary
param = ones(1,n_n); % record parameter
A_all = ones(2,n_n); % record brute param estimate
if unknowns3 == 1
    A_all = ones(3,n_n);
end
fval_all = ones(1,n_n); % record error
if linearmethod || extrapolate || brute
    Aloglog = ones(2,n_n); % record resulting param estimate
else
    Aloglog = ones(3,n_n);
end
hist_all = {}; % record all the histogram data
% for jjj = 1:length(sample_size)
jjjj = 1; % initialise for hyperparameter iteration

% hyperparameters
BinNo = 300; % the number of bin in histogram
data_portion = 0.3; % portion of data point to fit

% for BinNo = 20

for data_portion = [0.1,0.35,0.6]
pp = round(linspace(5,100,hyperlen));

% pp = linspace(0.01,0.8,hyperlen); % run through different data portion
% for BinNo = [10,20,30]
jjj = 1;
for BinNo = pp 
% for data_portion = pp
% for nsample1 = nsample2
% for data_portion = 0.6
    % BinNo = ceil(nsample1/500);
    lambda_record = zeros(1,nn); % for adding sample
    if average ~= 1
        single_traj = nsample1;
        nsample = sample_size;
    end
    % for j = [18,41,70,89] % [1,46,71,79]
    %     for boun_count = 1:boun_n
            if display_on
            disp(j)
            end
        for jj = 1:nsample %1:nsample
        for j = 1:length(rr) %1:nn
        
            
            
            fval_tol = fval_tol1;
            a_para = rr(j); % retrieve a parameter value 
            param(j) = a_para; % record parameter
            % epsilon = sigma; % retrieve epsilon parameter value  
            % c = cr_all(1,j); r = cr_all(2,j); % retrieve center and radius of the box
            % retrieve the time-series data
            if paramchange == 1
                if mod(j,single_traj)~=0
                    continue
                end
            end
            
            if average == 1
                xx = x_all{j}{jj}{1}; 
            else
                xx = zeros(1,single_traj); % compile each trajectories
                for ii = 1:single_traj
                    if paramchange == 1
                        xx(ii) = x_all{j-single_traj+ii}{1}(jj);
                    else
                        xx(ii) = x_all{j}{ii}{1}(jj);
                        if length(nsample2) > 1
                            xx(ii) = x_all{j}{ii}{jjj}(jj);
                        end
                    end
                end
            end

            if d2==1 % take the first coordinate for 2-D data
                xx = xx(:,1);
            end
            if mapType{1}(1)=='1'
                f = @(x) 3*((1-exp(-x))./(1+exp(-x))) - a_para; % the mapping
                % f = @(x) (3+4*a_para)*((1-exp(-x-a_para))./(1+exp(-x-a_para)))-4*a_para;
            elseif mapType{1}(1)=='x'
                f = @(x) 3*((1-exp(-x))./(1+exp(-x))) - a_para - sigma*abs(x); % the mapping
            else
                f = @(x) x+a_para-x^2;
            end
            if plot_on
                figure(2)
            end
            if smooth_on
                [h,x]=ksdensity(xx,'NumPoints',BinNo); width = x(2)-x(1); % smoothen the stationary dist curve
            else
                if uneven
                    hh = histogram(xx,min(xx)+(max(xx)-min(xx))*kk); 
                else
                    if edge_bin
                        [h,x] = histcounts(xx,linspace(min(xx),max(xx),BinNo+1)); width = x(2)-x(1);
                    else
                        [h,x] = histcounts(xx,BinNo); width = x(2)-x(1);
                    end
                end
                % title(num2str(a_para)); 
                % x = hh.BinEdges; h = hh.Values; width = hh.BinWidth; % take data of the distribution
                if all(width == 'nonuniform'), width = x(2)-x(1);end
                x = (x(2:end)+x(1:end-1))/2; % take center of the bins
            end
            
            xx = x;

            % only take data near left boundary
            if extol_maxh == 1
                extol = max(h/norm(h,1)/width)*1e-2; % exclude data where h value is small
            else
                extol = minheight/norm(h,1)/width; % exclude data where h value is small
            end

            if uneven == 1
                h = h/sum(h.*(kk(2:end)-kk(1:end-1)));
            else
                h = h/norm(h,1)/width; % so that area is 1
            end
            hist_all{j} = [x;h]; % record the data
            % h = data_all{j}(2,:); % retrieve bins
            % width = r*2;
            % h = h/norm(h,1)/width; % so that area is 1
            
            if boun_data % for left tail distribution
                if left_data == 1
                    bounmin = min(xx); % the lower boundary from data
                else
                    bounmin = max(xx); % the upper boundary from data   
                end
            else
                if left_data == 1
                    bounmin = min(x)-width; % the lower boundary from hist
                else
                    bounmin = max(x)+width; % the upper boundary from hist
                end
            end 
            % bounmin = c-r; % the lower boundary
            % bounmax = c+r; % upper boundary
            syms yy; fprime = matlabFunction(diff(f(yy))); % derivative of the map
        
            
            % take data within data_portion, h value sufficiently large
            % (avoid infinite when log) and not too far from edge
            if testboun
                if left_data
                    w = find(x<(bounmin+boundata)& abs(h)>extol); % only test up to epsilon away from boundary
                else
                    w = find(x>(bounmin-boundata)& abs(h)>extol); % only test up to epsilon away from boundary
                end
            else
                if uneven
                    if left_data
                        w = find(cumsum(h.*(kk(2:end)-kk(1:end-1)))<data_portion & abs(h)>extol ); % only test some data portion & log(xx-bounmin)<0
                    else
                        w = find(cumsum(h.*(kk(2:end)-kk(1:end-1)))>1-data_portion & abs(h)>extol) ; % only test some data portion& log(abs(xx-bounmin))<0
                    end
                else
                    if left_data
                        w = find(cumsum(h)*width<data_portion & abs(h)>extol ); % only test some data portion & log(xx-bounmin)<0
                    else
                        w = find(cumsum(h)*width>1-data_portion & abs(h)>extol); % only test some data portion & log(abs(xx-bounmin))<0
                    end
                end
            end
            
            if left_data == 1
                w = w(2:end); % get rid of end point 
            else
                w = w(1:end-1); % get rid of end point 
            end
            x1 = xx(w); % very tail of the distribtuion to test
            if bounleft == 1
                bounmin = x1(1)-width;
            end
            
            
            % take real lambda
            % c = cr_all(1,j); r = cr_all(2,j); % retrieve center and radius of the box
            % bounm = c-r; % the lower boundary
            if calc_real_lambda && abs(bounmin)<10
                bounm = fzero(@(x)(f(x)-sigma-x),bounmin); % find the actual boundary
                % find the boundary around the estimates
                if linear_map
                    bounm = -sigma;
                else
                if abs(bounm-bounmin)>1
                    for bountry = linspace(bounmin-1,bounmin+1,1000)
                        bounm = fzero(@(x)(f(x)-sigma-x),bountry);
                        if abs(bounm-bounmin)<1
                            break
                        end
                    end
                end
                end

                % bounm = fzero(@(x)(f(x)-sigma-x),1);
                if boundary_cheat
                    if linear_map == 1
                        bounmin = -sigma; % for linear map fixed boundary simulation
                    else
                        bounmin = bounm; % for non-linear map boundary
                    end
                end
                real_lambda = fprime(bounm); % real lambda
                target_lambda = fprime(bounmin); % for 1-exp(-x)/(1+exp(-x)) S map
                if linear_map == 1
                    real_lambda = a_para; % for linear map
                    target_lambda = 1-(sigma*(1-a_para))/abs(bounmin); % target lambda from boundary observed 
                end
                real_boun_all(1,j) = bounm;
                real_lambda_all(1,j) = real_lambda;
                target_lambda_all(1,j) = target_lambda;
                if boun_compare == 1
                    if jj == 1
                    boun_vec = linspace(bounm,bounmin,boun_n); 
                    end
                    boun_plus = boun_vec(boun_count); 
                    bounmin = boun_plus;
                end
            end
            boun = bounmin; % boundary of the minimal invariant set

            opts=optimset('Display','off');
            if isempty(x1),continue;end
            % figure(1)
            % brute or log-log linear regression
            % hyper case
            if brute 
                % hyper case
                opts=optimset('Display','off');
                if plot_on
                figure(2)
                plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
                end
                repeat = 0; % to know whether need to recalculate
                fval1 = 10; % initialise large fval
                count = 1;
                while fval1 > fval_tol % repeat if tolerance not met
                    A1 = optimvar('A1',2);
                    fun = log(abs(x1-bounmin));
                    fun1 = A1(1)*exp(A1(2)*fun.*fun); % hyper tail
                    weight = 1./linspace(1,length(h(w)),length(h(w))).^weightpower ; % to apply weightage for points near boundary is priotised
                    if weight_on
                        obj1 = sum(weight.*(fun1-h(w)).^2); % objective function
                    else
                        obj1 = sum((fun1-h(w)).^2); % objective function
                    end
                    lsqproblem1 = optimproblem("Objective",obj1); % optimization problem
                    if count==1
                        x0.A1 = [1,-1/10]; % initial point of a
                    else
                        x0.A1 = sol1.A1';
                    end
                    if repeat  == 1 % repeat the calculation with random initial
                        x0.A1 = [unifrnd(0,1),unifrnd(-10,0)];
                    end
                    % show(lsqproblem)
                    [sol1,fval1] = solve(lsqproblem1,x0,'Options',opts);
                    repeat = 1;
                    count = count +1;
                    if mod(count,10) == 0
                        % disp('doubled fval tol');
                        fval_tol = fval_tol*2; 
                    end
                end
                A_all(:,j) = sol1.A1; % record fit parameter
                fitted_hyper = evaluate(fun1,sol1); % evaluate function at solution
                lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                if plot_on
                    hhh = plot(x1,fitted_hyper,"Color",'red'); % plot hyper best fit
                    title(strcat('a = ',num2str(a_para),'$\lambda$ = ',num2str(lambda)),'interpreter','latex')
                    % xlim([c-r,max(x1)])
                end
                hold off
                if display_on
                    disp(fval1)
                end
            elseif brutequad
                opts=optimset('Display','off');
                % plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
                repeat = 0; % to know whether need to recalculate
                for brutequad1 = brutequad11
                fval1 = 10; % initialise large fval
                count = 1;
                while fval1 > fval_tol % repeat if tolerance not met
                    if unknowns3 == 1
                        A1 = optimvar('A1',3);
                    else
                        A1 = optimvar('A1',2);
                    end
                    A1(2).UpperBound = 0;
                    fun = log(abs(x1-boun));
                    if unknowns3 == 1
                        fun1 = A1(2)*(fun-A1(1)).*(fun-A1(1)) + A1(3); % hyper tail
                    else
                        if brutequad1 == 1
                            fun1 = A1(1) + A1(2)*fun.*fun; % hyper tail;
                        elseif brutequad2 == 1
                            fun1 = A1(2)*(fun.*fun-2*fun.*log(abs(fun)))+A1(1)*fun; % improved asymptotic
                        else
                            fun1 = A1(1)*fun + A1(2)*fun.*fun; % hyper tail;
                        end
                    end
                    weight = 1./linspace(1,length(h(w)),length(h(w))).^weightpower ; % to apply weightage for points near boundary is priotised
                    if weight_on
                        obj1 = sum(weight.*(fun1-log(h(w))).^2)/sum(weight); % objective function
                    else
                        obj1 = sum((fun1-log(h(w))).^2); % objective function
                    end
                    
                    lsqproblem1 = optimproblem("Objective",obj1); % optimization problem
                    % constraint = A1<=0;
                    % lsqproblem1 = optimproblem("Objective",obj1,"Constraints",constraint); % optimization problem
                    if count==1
                        if unknowns3 == 1
                            x0.A1 = [1,-1/10,0]; % initial point of a
                        else
                            x0.A1 = [1,-1/10]; % initial point of a
                        end
                    else
                        x0.A1 = sol1.A1';
                    end
                    if repeat  == 1 % repeat the calculation with random initial
                        if unknowns3 == 1
                            x0.A1 = [unifrnd(0,10),unifrnd(-10,0),unifrnd(-10,10)];
                        else
                            x0.A1 = [unifrnd(0,10),unifrnd(-10,0)];
                        end
                    end
                    % show(lsqproblem)
                    [sol1,fval1] = solve(lsqproblem1,x0,'Options',opts);
                    repeat = 1;
                    count = count +1;
                    if mod(count,10) == 0
                        if display_on 
                            disp('doubled fval tol'); 
                        end 
                        fval_tol = fval_tol*2; 
                    end
                end
                A_all(:,j) = sol1.A1; % record fit parameter
                logx = log(abs(x(w(1):end)-boun));
                p = log(abs(x(w)-boun));
                if plot_on
                    figure(2)
                    hhh = plot(p,log(h(w)),'.-'); hold on;
                    plot(logx,log(h(w(1):end)))
                    plot(p,evaluate(fun1,sol1),'.-');
                    title(strcat('a = ',num2str(a_para),', length of edge data: ',num2str(length(w))))
                    % xlim([c-r,max(x1)])
                    hold off
                end
                if display_on
                    disp(fval1)
                end
                if length(brutequad11) == 1
                    lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                elseif brutequad1 == 0
                    fval0 = fval1;
                elseif (abs(fval1) < abs(fval0) && brutequad1 == 1) || brutequad1 == 0
                    lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                end
                
                end
            elseif brutequadnolog
                opts=optimset('Display','off');
                % plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
                repeat = 0; % to know whether need to recalculate
                fval1 = 10; % initialise large fval
                count = 1;
                while fval1 > fval_tol % repeat if tolerance not met
                    A1 = optimvar('A1',2);
                    A1(2).UpperBound = 0;
                    fun = log(abs(x1-boun));
                    fun1 = exp(A1(1)*fun + A1(2)*fun.*fun); % hyper tail
                    weight = 1./linspace(1,length(h(w)),length(h(w))).^weightpower ; % to apply weightage for points near boundary is priotised
                    if weight_on
                        obj1 = sum(weight.*(fun1-h(w)).^2)/sum(weight); % objective function
                    else
                        obj1 = sum((fun1-h(w)).^2); % objective function
                    end
                    constraint = A1(2)<=0;
                    % lsqproblem1 = optimproblem("Objective",obj1); % optimization problem
                    lsqproblem1 = optimproblem("Objective",obj1,"Constraints",constraint); % optimization problem
                    if count==1
                        x0.A1 = [1,-0.1]; % initial point of a
                    else
                        x0.A1 = sol1.A1';
                    end
                    if repeat  == 1 % repeat the calculation with random initial
                        x0.A1 = [unifrnd(0,10),unifrnd(-10,0)];
                    end
                    % show(lsqproblem)
                    [sol1,fval1] = solve(lsqproblem1,x0,'Options',opts);
                    repeat = 1;
                    count = count +1;
                    if mod(count,10) == 0
                        if display_on 
                            disp('doubled fval tol'); 
                        end 
                        fval_tol = fval_tol*2; 
                    end
                end
                A_all(:,j) = sol1.A1; % record fit parameter
                logx = log(abs(x(w(1):end)-boun));
                p = log(abs(x(w)-boun));
                if plot_on
                    figure(2)
                    hhh = plot(x(w),h(w),'.-'); hold on;
                    plot(x,h)
                    plot(x(w),evaluate(fun1,sol1),'.-');
                    title(strcat('a = ',num2str(a_para),', length of edge data: ',num2str(length(w))))
                    % xlim([c-r,max(x1)])
                    % disp(fval1)
                end
                lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                hold off
            elseif linearmethod % approximate linear equation
                A = optimvar('A',2);
                fun = A(1)*log(-log(abs(xx(w)-boun))) + A(2); % linear equation
                obj1 = sum((fun-log(-log(h(w)))).^2); % objective function
                lsqproblem0 = optimproblem("Objective",obj1); % optimization problem
                x0.A = [2,-0.8]; % initial point of a
                % show(lsqproblem)
                [sol0,fval0] = solve(lsqproblem0,x0,'Options',opts);
                Aloglog(:,j) = sol0.A; % record fit parameter
                % plot the result
                loglogx = log(-log(xx(2)-boun)); % the value of intersection chosen for log-log|x-x0|
                p = log(-log(abs(x(w)-boun)));
                if plot_on
                    plot(p,log(-log(h(w))),'.-'); hold on;
                    plot(p,evaluate(fun,sol0),'.-');
                    plot(loglogx,sol0.A(1)*loglogx+sol0.A(2),'.');
                    plot(p,2*p+log(-1/(2*log(real_lambda))),'.-');
                    plot(loglogx,2*loglogx+log(-1/(2*log(real_lambda))),'.-');
                    hold off
                    legend('data','best fit','actual line','Location','best')
                end
                aa = (sol0.A(1)-2)*loglogx + sol0.A(2); % calculate the y-intercept for y=2x+c
            elseif extrapolate % technically this method only use last two data to extrapolate
                findind = find(abs(h) > 1e-6);
                loglogx = log(-log(xx(2)-boun)); % the value of intersection chosen for log-log|x-x0|
                p = log(-log(abs(x(w)-boun))); q = log(-log(h(w))); % data points for log-log scale
                plotx = [loglogx p];
                pq = interp1(p,q,plotx,'linear','extrap');
                if plot_on
                    plot(p,log(-log(h(w))),'.-'); hold on;
                    plot(plotx,pq);
                    plot(plotx,2*plotx+log(-1/(2*log(real_lambda))),'.-');
                    hold off
                    legend('data','best fit','actual line','Location','best')
                end
                aa = interp1(p,q,loglogx,'linear','extrap')-2*loglogx;
            else % exponential term in the log log scale
                repeat = 0; % to know whether need to recalculate
                fval0 = 1; % initialise large fval
                while fval0 > fval_tol % repeat if tolerance not met
                    A = optimvar('A',3);
                    pp = log(-log(abs(x1-boun)));
                    fun = 2*pp + A(2) + A(1)*exp(-A(3)*pp);
                    obj1 = sum((fun-log(abs((log(h(w)))))).^2); % objective function
                    lsqproblem0 = optimproblem("Objective",obj1); % optimization problem
                    if j == 1
                        x0.A = [2,-0.8,1]; % initial point of a
                    else
                        x0.A = sol0.A';
                    end
                    if repeat == 1
                        x0.A = [unifrnd(0,10),unifrnd(-10,0),unifrnd(0,10)];
                        fval_tol = fval_tol*2;
                    end
                    % show(lsqproblem)
                    [sol0,fval0] = solve(lsqproblem0,x0,'Options',opts);
                    repeat = 1;
                    if display_on
                        disp(fval0)
                    end
                end
                Aloglog(:,j) = sol0.A; % record fit parameter
                % plot the result
                findind = find(abs(h) > 1e-6);
                loglogx = log(-log(xx(2)-boun)); % the value of intersection chosen for log-log|x-x0|
                p = log(-log(abs(x(w)-boun)));
                if plot_on
                    plot(p,log(-log(h(w))),'.-'); hold on;
                    plot(p,evaluate(fun,sol0),'.-');
                    plot(p,2*p+log(-1/(2*log(real_lambda))),'.-');
                    scatter(loglogx,2*loglogx+sol0.A(2)+sol0.A(1)*exp(-sol0.A(3)*loglogx));
                    scatter(loglogx,2*loglogx+log(-1/(2*log(real_lambda))));
                    grid on
                    hold off
                    xlabel('log-log(|x-x_0|)')
                    ylabel('log-log(h)')
                    legend('data','best fit','actual line','Location','best')
                end
                aa = sol0.A(1)*exp(-sol0.A(3)*loglogx)+ sol0.A(2); % calculate the y-intercept for y=2x+c
                if display_on
                    disp(fval0)
                end
            end
            % record error
            if exist('fval0','var')==1,fval_all(j) = fval0; else, fval_all(j) = fval1; end
            % title(strcat('a = ',num2str(a_para),', length of edge data: ',num2str(length(w))))
            if writevideo, F(i) = getframe(gcf);writeVideo(vidfile1,F(i)); end
            if not(brute || brutequad ||brutequadnolog)
                lambda = round(exp(-1/2*exp(-aa)),3);
                % figure(2)
                % plot(xx,h,"Color","blue"); hold on; % plot stationary distribution
                % hold on
                % plot(x1,h(w),'.')
                % hold off
            end
            lambda_all(1,j) = lambda; % record lambda estimate

            
            if boun_compare
                bounmin_all(boun_count,j,jj) = bounmin; % record data boundary
                lambda_boun_all(boun_count,j,jj) = lambda; % record lambda estimate
            else
                bounmin_all(1,j) = bounmin; % record data boundary
            end
            % 
            % tic;
            
            % 
            % %nonhyper case
            % A2 = optimvar('A2',2);
            % fun2 = A2(1)*exp(A2(2)*log(abs(x1-bounmin))./abs(x1-bounmin)); % nonhyper tail
            % obj2 = sum((fun2-h(w)).^2); % objective function
            % lsqproblem2 = optimproblem("Objective",obj2); % optimization problem
            % x0.A2 = [1,0.0001]; % initial point of a
            % [sol2,fval2] = solve(lsqproblem2,x0);
            % 
            % fitted_nonhyper = evaluate(fun2,sol2); % evaluate function at solution
            % plot(x1,fitted_nonhyper,"Color",'black') % plot nonhyper best fit
            % text(x1(end),fitted_nonhyper(end),num2str(lambda));
            % if disp_on
            %     disp(a_para)
            %     disp(sol1.A1)
            %     disp(sol2.A2)
            % end
            % disp(num2str((fval1-fval2)/fval2))
            % legend('stationary dist','hyper fit','nonhyper fit','location','northwest')
            % title(strcat('a = ',num2str(a_para)))
            % xlim([c-r,c+r])
            % ylim([0 10])
            if plot_on
            hold off
            pause(.01)
            end
            % toc;
        % end % for boun_compare
            if boun_compare 
            boun_vec_all(:,j) = boun_vec; 
            end
        end
    if writevideo, close(vidfile1);end
    % plot lambda estimates and the real lambda
    if plot_on
    figure(3)
    if paramchange == 1
        hh1 = floor(linspace(single_traj,floor(nn/single_traj)*single_traj,floor(nn/single_traj)));
        hh = floor(linspace(single_traj/2,nn-single_traj/2,floor(nn/single_traj)));
        if mapType{1}(1)~='1'
            yyaxis right
            plot(rr(hh),s(hh1)); hold on;
            ylabel('standard deviation')
            xx = ones(1,nn);
            for iiii =1:nn
                xx(iiii) = x_all{iiii}{1}(1);
            end
            [~,index1]=find(isinf(xx)||isnan(xx)); index1 = unique(index1);
            % g = find((index1(2:end)-index1(1:end-1))>ceil(nn/10));
            % index2 = index1(g(1)); index3 = index1(g(1)+1);
            index2 = index1(1);
            plot([rr(index2),rr(index2)],[0,1],'k-');
            % plot([xtk(index3),xtk(index3)],[0,1],'k-'); hold off;
            yyaxis left
        end
        plot(rr(hh),lambda_all(hh1))
        ylabel('lambda estimates')
        if mapType{1}(1)~='1'
            continue
            legend('lambda estimates','standard deviation')
        end
        if calc_real_lambda
            hold on
            plot(rr(hh),real_lambda_all(hh1))
            plot(rr(hh),target_lambda_all(hh1))
            legend('lambda estimates','real lambda','target lambda','Location','Best')
            hold off
        end
    else
        plot(rr,lambda_all)
        if mapType{1}(1)~='1'
            yyaxis right
            plot(rr,s); hold on;
            ylabel('standard deviation')
            xx = ones(1,nn);
            for iiii =1:nn
                for iiiii = 1:single_traj
                    if x_all{iiii}{iiiii}{1}(1) < xx(iiii) || isnan(x_all{iiii}{iiiii}{1}(1))
                        xx(iiii) = x_all{iiii}{iiiii}{1}(1);
                    end
                end
            end
            [~,index1]=find(isinf(xx)|isnan(xx)); index1 = unique(index1);
            % g = find((index1(2:end)-index1(1:end-1))>ceil(nn/10));
            % index2 = index1(g(1)); index3 = index1(g(1)+1);
            index2 = index1(1);
            plot([rr(index2),rr(index2)],[0,1],'k-');
            % plot([xtk(index3),xtk(index3)],[0,1],'k-'); hold off;
            yyaxis left
        end
        if calc_real_lambda
            hold on
            plot(rr,real_lambda_all)
            plot(rr,target_lambda_all)
            legend('lambda estimates','real lambda','target lambda','Location','Best')
            hold off
        end
    end
    
    grid on
    xlabel('a'); ylabel('lambda')
    axis([min(param),max(param),0,1])
    
    if average ~= 1
        title(strcat('Lambda estimate for {}',num2str(sample_size),' data, data portion: ',num2str(data_portion),', Bin:',num2str(BinNo)))
    else
        title(strcat('Lambda estimate for {}',num2str(single_traj),' data, data portion: ',num2str(data_portion),', Bin:',num2str(BinNo)))
    end
    end
    % drawnow
    disp(['sum of error squared: ',num2str(norm(lambda_all-real_lambda_all))])
    disp(['Number of data:',num2str(length(x)),', BinNo: ',num2str(BinNo),'data portion',num2str(data_portion)])
    disp('Next sample')
    disp(jj)
    lambda_record = lambda_record + lambda_all;
    lambda_each{jj} = lambda_all;
    norm_all(jj,jjj,jjjj) = norm(lambda_all-real_lambda_all)/sqrt(length(rr));
    end
    lambda_record = lambda_record/nsample;
    if average == 0 && nsample ~= 0
        lambda_data = zeros(nsample,n_n);
        for i = 1:nsample
            lambda_data(i,:) = lambda_each{i};
        end
        
        xl = {};
        for i = 1:n_n
            if mod(i,10) == 1
                xl{i} = num2str(rr(i));
            else
                xl{i} = '';
            end
        end
        if plot_on_large
        figure
        boxplot(lambda_data,rr)
        hAx=gca;                                   % retrieve the axes handle
        xtk=hAx.XTick;                             % and the xtick values to plot() at...
        hold on
        plot(xtk,lambda_record)
        plot(xtk,real_lambda_all)
        xticklabels(xl)
        xlabel('$a$','Interpreter','latex')
        ylabel('$\lambda$','Interpreter','latex')
        legend('Averages of $\lambda$ estimates', 'True $\lambda$ values','Interpreter','latex')
        title(strcat('$\lambda$ estimates for$\;$',num2str(single_traj),' data, data portion: ',num2str(data_portion),', bins:',num2str(BinNo)),'Interpreter','latex')
        % exportgraphics(gcf,strcat('Lambda_estimates_for_',num2str(single_traj),'_data_data portion_',num2str(data_portion),'_bins_',num2str(BinNo),'.png'),'Resolution',300)
        end
    end
    disp(strcat('large step ',num2str(jjj)))
    lambda_all_all{jjj} = lambda_record;
    
    jjj = jjj + 1;
end
jjjj = jjjj +1;
end
toc;

%% Method from Kuehn, Malavolta, Rasmussen 
portion = 0.05;
fixportion = 1; % whether to fix the portion of data for choosing sigma1 and sigma2
sigma1 = 6e-2; sigma2 = sigma1; sigma0 = 1e-4;
lambda_all = ones(1,nn); % record all lambda
real_lambda_all = ones(1,nn); % record all real lambda
real_boun_all = ones(1,nn); % record all real boundary
target_lambda_all = ones(1,nn); % record all target lambda
lambda_each = {}; % record lambda for nsample jj
lambda_record = zeros(1,nn); % for adding sample
calc_real_lambda = 1;
if average ~= 1
    single_traj = nsample1;
    nsample = sample_size;
end
figure
for jj = 1:nsample
    for j = 1:nn % different parameters
        % disp(j)
        a_para = rr(j); % retrieve a parameter value
        param(j) = a_para; % record parameter

        % retrieve the time-series data
        if paramchange == 1
            if j < single_traj
                continue
            end
        end
        if average == 1
            xx = x_all{j}{jj}{1}; 
        else
            xx = zeros(1,single_traj); % compile each trajectories
            for ii = 1:single_traj
                if paramchange == 1
                    xx(ii) = x_all{j}{1};
                else
                    xx(ii) = x_all{j}{ii}{1}(jj);
                end
            end
        end

        % take the first coordinate for 2-D data
        if d2==1
            xx = xx(:,1);
        end

        % retrieve the mapping for derivative calculations
        if mapType{1}(1)=='1'
            f = @(x) 3*((1-exp(-x))./(1+exp(-x))) - a_para; % the mapping
        elseif mapType{1}(1)=='x'
            f = @(x) 3*((1-exp(-x))./(1+exp(-x))) - a_para - sigma*abs(x); % the mapping
        else
            f = @(x) x+a_para-x^2;
        end
        syms yy; fprime = matlabFunction(diff(f(yy))); % derivative of the map

        x = sort(xx);
        if fixportion == 1
            sigma1 = x(ceil(length(xx)*portion))-min(xx); sigma2 = sigma1;
        end
        m1 = min(xx); m2 = m1 + sigma1 + sigma0;
        stack1 = find(xx(1:end-1)<=m1+sigma1 & xx(1:end-1)>=m1);
        stack2 = find(xx(1:end-1)<=m2+sigma2 & xx(1:end-1)>=m2);
        x1 = xx(stack1+1); x2 = xx(stack2 +1);
        lambda = (sum(x2)/length(x2)-sum(x1)/length(x1))/(m2-m1+sigma0);
        lambda_all(1,j) = lambda; % record lambda estimate
        bounmin = min(xx);
        if calc_real_lambda && abs(bounmin)<10
            bounm = fzero(@(x)(f(x)-sigma-x),bounmin); % find the actual boundary
            % find the boundary around the estimates
            if abs(bounm-bounmin)>1
                for bountry = linspace(bounmin-1,bounmin+1,1000)
                    bounm = fzero(@(x)(f(x)-sigma-x),bountry);
                    if abs(bounm-bounmin)<1
                        break
                    end
                end
            end
            real_lambda = fprime(bounm); % real lambda
            target_lambda = fprime(bounmin); % for 1-exp(-x)/(1+exp(-x)) S map
            if linear_map == 1
                real_lambda = a_para; % for linear map
                target_lambda = 1-(sigma*(1-a_para))/abs(bounmin); % target lambda from boundary observed 
            end
            real_boun_all(1,j) = bounm;
            real_lambda_all(1,j) = real_lambda;
            target_lambda_all(1,j) = target_lambda;
        end
        % disp(lambda)
    end
    if paramchange == 1
        hh1 = floor(linspace(single_traj,floor(nn/single_traj)*single_traj,floor(nn/single_traj)));
        hh = floor(linspace(single_traj/2,nn-single_traj/2,floor(nn/single_traj)));
        plot(rr(hh),lambda_all(hh1))
        hold on
        plot(rr(hh),real_lambda_all(hh1))
        hold off
    else
        plot(rr,lambda_all)
        hold on
        plot(rr,real_lambda_all)
        hold off
        % plot(rr,target_lambda_all)
    end
    % legend('estimates','real','target')
    ylim([0,1])
    disp('Next sample')
    disp(jj)
    lambda_record = lambda_record + lambda_all;
    lambda_each{jj} = lambda_all;
    drawnow
end
lambda_record = lambda_record/nsample;
if average == 0
    figure
    lambda_data = zeros(nsample,n_n);
    for i = 1:nsample
        lambda_data(i,:) = lambda_each{i};
    end
    % to set nice x ticks (otherwise too much ticks)
    xl = {};
    for i = 1:n_n
        if mod(i,10) == 1
            xl{i} = num2str(rr(i));
        else
            xl{i} = '';
        end
    end
    boxplot(lambda_data,rr)
    boxplot(lambda_data,rr)
    hAx=gca;                        % retrieve the axes handle
    xtk=hAx.XTick;                  % and the xtick values to plot()
    % if mapType{1}(1)~='1'
    %     yyaxis right
    %     boxplot(s,rr); hold on;
    %     ylabel('standard deviation')
    %     [~,index1]=find(isnan(s)); index1 = unique(index1);
    %     g = find((index1(2:end)-index1(1:end-1))>ceil(nn/10));
    %     index2 = index1(g(1)); index3 = index1(g(1)+1);
    %     plot([xtk(index2),xtk(index2)],[0,1],'k-');
    %     plot([xtk(index3),xtk(index3)],[0,1],'k-'); hold off;
    %     yyaxis left
    % end
    hold on
    plot(xtk,lambda_record) % plot the average lambda approximations
    plot(xtk,real_lambda_all) % plot real lambda values
    xticklabels(xl)
    xlabel('$a$','Interpreter','latex')
    ylabel('$\lambda$','Interpreter','latex')
    legend('Averages of $\lambda$ estimates', 'True $\lambda$ values','Interpreter','latex')
end

%%
figure
boxplot(lambda_data_2,rr,'Colors','b','Symbol','+b'); hold on
boxplot(lambda_data_3,rr,'Colors','k','Symbol','+k'); hold on
boxplot(lambda_data_1,rr,'Colors','g','Symbol','+g'); hold on
plot(xtk,lambda_record_2,'b') % plot the average lambda approximations
plot(xtk,lambda_record_3,'k') % plot the average lambda approximations
plot(xtk,lambda_record_1,'g') % plot the average lambda approximations
plot(xtk,real_lambda_all,'r')
ylim([0,1])
legend('$a_2y^2 + a_1y, y = \log(x-x_-)$','$a_2(y^2-2ylog(-y)) + a_1y, y = \log(x-x_-)$','Interval method','True $\lambda$ values','Interpreter','latex')
xticklabels(xl)
xlabel('$a$','Interpreter','latex')
ylabel('$\lambda$','Interpreter','latex')
title('$\lambda$ estimations of $tanh$ map with uniform bounded noise, estimated boundary $x_-$','Interpreter','latex')

figure
boxplot(abs(lambda_data_1-real_lambda_all),rr,'Colors','b','Symbol','+b'); hold on
boxplot(abs(lambda_data_2-real_lambda_all),rr,'Colors','k','Symbol','+k'); hold on; 
boxplot(abs(lambda_data_9-real_lambda_all),rr,'Colors','g','Symbol','+g')
plot(xtk,abs(lambda_record_1-real_lambda_all),'b') % plot the average lambda approximations
plot(xtk,abs(lambda_record_2-real_lambda_all),'k');
plot(xtk,abs(lambda_record_9-real_lambda_all),'g') % plot the average lambda approximations
xticklabels(xl)
xlabel('$a$','Interpreter','latex')
ylabel('$|\lambda-\lambda_{estimated}|$','Interpreter','latex')
legend('$a_2y^2 + a_1y, y = \log(x-x_-)$','$a_2(y^2-2ylog(-y)) + a_1y, y = \log(x-x_-)$','Interval method','Interpreter','latex')
title('Errors of different methods for $tanh$ map, uniform noise, actual boundary $x_-$','Interpreter','latex')

