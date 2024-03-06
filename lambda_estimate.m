% Define different parameters and simulate data
tic; close all;
nsample = 1000; % length of time series data
sigma = .47; % bound of noise or 2*std
h = 0.01; % time step for continuous time
n_avg = 0; % the number of iteration to ensure convergence
transient = 100; % initial amount of iteration to get rid of transient
nn = 101; % number of parameters to run
% use different sample size for average == 1 
sample_size = 10; % number of different samples of time series data
x_all = {}; % to record all x values
limitx = ones(2,nn); % record the limit values of data
d2 = 0; % if use two dimensional map
plot_on = 1; % whether to plot the histogram
mapType = "Uniform"; % type of deterministic map "1duni"
negativestart = 0; % whether to start at negative initial point if less than 0
discrete = 1; %whether the system used is discrete
average = 0; % whether to use lots of points, instead of one traj
continuous_change = 0; % whether to use initial data continuously from previous parameter
nsample1 = nsample; % for single trajectory
rr = linspace(0.3,1.2,nn);

% initialise matrix for recording standard deviation
if d2==1
    if average == 0
        s = ones(2,sample_size,nn);
    else
        s = ones(2,nsample,nn);
    end
else
    if average == 0
        s = ones(sample_size,nn);
    else
        s = ones(nsample,nn);
    end
end

j = 1; % for different sample size
for n = sample_size % goes through different time series
    for k = 1:nn % number of parameters
        r = rr(k); % parameter
        if d2 == 0 % one-dimensional map
            % specify different intial point for different type of map 
            if mapType{1}(1)=='1' 
                x0 = 1;
            else
                x0 = 1;
            end
                
            % restart initial point i if continuous change is on or for first iteration
            if continuous_change == 0 || k == 1 
                x = ones(1,n)*x0; % start all points from x0
                for i = 1:transient % run the system a few round to ensure it is in minimal invariant set
                    x = saddle_discrete(x,r,sigma,"dist",mapType);
                end
                if any(x<0) && negativestart && mapType{1}(1)=='1' % restart for points going over to negative side
                    x0 = -1;
                    x = ones(1,n)*x0; % start all points from x0
                    for i = 1:transient
                        x = saddle_discrete(x,r,sigma,"dist",mapType);
                    end
                end
            end
            if average == 0
                xx = ones(nsample,n); % recording single traj data for calc std
            end

            % start the iteration of mappings
            for i = 1:n_avg+nsample
                x = saddle_discrete(x,r,sigma,"dist",mapType);
                if i>n_avg % take the next few iterations as different samples
                    x_all{k}{i-n_avg}{j} = x;

                    % record the standard deviations
                    if average == 1
                        s(i-n_avg,k) = std(x);
                    else
                        xx(i-n_avg,:) = x;
                    end
                end
            end
            if average == 0
                s(:,k) = std(xx)';
            end

        else % two-dimensional map
            x0 = 1; 
            x = ones(n,2)*x0; % start all points from x0
            if average == 0
                xx = ones(2,nsample,n); % recording single traj data for calc std
            end
            for i = 1:n_avg+nsample
                rrrr = unifrnd(0,1,[n,1]);
                theta = unifrnd(-pi,pi,[n,1]);
                % randnum = normrnd(0,sigma,n,2); % normally distributed
                randnum = sigma*rrrr.*[cos(theta),sin(theta)]; % uniformly distributed
                % x = henon(x,r,0.3,2,0); % deterministic
                x = henon(x,r,0.3,2,0)+randnum;
                if i>n_avg % take the time series as data
                    x_all{k}{i-n_avg}{j} = x;

                    % record standard deviation
                    if average == 1
                        s(:,i-n_avg,k) = std(x)';
                    else
                        xx(:,i-n_avg,:) = x';
                    end
                end
                if average == 0
                    s(1,:,k) = std(xx(1,:,:))';
                    s(2,:,k) = std(xx(2,:,:))';
                end
            end
        end

        % plot the histogram data (only one set of data for multiple samples)
        if plot_on
            if average == 1
                histogram(x)
            else
                xx = ones(1,nsample);
                for ii = 1:nsample
                    xx(ii) = x_all{k}{ii}{1}(1);
                end
                histogram(xx)
            end
            title(num2str(r))
            drawnow
        end
        disp(k)
    end
    disp(n)
    j = j + 1;
end

%% Find the optimal tail distribution from the stationary distribution using least square model
tic; clear x0
writevideo = 0; % whether to record all graphical output in a video
if writevideo
    vidfile1 = VideoWriter('exponential_loglog','MPEG-4');
    vidfile1.FrameRate = 2;
    open(vidfile1);
end

% important hyperparameters: please change to suitable bin (less bins for less data)
BinNo = 200; % the number of bin in histogram
data_portion = 0.1; % portion of data point to fit

% parameter values and settings
brute = 1; brutequad = 0; % which method to use
disp_on = 0; % whether to display result
plot_on = 1; % whether to plot the figures
display_on = 0; % whether to display fval and other output related message
weight_on = 0; weightpower = 1; % whether to introduce weight to the least square objective function
smooth_on = 0; % whether to smoothen the histrogram data
fval_tol = 3e-1; % tolerance of function evaluation
left_data = 1; % whether to test left or right data
uneven = 0; % whether to use uneven bins
if uneven == 1
    k = (atan(linspace(-1.3,1.3,BinNo+1))/atan(1.3)+1)/2; % uneven histogram
end
calc_real_lambda = 1; % whether to calculate real lambda

% initialise matrices and vectors
lambda_all_all = {}; % record all lambda
lambda_each = {}; % record lambda for nsample jj
norm_all = zeros(1,length(sample_size)); % record the norm of real lambda and estimated
bounmin_all = ones(1,nn); % record data boundary
lambda_all = ones(1,nn); % record lambda estimate
real_boun_all = ones(1,nn); % record real boundary point
real_lambda_all = ones(1,nn); % record real lambda
target_lambda_all = ones(1,nn); % record targeted lambda given data boundary
param = ones(1,nn); % record parameter
A_all = ones(2,nn); % record brute param estimate
fval_all = ones(1,nn); % record error
hist_all = {}; % record all the histogram data
jjj = 1;

% for linear map with bounded noise
linear_map = 0; % deterministic map is linear
linear_cheat = 0; % whether to have known boundary

% start the approximations
for BinNo = 200
% for BinNo = linspace(10,500,50) % run through different number of bins
% for data_portion = linspace(0.01,0.8,50) % run different data portions
    lambda_record = zeros(1,nn); % for adding sample
    if average ~= 1
        single_traj = nsample1;
        nsample = sample_size;
    end

    % run through different sample of time-series
    for jj = 1:nsample
        for j = 1:nn % different parameters
            % disp(j)
            a_para = rr(j); % retrieve a parameter value
            param(j) = a_para; % record parameter

            % retrieve the time-series data
            if average == 1
                xx = x_all{j}{jj}{1}; 
            else
                xx = zeros(1,single_traj); % compile each trajectories
                for ii = 1:single_traj
                    xx(ii) = x_all{j}{ii}{1}(jj);
                end
            end

            % take the first coordinate for 2-D data
            if d2==1
                xx = xx(:,1);
            end

            % retrieve the mapping for derivative calculations
            if mapType{1}(1)=='1'
                f = @(x) 3*((1-exp(-x))./(1+exp(-x))) - a_para; % the mapping
            else
                f = @(x) x+a_para-x^2;
            end
            figure(1)
            if smooth_on
                [h,x]=ksdensity(xx,'NumPoints',BinNo); width = x(2)-x(1); % smoothen the stationary dist curve
            else
                if uneven
                    hh = histogram(xx,min(xx)+(max(xx)-min(xx))*kk); 
                else
                    hh = histogram(xx,BinNo);
                end
                title(['a = ',num2str(a_para)]); 
                x = hh.BinEdges; h = hh.Values; width = hh.BinWidth; % take data of the distribution
                if all(width == 'nonuniform'), width = x(2)-x(1);end
                x = (x(2:end)+x(1:end-1))/2; % take center of the bins as x-data
            end

            % heights of histogram
            if uneven == 1
                h = h/sum(h.*(kk(2:end)-kk(1:end-1)));
            else
                h = h/norm(h,1)/width; % so that area is 1
            end
            hist_all{j} = [x;h]; % record the data
                
            % get the leftmost or rightmost data as boundary estimation
            if left_data == 1
                bounmin = min(x)-width; % the lower boundary from hist
                if linear_cheat
                    bounmin = -sigma; % for linear map fixed boundary simulation
                end
            else
                bounmin = max(x)+width; % the upper boundary from hist
            end
            syms yy; fprime = matlabFunction(diff(f(yy))); % derivative of the map
        
            % determine which data to include in approximations
            extol = max(h)*1e-2; % exclude data where h value is small
            % take data within data_portion, h value sufficiently large (avoid infinite when log) and not too far from edge
            if uneven
                if left_data == 1
                    w = find(cumsum(h.*(kk(2:end)-kk(1:end-1)))<data_portion & abs(h)>extol & log(x-bounmin)<0); % only test some data portion
                else
                    w = find(cumsum(h.*(kk(2:end)-kk(1:end-1)))>1-data_portion & abs(h)>extol & log(abs(x-bounmin))<0); % only test some data portion
                end
            else
                if left_data == 1
                    w = find(cumsum(h)*width<data_portion & abs(h)>extol ); % only test some data portion & log(x-bounmin)<0
                else
                    w = find(cumsum(h)*width>1-data_portion & abs(h)>extol & log(abs(x-bounmin))<0); % only test some data portion
                end
            end
            if left_data == 1
                w = w(2:end); % get rid of end point 
            else
                w = w(1:end-1); % get rid of end point 
            end
            x1 = x(w); % very tail of the distribtuion to test
            
            % calculate real lambda from derivative of function
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
                % bounm = fzero(@(x)(f(x)-sigma-x),1);
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

            opts=optimset('Display','off');
            if isempty(x1),continue;end % if no data then proceed to next

            % start of approximation
            if brute % optimzation of the asymptotic equation
                opts=optimset('Display','off');
                figure(2)
                plot(x,h,"Color","blue"); hold on; % plot stationary distribution
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
                if plot_on
                    hhh = plot(x1,fitted_hyper,"Color",'red'); % plot hyper best fit
                    title(strcat('a = ',num2str(a_para)))
                    % xlim([c-r,max(x1)])
                end
                lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                hold off
                if display_on
                    disp(fval1)
                end
            elseif brutequad == 1
                opts=optimset('Display','off');
                figure(2)
                % plot(x,h,"Color","blue"); hold on; % plot stationary distribution
                repeat = 0; % to know whether need to recalculate
                fval1 = 10; % initialise large fval
                count = 1;
                while fval1 > fval_tol % repeat if tolerance not met
                    A1 = optimvar('A1',2);
                    fun = log(abs(x1-bounmin));
                    fun1 = A1(1)*fun + A1(2)*fun.*fun; % hyper tail
                    weight = 1./linspace(1,length(h(w)),length(h(w))).^weightpower ; % to apply weightage for points near boundary is priotised
                    if weight_on
                        obj1 = sum(weight.*(fun1-log(h(w))).^2)/sum(weight); % objective function
                    else
                        obj1 = sum((fun1-log(h(w))).^2); % objective function
                    end
                    constraint = A1<=0;
                    % lsqproblem1 = optimproblem("Objective",obj1); % optimization problem
                    lsqproblem1 = optimproblem("Objective",obj1,"Constraints",constraint); % optimization problem
                    if count==1
                        x0.A1 = [1,-1/10]; % initial point of a
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
                logx = log(abs(x(2:end)-bounmin));
                p = log(abs(x(w)-bounmin));
                if plot_on
                    hhh = plot(p,log(h(w)),'.-'); hold on;
                    plot(p,evaluate(fun1,sol1),'.-');
                    title(strcat('a = ',num2str(a_para),', length of edge data: ',num2str(length(w))))
                    % disp(fval1)
                end
                lambda = round(exp(1/2/sol1.A1(2)),3); % lambda estimate
                hold off
            end

            % record error
            if exist('fval0','var')==1,fval_all(j) = fval0; else, fval_all(j) = fval1; end
            title(strcat('a = ',num2str(a_para),', length of edge data: ',num2str(length(w))))
            if writevideo, F(i) = getframe(gcf);writeVideo(vidfile1,F(i)); end
            lambda_all(1,j) = lambda; % record lambda estimate
            bounmin_all(1,j) = bounmin; % record data boundary
            hold off
            pause(.01)
            % toc;
        end
        if writevideo, close(vidfile1);end
    
        % plot lambda estimates and the real lambda
        figure
        plot(param, lambda_all,'.-'); 
        if calc_real_lambda
            hold on
            plot(param, real_lambda_all)
            plot(param, target_lambda_all,'.-')
            legend('lambda estimates','real lambda','target lambda','Location','Best')
            hold off
        end
        grid on
        xlabel('a'); ylabel('lambda')
        axis([min(param),max(param),0,1])
        disp(['sum of error squared: ',num2str(sum(abs(lambda_all-target_lambda_all)))])
        disp(['Number of data:',num2str(length(x)),', BinNo: ',num2str(BinNo),'data portion',num2str(data_portion)])
        if average ~= 1
            title(strcat('Lambda estimate for {}',num2str(sample_size),' data, data portion: ',num2str(data_portion),', Bin:',num2str(BinNo)))
        else
            title(strcat('Lambda estimate for {}',num2str(single_traj),' data, data portion: ',num2str(data_portion),', Bin:',num2str(BinNo)))
        end
        drawnow
        disp('Next sample')
        disp(jj)
        lambda_record = lambda_record + lambda_all;
        lambda_each{jj} = lambda_all;
    end

    % plot the box plot for different sets of time series
    if average == 0
        figure
        lambda_data = zeros(100,101);
        for i = 1:nsample
            lambda_data(i,:) = lambda_each{i};
        end
        % to set nice x ticks (otherwise too much ticks)
        xl = {};
        for i = 1:101
            if mod(i,10) == 1
                xl{i} = num2str(rr(i));
            else
                xl{i} = '';
            end
        end
        boxplot(lambda_data,rr)
        hAx=gca;                        % retrieve the axes handle
        xtk=hAx.XTick;                  % and the xtick values to plot()
        if mapType{1}(1)~='1'
            yyaxis right
            boxplot(s,rr); hold on;
            ylabel('standard deviation')
            [~,index1]=find(isnan(s)); index1 = unique(index1);
            g = find((index1(2:end)-index1(1:end-1))>ceil(nn/10));
            index2 = index1(g(1)); index3 = index1(g(1)+1);
            plot([xtk(index2),xtk(index2)],[0,1],'k-');
            plot([xtk(index3),xtk(index3)],[0,1],'k-'); hold off;
            yyaxis left
        end
        
        hold on
        plot(xtk,lambda_record) % plot the average lambda approximations
        plot(xtk,real_lambda_all) % plot real lambda values
        xticklabels(xl)
        xlabel('a')
        ylabel('lambda')
        legend('Averages of lambda estimates', 'True values of lambda')
        title(strcat('Lambda estimates for {}',num2str(single_traj),' data, data portion: ',num2str(data_portion),', bins:',num2str(BinNo)))
        % exportgraphics(gcf,strcat('Lambda_estimates_for_',num2str(single_traj),'_data_data portion_',num2str(data_portion),'_bins_',num2str(BinNo),'.png'),'Resolution',300)
    end
    disp(strcat('large step ',num2str(jjj)))
    
    lambda_record = lambda_record/nsample;
    lambda_all_all{jjj} = lambda_record;
    norm_all(jjj) = norm(lambda_record-real_lambda_all);
    jjj = jjj + 1;
end
toc;
