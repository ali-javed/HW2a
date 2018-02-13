function comparederivs(allx,f,truedf, true2df, true3df)
% COMPAREDERIVS: % Compares true, 1st-order, and 2nd-order slope approximations
% function comparederivs(allx,f,truedf, true2df, true3df)
%
% INPUTS:
%   allx: is a scalar or vector of what domain values to look at
%   f: is a function handle of the function to differentiate
%   truedf: is a function handle for the analytic derivative function of f
%   true2df: is a function handle for the analytic 2nd derivative function of f
%   true3df: is a function handle for the analytic 3rd derivative function of f
%
% OUTPUTS:
%   This program produces one plot for each value in the vector allx
%   The plots show the error in the approximations as a function of step size
%
% sample call: comparederivs([0 5],@exp, @exp, @exp, @exp)

% Maggie Eppstein, 2/10/08; documentation improved 2/15/11
% Edited by Ali Javed and Zoe Koch 2/12/18

% THIS IS A GOOD EXAMPLE OF A WELL-DOCUMENTED FILE;  NOTE THE FOLLOWING:
%   a) contents and consistent organization of function headers; first
%   comment line for LOOKFOR command, 1st contiguous comment block for HELP
%   command; always define inputs/outputs, including size constraints or
%   other pre-/post-conditions;
%   b) in-line comments should always be at one level of abstraction higher
%   than the code itself;
%   c) use of full-line UPPER-CASE in-line comments to give a high-level
%   description of what each logically-related code block does; you can read through these alone to
%   get a good understanding of what the code does, without even looking at
%   the code;
%   d) additional lower-case comments at the ends of potentially confusing
%   lines for clarification;


% FOR EACH X-VALUE SPECIFIED, PLOT THE APPROXIMATION ERRORS AS A FUNCTION OF STEPSIZE
h=logspace(-20,-1,20); %logarithmically-spaced step sizes to try
for xi=1:length(allx) %each x gets its own plot
    x=allx(xi);
    
    % FIND OPTIMAL H-VALUES BY ITERATING THROUGH THE APPROXIMATION FUNCTION
    % UNTIL WE'VE CONVERGED TO WITHIN A CERTAIN TOLERANCE 
    %backward differencing method
     h_opt_backwards = nthroot(eps,3);
     
     %start loop
     while true
         
         h_temp = h_opt_backwardsFunc(x,true2df, h_opt_backwards);
         %if difference is less than a threshold break
         if h_temp-h_opt_backwards < .0000001
             break
         end
         h_opt_backwards = h_temp;
         
     end
     
     %for central differencing method
     h_opt_central = nthroot(eps,3);
     
     while true
         h_temp = h_opt_centralFunc(x,true3df, h_opt_central);
         %if difference is less than a threshold, break
         if h_temp-h_opt_central < .0000001
             break
         end
         h_opt_central = h_temp;
         
     end
      
    
    
   

    % COMPUTE TRUE DERIVATIVE AND ITS APPROXIMATIONS
    df=truedf(x); % compute true derivative at x
    df1=backdiff(x,f,h); %approximate with 1st order backwards difference
    df2=centraldiff(x,f,h); %approximate with 2nd order central difference
    
    % COMPUTE ANALYTICAL ERROR ESTIMATES
    analytical_backward_err = analytical_BKW(x,h,true2df);
    analytical_ctrl_err = analytical_CN(x,h,true3df);
    
    % COMPUTE APPROXIMATION ERRORS BY COMPARING TO TRUE DERIVATIVES
    err1=abs(df-df1);
    err2=abs(df-df2);
    err3=abs(analytical_backward_err);
    err4=abs(analytical_ctrl_err);

    % PLOT THE APPROXIMATION ERRORS AS A FUNCTION OF STEPSIZE
    fig = figure
    
    loglog(h,err1,'ro',h,err2,'b*',h,err3,'g+',h,err4,'mx');
    
    % LINEAR REGRESSION OF LOG-LOG RELATIONSHIPS (ONLY IN REGION GOVERNED BY TRUNCATION ERROR)
    coef1=polyfit(log(h(end-4:end)),log(err1(end-4:end)),1);
    coef2=polyfit(log(h(end-4:end)),log(err2(end-4:end)),1);
    
    hold on % add the best-fit lines to the plot
    loglog(h(end-4:end),exp(coef1(2))*h(end-4:end).^coef1(1),'r-');
    loglog(h(end-4:end),exp(coef2(2))*h(end-4:end).^coef2(1),'b-');
    
    % LABEL THE PLOT
    set(gca,'fontsize',14) % be kind to the instructor's aging eyes!
    xlabel('h')
    ylabel('error')
    legend('back diff','central diff','analytical backward','analytical central',...
        ['slope = ',num2str(coef1(1),4)],['slope = ',num2str(coef2(1),4)],...
        'Location','BestOutside');
    title(['at x = ',num2str(x)])
    
    %create file name and save the image
    f_name = sprintf('../figures/fig%0.2f.png',x)
    saveas(fig,f_name)
    
    %ALLOW USER TO VIEW EACH PLOT BEFORE MOVING ON TO THE NEXT
    if xi<length(allx)
        disp('Hit any key to continue...')
        pause 
    end
end
figure(gcf)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: the following functions are placed here for convenience for this demo code, 
% but cannot be called from outside this file; in general, these should be
% in their own files (e.g. in a directory for your personal library "toolbox" that 
% you add to the Matlab path to access your own handy utility functions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df1=backdiff(x,f,h) 
% BACKDIFF: 1st order backwards difference approximation to first derivative
% function df1=backdiff(x,f,h) 
%
% INPUTS:
% x: location(s) of where in domain to approximat the derivative
% f: handle of function to approximate derivative of
% h: stepsize(s) to use in approximation
% SIZE CONSTRAINTS: at least one of x or h must be a scalar, but the other
% can be of any dimension (scalar, vector, matrix)
%
% OUTPUTS:
% df1: 1st order approximation to first deriv (slope) of f at x 
%     (same size as largest of x or h)
% 
% SAMPLE CALLS:
%   df1=backdiff(0,@sin,[1e-3 1e-2 1e-1]) % can call with vector of stepsizes
%   df1=backdiff(0:.5:3,@sin,1e-3) % or can call with vector of domain values

% AUTHOR: Maggie Eppstein, 2/15/2011

df1=(f(x)-f(x-h))./h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df2=centraldiff(x,f,h)
% CENTRALDIFF: 2nd order central difference approximation to first derivative
%
% INPUTS:
% x: location(s) of where in domain to approximat the derivative
% f: handle of function to approximate derivative of
% h: stepsize(s) to use in approximation
% SIZE CONSTRAINTS: at least one of x or h must be a scalar, but the other
% can be of any dimension (scalar, vector, matrix)
%
% OUTPUTS:
% df2: 2nd order approximation to first deriv (slope) of f at x 
%     (same size as largest of x or h)
%
% SAMPLE CALLS:
%   df2=backdiff(0,@sin,[1e-3 1e-2 1e-1]) % can call with vector of stepsizes
%   df2=backdiff(0:.5:3,@sin,1e-3) % or can call with vector of domain values

% AUTHOR: Maggie Eppstein, 2/15/2011

df2=(f(x+h)-f(x-h))./(2*h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h_opt_backwards = h_opt_backwardsFunc(x,true2df, h)
% h_opt_backwardsFunc: Finds optimal step size (h) value using backwards differencing
% method i.e. between x-h and x
%
% INPUTS:
% x: location(s) of where in domain to approximat the error
% h: stepsize(s) to use in approximation
% true2df: second derivative of required function
% SIZE CONSTRAINTS:x must be a sclare and h can be either a scalar or a vector
%
% OUTPUTS:
%h_opt_backwards: this this is step size that will provide the maximum
%upper bound error 
% 
% SAMPLE CALLS:
%   h_opt_backwardsFunc(pi/4,@cos, 0.0000001);

% AUTHOR: Zoe Koch and Ali Javed, 2/13/2018
    M2 = max(true2df(linspace(x-h,x,5)));
    h_opt_backwards = 2*sqrt(eps/M2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function h_opt_central = h_opt_centralFunc(x,true3df, h)
% h_opt_centralFunc: Finds optimal step size (h) value using central differencing
% method i.e. between x-h and x+h
%
% INPUTS:
% x: location(s) of where in domain to approximat the error
% h: stepsize(s) to use in approximation
% true2df: second derivative of required function
% SIZE CONSTRAINTS:x must be a sclare and h can be either a scalar or a vector
%
% OUTPUTS:
%h_opt_central: this this is step size that will provide the maximum
%upper bound error 
% 
% SAMPLE CALLS:
%   h_opt_h_opt_central(pi/4,@cos, 0.0000001);

% AUTHOR: Zoe Koch and Ali Javed, 2/13/2018
    M3 = max(true3df(x-h:.001:x+h) );
    h_opt_central = nthroot(12*eps/M3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function analytical_backward_err = analytical_BKW(x,h,true2df)
% analytical_backward_err: Estimates the error using backward differencial
% in the range x, x-h

% INPUTS:
% x: location(s) of where in domain to approximat the error
% h: stepsize(s) to use in approximation
% true2df: second derivative of required function
% SIZE CONSTRAINTS:x must be a sclare and h can be either a scalar or a vector
%
% OUTPUTS:
%analytical_backward_err: The maximum backward difference error
% 
% SAMPLE CALLS:
%   analytical_BKW(pi/4,0.0000001,@cos);

% AUTHOR: Zoe Koch and Ali Javed, 2/13/2018
  

%declare empty array to return
analytical_backward_err = [];
for i=1:length(h)
    
    %h* max(f''(x-a:x))/2 + 2*eps/h
    analytical_backward_err_1 = h(i) * max(true2df(...
    linspace(x-h(i),x,5)))./2 + 2*eps./h(i);  
    
    %save values in vector
    analytical_backward_err = [analytical_backward_err analytical_backward_err_1];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function analytical_central_err = analytical_CN(x,h,true3df)
analytical_central_err = [];
% analytical_central_err: Estimates the error using central differencial
% in the range x-h, x+h

% INPUTS:
% x: location(s) of where in domain to approximat the error
% h: stepsize(s) to use in approximation
% true2df: second derivative of required function
% SIZE CONSTRAINTS:x must be a sclare and h can be either a scalar or a vector
%
% OUTPUTS:
%analytical_backward_err: The maximum backward difference error
% 
% SAMPLE CALLS:
%   analytical_CN(pi/1,0.0000001,@exp);

% AUTHOR: Zoe Koch and Ali Javed, 2/13/2018

for i=1:length(h)
    %h^2 * max(f'''(x-a:x+a))/12 + eps/h)
    analytical_forward_err_1 = h(i).^2 * max(true3df(linspace(x-h(i),x+h(i),5)))./12 + eps./h(i);
     %save values in vector
    analytical_central_err = [analytical_central_err analytical_forward_err_1];
end


