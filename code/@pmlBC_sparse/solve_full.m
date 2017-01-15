function scattValuesVec = solve_full(obj,k,incfun,precond,verbose)
% Solve the scattering problem to get values of scattered wave on mesh of
% B(0,R). 
%  Inputs:
%    - k - wave number
%    - incfun - the incident function OR a vector of -V*incfun values.
%               If neither incfun nor precond passed or incfun = [],  
%               uses default incfun. 
%    - precond - preconditioner function x -> M*x for inverting T(k).
%                Uses best preconditioner available if none is passed or 
%                if precond = [].
%    - verbose - only shows the gmres output if true. Default false.
    if nargin < 5, verbose = 0; end
    if nargin < 3 || isempty(incfun)
        RHS = obj.RHSfromFun_full(k);
    elseif isvector(incfun)
        RHS = incfun;        
    else
        RHS = obj.RHSfromFun_full(incfun);
    end
    if nargin < 4 || isempty(precond)
        precond = @(x) obj.apply_pc_full(x,k);
    end
    
    T = @(x) obj.apply_Tfull(x,k); % solving T(k)*scattValuesVec = RHS
    restart = []; %10; 10 gives terrible results
    tol = 1e-12;
    maxiter = obj.Nt*obj.Nr; % size of matrix
    if verbose
        scattValuesVec = gmres(T,RHS,restart,tol,maxiter,precond);
    else
        [scattValuesVec,flag] = gmres(T,RHS,restart,tol,maxiter,precond);
    end
