function test_Tstar

%% Check that I apply Dr' correctly to vectors

a = rand; b = 2*a; N = 10; x = rand(N+1,1);
[Dr,r] = cheb(N,a,b);

[abserr,relerr] = get_err(Dr*x, chebfft_ah(x,a,b));
fprintf('Error between mult by Dr and chebfft_ah on vector: %4.2e, %4.2e\n', abserr, relerr);

[abserr,relerr] = get_err(Dr'*x, chebfft_star_ah(x,a,b));
fprintf('Error between mult by Dr_star and chebfft_star_ah on vector: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply Dr' correctly to matrices

x = rand(N+1,2);
[abserr,relerr] = get_err(Dr*x, chebfft_ah(x,a,b));
fprintf('Error between mult by Dr and chebfft_ah: %4.2e, %4.2e\n', abserr, relerr);

[abserr,relerr] = get_err(Dr'*x, chebfft_star_ah(x,a,b));
fprintf('Error between mult by Dr_star and chebfft_star_ah: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply Dr_I', Dr_J', Drr_I', and Drr_J' right

Nt = 20; Nrs = [10, 15]; Rs = [1,2];
s        = scattResComp2d       (Nt,Nrs,Rs);
s_sparse = scattResComp2d_sparse(Nt,Nrs,Rs);

x = rand(s.Nr,2);

[abserr,relerr] = get_err(s.Dr_I'*x, s_sparse.apply_Dr_I_star(x));
fprintf('Error between Dr_I adjoint and apply_Dr_I_star: %4.2e, %4.2e\n', abserr, relerr);

[abserr,relerr] = get_err(s.Dr_J'*x, s_sparse.apply_Dr_J_star(x));
fprintf('Error between Dr_J adjoint and apply_Dr_J_star: %4.2e, %4.2e\n', abserr, relerr);

[abserr,relerr] = get_err(s.Drr_I'*x, s_sparse.apply_Drr_I_star(x));
fprintf('Error between Drr_I adjoint and apply_Drr_I_star: %4.2e, %4.2e\n', abserr, relerr);

[abserr,relerr] = get_err(s.Drr_J'*x, s_sparse.apply_Drr_J_star(x));
fprintf('Error between Drr_J adjoint and apply_Drr_J_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply Aorig' correctly

Vs = {@(r,t) 0*r, @(r,t) 0*r + 1};
coords = 'polar';
dtn  = DtNBC       (Nt,Nrs,Vs,coords,Rs);
dtns = DtNBC_sparse(Nt,Nrs,Vs,coords,Rs);

x = rand(s.Nr*s.Nt,1);
VvaluesVec = s.valuesVecFromFunCellArray(Vs,coords);

[abserr,relerr] = get_err(dtn.Aorig'*x, s_sparse.apply_Aorig_star(x,conj(VvaluesVec)));
fprintf('Error between dtn.Aorig adjoint and apply_Aorig_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that dtn.A = dtn.Aorig + kron(I,localPlacement)*ArowChange

I = eye(dtn.Nt);

[abserr,relerr] = get_err(dtn.A, ...
    dtn.Aorig + kron(I,dtns.localPlacement)*dtns.ArowChange);
fprintf('Error between dtn.A and dtn.Aorig + change: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply dtnA' correctly

[abserr,relerr] = get_err(dtn.A'*x, dtns.apply_A_star(x));
fprintf('Error between dtn.A adjoint and apply_A_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply rat.A21' correctly

ell = ellipse(50-10i, 0, 30, 15, [], []);
br = 1;
N = 40;

rat  = ratApproxDtNBC       (dtn ,ell,br,N);
rats = ratApproxDtNBC_sparse(dtns,ell,br,N);

x = rand(size(rat.A21,1),1);

[abserr,relerr] = get_err(rat.A21'*x, rats.apply_A21_star(x));
fprintf('Error between rat.A21 adjoint and apply_A21_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply rat.A12' correctly

x = rand(dtn.Nt*dtn.Nr,1);

[abserr,relerr] = get_err(rat.A12'*x, rats.apply_A12_star(x));
fprintf('Error between rat.A12 adjoint and apply_A12_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply rat.Tfull(k)' correctly

x = rand(dtn.Nt*dtn.Nr + sum(rats.ratf_lens),1);
k = 2.5;

[abserr,relerr] = get_err(rat.Tfull(k)'*x, rats.apply_Tfull_star(x,k));
fprintf('Error between rat.Tfull(k) adjoint and apply_Tfull_star: %4.2e, %4.2e\n', abserr, relerr);

%% Check that I apply rat.Tfull(k)'*rat.Tfull(k) correctly

[abserr,relerr] = get_err(rat.Tfull(k)'*rat.Tfull(k)*x, rats.apply_Tfull_star_Tfull(x,k));
fprintf('Error between rat.Tfull(k) adjoint rat.Tfull(k) and apply_Tfull_star_Tfull: %4.2e, %4.2e\n', abserr, relerr);

%% Check that solving and solving normal equations give same thing

b = rand(length(rat.A),1);
Tk = rat.Tfull(k);

x        = Tk\b;
x_normal = (Tk'*Tk)\(Tk'*b);
[abserr,relerr] = get_err(x,x_normal);
fprintf('Error between sol to Tk*x = b and normal equations sol: %4.2e, %4.2e\n', abserr, relerr);

%% Check that pcg gives right solution

x_pcg = pcg(Tk'*Tk, Tk'*b,1e-4,1000,[],[],rand(size(b)));
[abserr,relerr] = get_err(x,x_pcg);
fprintf('Error between sol to Tk*x = b and pcg sol of normal eqtns: %4.2e, %4.2e\n', abserr, relerr);

%% Check that pcg with function instead of matrix gives right solution

x_pcg = pcg(@(x) rats.apply_Tfull_star_Tfull(x,k), rats.apply_Tfull_star(b,k),1e-4,1000);
[abserr,relerr] = get_err(x,x_pcg);
fprintf('Error between sol to Tk*x = b and pcg iterative sol of normal eqtns: %4.2e, %4.2e\n', abserr, relerr);



