function test_all(verbose)

if nargin == 0, verbose = 0; end

if verbose
    fprintf('\nAbout to run all tests with errors and figures shown...\n\n');
else
    fprintf('\nAbout to run all tests with errors and figures shown\n');
    fprintf('only if errors larger than expected...\n\n');
end

disp('***Running test_laplacianEigs2D...'); 
test_laplacianEigs2D(verbose);
disp('***Running test_scattResComp2d_plotting...'); 
test_scattResComp2d_plotting(verbose);
disp('***Running test_consistency...');
test_consistency(verbose);
disp('***Running test_sparse_versions...');
test_sparse_versions(verbose);
disp('***Running test_compute_truescatt...');
test_compute_truescatt(verbose);
disp('***Running test_axisymm2D...');
test_axisymm2D(verbose);
disp('***Running test_bumpSum2D...');
test_bumpSum2D(verbose);
disp('***Running test_singleBump2D...');
test_singleBump2D(verbose);
disp('***Running test_sqrBump2D...');
test_sqrBump2D(verbose);
disp('***Running test_ratApprox...');
test_ratApprox(verbose);
disp('***Running test_preconditioners...');
test_preconditioners(verbose);
disp('***Running test_chebBasis...');
test_chebBasis(verbose);
disp('***Running test_axisymm_objects...');
test_axisymm_objects(verbose);

