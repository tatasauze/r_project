function test_conversion_matlab()
% Test script to run MATLAB implementation of thd2arma

% Add the E4withSubspaces directory to the path
addpath('E4withSubspaces');

% Initialize E4 environment
e4init();

% Load test data
load('test_data.mat');

fprintf('Running MATLAB implementation...\n');

% Run thd2arma in MATLAB
[F, A, V, G] = thd2arma(theta, din);

fprintf('MATLAB implementation completed.\n');

% Save results for comparison
save('matlab_results.mat', 'F', 'A', 'V', 'G');

fprintf('Results saved to matlab_results.mat\n');

% Print results
fprintf('MATLAB Results:\n');
fprintf('F Matrix:\n');
disp(F);
fprintf('A Matrix:\n');
disp(A);
fprintf('G Matrix:\n');
disp(G);
fprintf('V Matrix:\n');
disp(V);

% Compare with R results if available
if exist('r_results.mat', 'file')
    load('r_results.mat');
    fprintf('\nComparison with R results:\n');
    
    if exist('r_F', 'var') && exist('r_A', 'var') && exist('r_G', 'var') && exist('r_V', 'var')
        F_diff = max(max(abs(F - r_F)));
        A_diff = max(max(abs(A - r_A)));
        G_diff = max(max(abs(G - r_G)));
        V_diff = max(max(abs(V - r_V)));
        
        fprintf('Max difference in F matrix: %e\n', F_diff);
        fprintf('Max difference in A matrix: %e\n', A_diff);
        fprintf('Max difference in G matrix: %e\n', G_diff);
        fprintf('Max difference in V matrix: %e\n', V_diff);
        
        tolerance = 1e-10;
        if F_diff < tolerance && A_diff < tolerance && G_diff < tolerance && V_diff < tolerance
            fprintf('SUCCESS: MATLAB and R implementations produce equivalent results!\n');
        else
            fprintf('FAILURE: MATLAB and R implementations produce different results!\n');
        end
    end
end

end