function likelihood = evaluate_likelihood(x, mu, Sigma)
    % Evaluate the likelihood of x-vector based on multivariate Gaussian distribution
    % Inputs:
    %   x: The vector at which to evaluate the likelihood
    %   mu: The mean vector of the Gaussian distribution
    %   Sigma: The covariance matrix of the Gaussian distribution
    % Output:
    %   likelihood: The likelihood of the given x-vector
    
    % Dimensionality of the multivariate Gaussian distribution
    n = length(mu);
    
    % Calculate the exponent term in the Gaussian distribution
    exponent = -0.5 * (x - mu)' * inv(Sigma) * (x - mu);
    
    % Calculate the normalization constant (assuming Sigma is positive definite)
    normalization = (2*pi)^(-0.5*n) * det(Sigma) ^ (-0.5);

    % Verbose
    % disp(num2str(exponent)+ " , "+num2str(normalization)) 
    
    % Calculate the likelihood
    likelihood = normalization * exp(exponent);
end