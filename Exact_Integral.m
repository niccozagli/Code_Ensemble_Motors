function pdf = Exact_Integral(x,y,yMF,zMF,SigmaFluct,N,f0,f_grid)
    nyMF = N * yMF; nzMF = N *zMF;
    Sigma = N*SigmaFluct;
    [q1,q2] = meshgrid(x,y);
    % Evaluate invariant distribution population
    q1c = q1 - nyMF; q2c = q2 - nzMF;

    A = inv(Sigma);
    G = 1/(2*pi*sqrt(det(Sigma)))* exp(-1/2* ( A(1,1)*q1c.^2 + 2*A(1,2)*q1c.*q2c + A(2,2)*q2c.^2));
    
    % Evaluate exact distribution of force
    mu_Force = 11/20*f0*q2; s2_Force = (1/3*q1 + 27/400*q2)*f0^2;
    gForce = @(f) 1 ./sqrt(2*pi*s2_Force) .* exp(-1/2 * (f- mu_Force).^2 ./s2_Force );
    I = @(f) trapz(y,trapz(x,G .* gForce(f),2));
    
    pdf = zeros(length(f_grid),1);
    for i = 1: length(f_grid)
        f = f_grid(i);
        pdf(i) = I(f);
    end

end

