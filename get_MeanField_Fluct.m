function [yMF,zMF,SigmaFluct] = get_MeanField_Fluct(k)
    % Chemical rates
    k1  = k(1);
    km1 = k(2);
    k2  = k(3);
    km2 = k(4);
    k3  = k(5);

    % Mean field
    G =  ( km1*(km2 + k3) + k2*k3 ) / (k2+km2+k3);
    yMF = k1/(k1+G)* (km2+k3)/(k2+km2+k3);
    zMF = k1/(k1+G)* (k2)/(k2+km2+k3);
    
    % Drift and Diffusion for Van Kampen FP
    J= zeros(2,2);
    J(1,1)= -(k1+ km1+ k2);
    J(1,2)= -k1+ km2;
    J(2,1)= k2;
    J(2,2)= -(km2+ k3);
    
    B = zeros(2,2);
    B(1,1)= k1*(1- yMF- zMF)+ (km1+ k2)*yMF+ km2*zMF; 
    B(2,2)= k2*yMF+ (km2+ k3)*zMF;                        
    B(1,2)= -k2*yMF- km2*zMF;   
    B(2,1) = B(1,2);
    
    % Covariance matrix: analytical formula for 2D
    d = det(J); t = trace(J);
    temp = J - eye(2)*t;
    SigmaFluct = ( -temp*B*temp' - d*B )/(2*d*t);

end

