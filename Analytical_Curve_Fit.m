function Iformula4 = Analytical_Curve_Fit(k,f0,N)

    % Get the mean field quantities and Fluctuations in Van Kampe approx
    [yMF,zMF,SigmaFluct] = get_MeanField_Fluct(k);
    % Approximate formula to 4th order
    [~,~,Iformula4] = Approximate_Formula(yMF,zMF,SigmaFluct,N,f0); 

end

