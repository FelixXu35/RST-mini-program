function yOut = evolve_Schor(yIn, protentialFunc, position, energy, stepLength)
%evolve_Schor - This function is used to evolve the Schordinger equation
%
% Syntax: yOut = evolve_Schor(yIn)
%
% This funciton use predictor and corrector method to evolve the schordinger equation
% 
    predictWave = yIn + f_matix(protentialFunc(position), energy) * yIn * step;
    predictDev = f_matix(protentialFunc(position + stepLength), energy) * predictWave * step;
    yOut = yIn + (predictDev + f_matix(protentialFunc(position), energy) / 2 * yIn * step);
end