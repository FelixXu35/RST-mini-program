function matrix = f_matrix(potenetialEnergy, energy)
%f_matrix - This m file contain the form of f matrix
%
% Syntax: matrix = f_matrix(position)
%
% This m file contain the form of f matrix, which is an important part of evolutation 
% This m file gives the f matrix which exactly corresponds to the Schordinger equation
% Input the potential energy at the given position, the energy of the partical
% Output the f matrix
    matrix = [0, 1; 2 * (potenetialEnergy - energy), 0];
end