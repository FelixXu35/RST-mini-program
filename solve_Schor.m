function waveFunc = solve_Schor(lowerLim, upperLim, lowerBC, upperBC, numOfResult, potenetialFunc)
%solve_Schor - This is the function used to solve the Schordinger equation
%
% Syntax: waveFunc = solve_Schor(lowerLim, upperLim, boundryComdition, potenetialFunc)
%
% Wretten by Xiaotian Xu, 10 March 2020.
% inpit
% output
% This function uses second order central difference scheme to represent the second order derivative of wave function.
% This function uses backward difference scheme to represent the deritive of wave function.
% This function uses Predictor-Corrector method to evolve the differential equation and uses shooting method to solve it.

    %% error detection
    if lowerLim > upperLim
        error('Upper limit should be greater than lower limit.')
    end
    if mod(numOfResult) ~= 0
        error('The number of result should be an integar')
    end

    %% Initialization
    stepLength = 1e-3; 
    numOfSteps = (upperLim - lowerLim) / stepLength + 1;
    position = lowerLim: stepLength: upperLim;
    waveFunc = zeros(numOfResult, numOfSteps);
    newFunc = zeros(numOfResult, numOfSteps);
    oldFunc = zeros(numOfResult, numOfSteps);
    newFunc(1, 1) = lowerBC;
    newFunc(2, 1) = upperBC;
    energy = min(potenetialFunc(position))

    %% find the result
    for resultIndex = 1: numbOfResult
        while signal == false
            oldFunc = newFunc;
            %% Evolve the differential equation
            for stepIndex = 1: (numOfSteps - 1)
                newFunc((2 * resultIndex - 1): (2 * resultIndex), stepIndex + 1) =...
                    evolve_Schor(newFunc((2 * resultIndex - 1): (2 * resultIndex), stepIndex), potenetialFunc, position(stepIndex), energy, stepLength);
            end
            if (newFunc(2 * resultIndex - 1, numOfSteps) - upperBC) * (oldFunc(resultIndex, numOfSteps) - upperBC) < 0
                signal = true;
            end
        end
        signal = false;
    end
    
end