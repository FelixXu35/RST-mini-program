function waveFunc = solve_Schor_single(lowerLim, upperLim, lowerBC, upperBC, potenetialFunc)
%solve_Schor_single - This function can find the wave function at the ground state
%
% Syntax: waveFunc = solve_Schor_single(lowerLim, upperLim, lowerBC, upperBC, )
%
% Long description

    %% error detection
    if lowerLim > upperLim
        error('Upper limit should be greater than lower limit.')
    end

    %% Initialization
    stepLength = 1e-3; % choose 1000 pts in one unit of length
    numOfSteps = (upperLim - lowerLim) / stepLength + 1;
    position = lowerLim: stepLength: upperLim;
    waveFunc = zeros(2, numOfSteps);
    newFunc = zeros(2, numOfSteps);
    oldFunc = zeros(2, numOfSteps);
    energy = min(potenetialFunc(position));
    newFunc(1, 1) = lowerBC;
    newFunc(1, numOfSteps) = upperBC;
    
    %% find the result
    while signal == false
        oldFunc = newFunc;
        %% Evolve the differential equation
        for stepIndex = 1: (numOfSteps - 1)
            newFunc((2 * resultIndex - 1): (2 * resultIndex), stepIndex + 1) =...
                evolve_Schor(newFunc((2 * resultIndex - 1): (2 * resultIndex), stepIndex), potenetialFunc, position(stepIndex), energy, stepLength);
        end
        if (newFunc(2 * resultIndex - 1, numOfSteps) - upperBC) * (oldFunc(resultIndex, numOfSteps) - upperBC) < 0
            signal = true;
        else
            energy = energy + 1; % if the result is not in this interval, search next interval
            newFunc()
        end
    end    
end