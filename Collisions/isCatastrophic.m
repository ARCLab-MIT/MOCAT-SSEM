function [overThreshold] = isCatastrophic(mass1, mass2, vels)
    % Determines if collision is catastrophic or non-catastrophic based on
    %   40 j/g threshold from Jonhson et al. 2001.
    %   mass1 is mass of first object in kg
    %   mass2 is the mass of the second object in kg
    %   vels is the array of relative orbital velocities to use to evaluate
    %   the collisions in each shell
    % Returns shell-wise list of bools (true if energy>=40 J/g, else false)
    
    if mass1 <= mass2
        smaller_mass = mass1;
    else
        smaller_mass = mass2;
    end
    smaller_mass_g = smaller_mass * (1000/1); % kg -> grams
    energy = 1/2 * smaller_mass * vels.^2;
    overThreshold = energy/smaller_mass_g > 40; %