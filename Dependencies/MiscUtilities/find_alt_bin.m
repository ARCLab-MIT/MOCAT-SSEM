function shell_index = find_alt_bin(altitude, scen_properties)
    % Takes an altitude and returns the relevant altitude bin
    lower = scen_properties.R02(1:end-1).';
    upper = scen_properties.R02(2:end).';
    shell_logic = horzcat(lower < altitude, altitude <= upper);
    shell_logic_sum = sum(shell_logic, 2);
    if any(shell_logic_sum==2)
         shell_index = find(shell_logic_sum == 2);
    else
        shell_index = NaN;
    end
end
