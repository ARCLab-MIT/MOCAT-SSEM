%file_path = ".\Launch Traffic\rand_sats.csv";
%bin_Lambda_Data(file_path, scen_properties)

%%

function T2 = bin_Lambda_Data(file_path, scen_properties)
    T = readtable(file_path);
    speciesName_list = unique(T.Type);
    T.AltBin(:) = 0;

    variable_names_types = [["Launch_Date", "datetime"]; ...
			    ["Inclination_deg", "double"]; ...
			    ["Mass_kg", "double"]; ...
			    ["Sats", "double"]; ...
			    ["Alt_km", "double"]; ...
			    ["Type", "string"];...
                ["AltBin", "double"]; ...
                ["AltBinAlt", "double"]; ...
                ["MassBin", "string"]];

    T2 = table('Size',[0,size(variable_names_types,1)],... 
	            'VariableNames', variable_names_types(:,1),...
	            'VariableTypes', variable_names_types(:,2));

    for k = 1:scen_properties.N_shell-1 % For each shell
        % Find rows in bin
        bin_min = scen_properties.R02(k);
        bin_max = scen_properties.R02(k+1);
        %disp(["Min Alt", bin_min, "Max Alt", bin_max])
        rows = T(T.Alt_km >= bin_min & T.Alt_km < bin_max, :);
        rows.AltBin(:) = k;
        rows.AltBinAlt(:) = scen_properties.HMid(k);
        for s = 1:length(speciesName_list) % For each species
            speciesName = string(speciesName_list(s));
            spec_rows = rows( strcmp(rows.Type,speciesName), :);
            spec_rows.MassBin(:) = "None";
            %disp(["    Species Name", speciesName])
            %spec_rows;
            for mass_i = 1:length(scen_properties.species_cell.(speciesName)) % For each mass bin
                mass_lb = scen_properties.species_cell.(speciesName)(mass_i).species_properties.mass_lb;
                mass_ub = scen_properties.species_cell.(speciesName)(mass_i).species_properties.mass_ub;
                %disp(["Mass lb", mass_lb, "Mass ub", mass_ub])
                mass_rows = spec_rows(spec_rows.Mass_kg >= mass_lb ...
                                      & spec_rows.Mass_kg < mass_ub, :);
                %disp(["Min Alt", bin_min, "Max Alt", bin_max])
                if height(mass_rows) > 0
                    %disp(mass_rows)
                    mass_rows.MassBin(:) = scen_properties.species_cell.(speciesName)(mass_i).species_properties.sym_name;
                    T2 = [T2; mass_rows];
                end
            end % mass loop
        end % species loop
    end % alt loop
    
    T2.Launch_T = years(T2.Launch_Date - scen_properties.start_date);
    %lambda_table = T2(strcmp(T2.MassBin, 'Su_250kg'),:);
    %lambda_fit2(lambda_table, scen_properties)
    disp("Done")

end