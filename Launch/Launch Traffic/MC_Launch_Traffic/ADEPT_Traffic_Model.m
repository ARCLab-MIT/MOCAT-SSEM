function [x0,FLM_steps] = ADEPT_Traffic_Model(file_path, scen_properties)
% Takes a set of objects in the modified ADEPT format and a scen_properties
% object and bins those objects into the relevant species based on
% altitude, time, and mass to create initial population and future launch.

    T = readtable(file_path, 'TextType','string');
    func = @(x) datetime(x,'convertfrom','juliandate');
    epoch_start_datime = varfun(func,T,'InputVariables','epoch_start');
    T.epoch_start_datime_matlab = epoch_start_datime.(1);
    
    if ~ismember('obj_class', T.Properties.VariableNames)
        T = define_obj_classes(T);
        T = T(T.obj_class ~= "Unknown",:);
    end

    % Altitude
    T.apogee = T.sma .* (1 + T.ecc);
    T.perigee = T.sma .* (1 - T.ecc);
    T.alt = (T.apogee + T.perigee)/2 - scen_properties.re;
    
    % Species Type
    species_dict = dictionary("Non-station-keeping Satellite", ...
                              "Sns", "Rocket Body", "B", ...
                              "Station-keeping Satellite", "Su", ...
                              "Coordinated Satellite", "S",...
                              "Debris", "N", "Candidate Satellite", "C");
    T.species_class = species_dict(T.obj_class);
    tabulate(T.obj_class)
    obj_classes = unique(T.obj_class);
    
    T_new = T(T.adept_id == -100, :); % empty table
    T_new.species =  strings(height(T_new),1);


    cluster_dict = dictionary("Non-station-keeping Satellite", 1, ...
                              "Station-keeping Satellite", 2, ...
                              "Coordinated Satellite", 3,...
                              "Debris", 2, "Rocket Body", 1, "Candidate Satellite", 1);
    %spec_prop_analysis(T, species_dict, cluster_dict, 1e5, "D:\Dropbox (MIT)\ARCLab Miles\Our Papers and Presentations\AMOS23\Species Clusters")

    % Mass Bin
    disp("Assigning objects to corresponding mass classes.")
    for ii = 1:length(obj_classes)
        obj_class = obj_classes(ii);
        species_class = species_dict(obj_class);
        try
            species_cell = getfield(scen_properties.species_cell,species_class);
        catch
            warning(species_class + " not found in MOCAT model. Omitting from mass class determination")
            continue
        end

        disp([obj_class, species_class])
        T2 = T(T.obj_class == obj_class,:);
        massbin_func = @(row) find_mass_bin(row, scen_properties, species_cell);
        T2species = rowfun(massbin_func,T2, 'InputVariables', "mass");
        T2.species = T2species{:,:};
        T_new = vertcat(T_new, T2);
    end

    % Add slot altitudes.
    disp("Assigning objects to corresponding altitude bins.")
    altbin_func = @(row) find_alt_bin(row, scen_properties);
    Taltbins = rowfun(altbin_func,T_new, 'InputVariables', "alt");
    T_new.alt_bin = Taltbins{:,:};
    T = T_new;
    
    % Remove any species from T that aren't in scen_properties.
    species_cell_names = string(fieldnames(scen_properties.species_cell));
    T = T(ismember(T.species_class, species_cell_names),:);

    % X0
    % We now need to pivot to a table with altitude bins on the row axis and
    % species on the cols axis. Then check if any species have no launches,
    % add them back and reorder to match the scen_properties.
    species_names = [];
    for ii = 1:length(scen_properties.species)
        species_names = [species_names scen_properties.species(ii).species_properties.sym_name];
    end

    x0 = T(T.epoch_start_datime_matlab < scen_properties.start_date,:);
    if ~isempty(x0)
        x0 = pivot(x0,Columns="species",Rows="alt_bin",Method="sum",DataVariable="weight");
        x0 = x0(~isnan(x0.alt_bin),:);
        x0.Properties.RowNames = string(x0{:,"alt_bin"}); % Make alt_bin the rowname
        x0 = removevars(x0,{'alt_bin'});
    end

    emptyx0 = array2table(zeros(scen_properties.N_shell,length(species_names)), RowNames = string(1:scen_properties.N_shell), VariableNames=cellstr(species_names));
    
    rowLabels_base = categorical(emptyx0.Properties.RowNames);
    colLabels_base = categorical(emptyx0.Properties.VariableNames);
    
    rowLabels_step = categorical(x0.Properties.RowNames);
    colLabels_step = categorical(x0.Properties.VariableNames);
    
    % Find matching rows and columns
    [commonRows, baseRowIndices] = ismember( rowLabels_step, rowLabels_base);
    [commonCols, colIndices] = ismember( colLabels_step, colLabels_base);
    
    if ~isempty(x0)
        % Update matching cells in FLM_step_base
        emptyx0(baseRowIndices,  colIndices(colIndices>0)) = x0(commonRows, commonCols(commonCols>0));
        
        % Sort cols to match order of species_list
        x0 = movevars(emptyx0, species_names);
    end
    
    %% Future Launch Model
    
    disp("Assembling future launch information. This may take a while...")
    FLM = T(T.epoch_start_datime_matlab >= scen_properties.start_date, :);
    time_steps = scen_properties.start_date + years(scen_properties.scen_times);
    FLM_steps = cell(length(time_steps)-1, 1);
    for ii = 1:length(time_steps)-1
        if mod(ii,10) == 0
            disp("Start Time: " + string(time_steps(ii)) + " Stop Time: " + string(time_steps(ii+1)))
        end
        FLM_step = T(T.epoch_start_datime_matlab > time_steps(ii)  & T.epoch_start_datime_matlab <= time_steps(ii+1), :);
        
        % If there is data, pivot table it, otherwise make a blank table.
        if ~isempty(FLM_step)
            FLM_step = pivot(FLM_step,Columns="species",Rows="alt_bin",Method="sum",DataVariable="weight");
            FLM_step = FLM_step(~isnan(FLM_step.alt_bin),:); % Remove NaN reow.
            FLM_step.Properties.RowNames = string(FLM_step{:,"alt_bin"});
            FLM_step = removevars(FLM_step,{'alt_bin'});
        else
            FLM_step = array2table(zeros(scen_properties.N_shell,length(species_names)), RowNames = string(1:scen_properties.N_shell), VariableNames=cellstr(species_names));
        end

        FLM_step_base = array2table(zeros(scen_properties.N_shell,length(species_names)), RowNames = string(1:scen_properties.N_shell), VariableNames=cellstr(species_names));

        rowLabels_base = categorical(FLM_step_base.Properties.RowNames);
        colLabels_base = categorical(FLM_step_base.Properties.VariableNames);
        
        rowLabels_step = categorical(FLM_step.Properties.RowNames);
        colLabels_step = categorical(FLM_step.Properties.VariableNames);
        
        % Find matching rows and columns
        [commonRows, baseRowIndices] = ismember( rowLabels_step, rowLabels_base);
        [commonCols, colIndices] = ismember( colLabels_step, colLabels_base);
        
        % Update matching cells in FLM_step_base
        FLM_step_base(baseRowIndices,  colIndices(colIndices>0)) = FLM_step(commonRows, commonCols(commonCols>0));

        %disp(FLM_step_base)
        %disp(FLM_step)
        FLM_steps(ii) = {FLM_step_base};
    end


%     for ii = 1:length(species_list)
%         species = species_list(ii).species_properties.sym_name;
%         disp("Species" + species)
%         T_spec = T(T.species == species,:);
%         [G, ID] = findgroups(T_spec.alt_bin);
%         counts = splitapply(@sum,T_spec.weight,G);
%         x0_agg(:,ii) = counts;
%         pivot(T_spec,Columns=alt_bin,Rows=alt_bin)
%     end
% 
%     altitudes = T.alt_bin;
%     names = T.species;
%     unique_ids = uint64((1:length(names)) + (max(altitudes)+1) * (uint64(altitudes) - min(altitudes)));
%     result = accumarray(unique_ids', weight, [], @sum);
%     unique_altitudes = unique(altitudes);
%     unique_names = unique(names);
% 
%     unique_names = unique(names);
%     unique_altitudes = unique(altitudes);
%     
%     unique_ids = zeros(size(names));
%     for i = 1:length(unique_names)
%         for j = 1:length(unique_altitudes)
%             mask = strcmp(names, unique_names{i}) & altitudes == unique_altitudes(j);
%             unique_ids(mask) = i + (j - 1) * length(unique_names);
%         end
%     end
% 
% 
%     final_matrix = zeros(length(unique_altitudes), length(unique_names));
% 
%     for i = 1:length(unique_names)
%         indices = find(strcmp(names, unique_names{i}));
%         for j = 1:length(indices)
%             final_matrix(find(unique_altitudes == altitudes(indices(j))), i) = result(indices(j));
%         end
%     end
% 
%     groupfilter(T_spec,"alt_bin")

    %TODO Aggregate into bins per timestep.
    % for bin  range o to <= 1
    %G = findgroups(T{:, "obj_class"})
    %T_split = splitapply( @(varargin) varargin, T , G);

end
% %%
% % Binning
% count_list = {};
% pop_list = {};
% speciesName_list = {};
% if isfield(scen_properties.species_cell, 'S')
%     s_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.S)); 
%     count_list{length(count_list)+1} = s_count;
%     %pop_list{length(pop_list)+1} = popMC_S;
%     speciesName_list{length(speciesName_list)+1} = "S";
% end
% if isfield(scen_properties.species_cell, 'D')
%     d_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.D));
%     count_list{length(count_list)+1} = d_count;
%     pop_list{length(pop_list)+1}  = popMC_D;
%     speciesName_list{length(speciesName_list)+1}= "D";
% end
% if isfield(scen_properties.species_cell, 'N')
%     n_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.N)); 
%     count_list{length(count_list)+1} = n_count;
%     %pop_list{length(pop_list)+1}  = popMC_N;
%     speciesName_list{length(speciesName_list)+1} = "N";
% end
% if isfield(scen_properties.species_cell, 'Su')
%     su_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.Su)); 
%     count_list{length(count_list)+1} = su_count;
%     %pop_list{length(pop_list)+1}  = popMC_Su;
%     speciesName_list{length(speciesName_list)+1} = "Su";
% end
% if isfield(scen_properties.species_cell, 'B')
%     b_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.B));
%     count_list{length(count_list)+1} = b_count;
%     %pop_list{length(pop_list)+1}  = popMC_B;
%     speciesName_list{length(speciesName_list)+1} = "B";
% end
% if isfield(scen_properties.species_cell, 'U')
%     u_count = zeros(scen_properties.N_shell, length(scen_properties.species_cell.U)); 
%     count_list{length(count_list)+1} = u_count;
%     %pop_list{length(pop_list)+1}  = popMC_U;
%     speciesName_list{length(speciesName_list)+1} = "U";
% end
% 
% for k = 1:scen_properties.N_shell-1
%     for i = 1:length(count_list)
%         %ind = count_list{i};
%         pop = pop_list{i};
%         speciesName = speciesName_list{i};
%         for ind = 1:length(pop)
%             a_t_shell = pop{ind}.a * scen_properties.re - scen_properties.re;
%             if a_t_shell>=scen_properties.R02(k) && a_t_shell<scen_properties.R02(k+1)
%                 species_cell = getfield(scen_properties.species_cell,speciesName);
%                 for mass_i = 1:length(species_cell)
%                     if (pop{ind}.mass > species_cell(1,mass_i).species_properties.mass_lb && ...
%                         pop{ind}.mass < species_cell(1,mass_i).species_properties.mass_ub )
% %                         disp(["Obj Mass " + num2str(pop{ind}.mass), ...
% %                               " Species mass " + num2str(species_cell(mass_i).species_properties.mass), ...
% %                               " mass lb " + num2str(species_cell(mass_i).species_properties.mass_lb), ...
% %                               " mass ub " + num2str(species_cell(mass_i).species_properties.mass_ub)])
% %                         disp("    Species Found")
%                         break %Found right mass bin
%                     else
%                        continue % Keep looking
%                     end
%                 end
%                 count_list{i}(k, mass_i) = count_list{i}(k, mass_i)+1; % count the number of active satellites in each shell
%             end % assign to right shell and species
%         end % sat loop
%     end % species loop
% end % alt loop

function species_name = find_mass_bin(mass, scen_properties, species_cell)
    for mass_i = 1:length(species_cell)
        if (mass >= species_cell(1,mass_i).species_properties.mass_lb && ...
            mass < species_cell(1,mass_i).species_properties.mass_ub )
%             disp(["Obj Mass " + num2str(row.mass), ...
%                   " Species mass " + num2str(species_cell(mass_i).species_properties.mass), ...
%                   " mass lb " + num2str(species_cell(mass_i).species_properties.mass_lb), ...
%                   " mass ub " + num2str(species_cell(mass_i).species_properties.mass_ub)])
%             disp("    Species Found")
            break %Found right mass bin
        else
           continue % Keep looking
        end
    end
    species_name = species_cell(1,mass_i).species_properties.sym_name;
end

function shell_index = find_alt_bin(altitude, scen_properties)
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

function spec_prop_analysis(T, species_dict, cluster_dict, max_size, folder)
    %clusters = [2, 1, 4, 3];
    obj_classes = unique(T.obj_class);
    for ii = 1:length(obj_classes)
        %figure()
        obj_class = obj_classes(ii);
        species_class = species_dict(obj_class);
        T2 = T(T.obj_class == obj_class,:);
        
        if nargin < 4
            max_size = 1e6;
        end
        if height(T2) > max_size
            T2 = datasample(T2, max_size, 'Replace', false);
        end

        T2.lifetime = T2.epoch_end - T2.epoch_start;
        % histogram(T2.mass, 'Normalization','probability', 'DisplayStyle','stairs')
        % xlabel('Mass [kg]')
        % ylabel('Count')
        % title(obj_class + " Mass")

        % figure()
        % A = sort(T2.mass);
        % [B,TFrm,TFoutlier,L,U,C] = rmoutliers(A, "mean");
        % plot(A)
        % hold on
        % plot(find(~TFrm),B,"o-")
        % yline([L U C],":",["Lower Threshold","Upper Threshold","Center Value"])
        % legend("Original Data","Cleaned Data")
        % title(obj_class + " mass outlier rejection")
        % 
        % figure()
        % obj_class = obj_classes(ii);
        % species_class = species_dict(obj_class);
        % T2 = T(T.obj_class == obj_class,:);
        % histogram(T2.area, 'Normalization','probability', 'DisplayStyle','stairs')
        % xlabel('Area [m^2]')
        % ylabel('Count')
        % title(obj_class + " Area")
        % 
        % figure()
        % A = sort(T2.area);
        % [B,TFrm,TFoutlier,L,U,C] = rmoutliers(A, "mean");
        % plot(A)
        % hold on
        % plot(find(~TFrm),B,"o-")
        % yline([L U C],":",["Lower Threshold","Upper Threshold","Center Value"])
        % legend("Original Data","Cleaned Data")
        % title(obj_class + " area outlier rejection")

        % K Means Analysis
        rng(1); 
        sum_means = [];
        sum_meds = [];
        cluster_count_range = 1:6;
        for cluster_count = cluster_count_range
            cluster_count
            [idx,C, sumd_mean] = kmeans(T2{:, ["area", "mass", "size", "lifetime"]},cluster_count);
            [idx,C2, sumd_med] = kmedoids(T2{:, ["area", "mass", "size", "lifetime"]},cluster_count);
            sum_means(cluster_count) = sum(sumd_mean);
            sum_meds(cluster_count) = sum(sumd_med);
        end
    
        figure()
        hold on
        plot(cluster_count_range, sum_means, "DisplayName", "Means")
        plot(cluster_count_range, sum_meds, "DisplayName", "Medians")
        xlabel('Clusters')
        ylabel('Within-Cluster Sum of Squares')
        title(obj_class + " Clustering Count")
        legend()
        hold off
        saveas(gcf,fullfile(folder, obj_class + ".png"))
        saveas(gcf,fullfile(folder, obj_class + ".fig"))

        figure()
        hold on
        [idx,C, sumd] = kmeans(T2{:, ["area", "mass", "size", "lifetime"]},cluster_dict(obj_class));
        [idx,C2, sumd] = kmedoids(T2{:, ["area", "mass", "size", "lifetime"]},cluster_dict(obj_class));
        scatter(T2.area, T2.mass, ".k", "DisplayName", "FLM Data")
        scatter(C(:,1), C(:,2), "ro", "DisplayName", "K-Means Centroids")
        scatter(C2(:,1), C2(:,2), "bo", "DisplayName","M-Medians Centroids")
        title(obj_class + " Clustering" + " (clusters = " + string(cluster_dict(obj_class)) + ")" )
        legend()
        xlabel('Area [m^2]')
        ylabel('Mass [kg]')
        hold off
        saveas(gcf,fullfile(folder, obj_class + "_Clusters_" + string(cluster_dict(obj_class)) + ".png"))
        saveas(gcf,fullfile(folder, obj_class + "_Clusters_" + string(cluster_dict(obj_class)) + ".fig"))

        disp("Done with " + obj_class + ".")
        disp("Mass med: " + median(T2.mass) + " mean: " + mean(T2.mass))
        disp("Area med: " + median(T2.area) + " mean: " + mean(T2.area))

        C_table = array2table(C, "VariableNames", {'area [m^2]', 'mass [kg]', 'size [m]', 'lifetime [days]'});
        C2_table = array2table(C2, "VariableNames", {'area [m^2]', 'mass [kg]', 'size [m]', 'lifetime [days]'});
        writetable(C_table, fullfile(folder, obj_class + "_kmeans_" + string(cluster_dict(obj_class)) + "_clusters.csv"))
        writetable(C2_table, fullfile(folder, obj_class + "_kmedians_" + string(cluster_dict(obj_class)) + "_clusters.csv"))
        
    end
end

function T = define_obj_classes(T)
    % Classifies the classes within the table T and adds them as a new 
    % column labeled "obj_type" or overwrites it if the column exists.
    
    obj_type_dict = containers.Map([1, 2, 3], {'Rocket Body', 'Satellite', 'Debris'});

    T.obj_class = repmat("Unknown", height(T), 1);

    % Rocket Bodies
    T.obj_class(T.obj_type == 1) = "Rocket Body";

    % Satellites 
    T.obj_class((T.obj_type == 2) & (T.stationkeeping ~= 0 & T.stationkeeping <5)) = "Station-keeping Satellite";
    T.obj_class((T.obj_type == 2) & (T.stationkeeping == 0)) = "Non-station-keeping Satellite";
    T.obj_class((T.obj_type == 2) & (T.stationkeeping==5)) = "Coordinated Satellite";
    T.obj_class((T.obj_type == 2) & (T.stationkeeping==6)) = "Candidate Satellite";

    % Debris
    T.obj_class(T.obj_type == 3 | T.obj_type == 4) = "Debris";

    unclassed_rows = sum(T.obj_class == "Unknown");
    if unclassed_rows > 0
        fprintf('\t%d Unclassified rows remain.\n', unclassed_rows);
        % warning('%d Unclassified rows remain.', unclassed_rows);
    end
end
