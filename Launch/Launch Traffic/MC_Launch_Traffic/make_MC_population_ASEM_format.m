function combined_objects = make_MC_population_ASEM_format(constellationFile,constellationSheet, initialPopFile, randomizedRepeatedLaunchFile, scenarioYrs, time0, save_table, folder, recompute_FLC, candidate_table)
    
    % This file combines a given starting population, randomized future
    % traffic, and large constellation model. 

    % Note that the future launch traffic is filtered with a set of rules
    % specific to the starting dataset which should be revised for your
    % files.
    
    % constellationFile is a filepath string for a constellation file in the .xlsx format.
    % constellationSheet is a string for which sheet of the constellation file
    % to use.
    % scenarioYrs is number of years to model future launchs into the future.
    % initialPopFile is starting population
    % randomizedRepeatedLaunchFile is the set of repeating traffic to
    % repeat and randomize over the same amount of time as the file covers
    % scenarioYrs is how many years are in the scenario
    % time0 is the starting datetime
    % save_table is bool to save the output of the function
    % folder is where to save it.
    % candidate is a structure of a candidate constellation.
    % TODO: Switch to taking a table and assemble structure to table outside the method. 

    %% x0 Starting Population
    disp("Loading starting population...")
    initial_x0_data_struct = load(initialPopFile);
    x0_ASEM_format_table = matsats2ASEM(initial_x0_data_struct.initial_data);
    %x0_ASEM_format_table{:,"area"} = pi*(x0_ASEM_format_table.size/2).^2;
    disp("x0: start: " + string(1) + " end: " + string(height(x0_ASEM_format_table)))
    x0_ASEM_format_table = x0_ASEM_format_table(x0_ASEM_format_table.mass < 13055, :); %remove ISS modules.

    %% Repeated Randomized Launch
    disp("Loading Repeated Randomized Launch...")
    rrlf = load(randomizedRepeatedLaunchFile);
    rrlf_ASEM_format_table = matsats2ASEM(rrlf.launch_data);

    % IMPORTANT!
    % These rules are used to clean up the incoming historical launch data.
    % You should change or remove these if not suitable for your
    % application!
    
    rrlf_ASEM_format_table = rrlf_ASEM_format_table(rrlf_ASEM_format_table.mass <= 7000,:); % Remove ISS, dragon, progress, etc.
    rrlf_ASEM_format_table{rrlf_ASEM_format_table.mass <= 10,"stationkeeping"} = 0; % Assume cubesats can't stationkeep
    rrlf_ASEM_format_table = rrlf_ASEM_format_table(~ismember(rrlf_ASEM_format_table.mass, [1250, 750, 700, 148, 260]),:); % Coarse and inaccurate way to remove LLCs.
    %rrlf_ASEM_format_table{:,"area"} = pi*(rrlf_ASEM_format_table.size/2).^2;
    
    % Repeat all data in file each year of sim.
    cycle_launches = {};
    data_duration = years(days(max(rrlf_ASEM_format_table.epoch_start) - min(rrlf_ASEM_format_table.epoch_start)));
    repetitions = ceil(scenarioYrs/data_duration);
    for repitition = 0:repetitions
        rrlf_ASEM_format_table_copy = rrlf_ASEM_format_table;
        start_date = time0 + years(data_duration * repitition);
        end_date = start_date + years(data_duration);
        disp("Regenerating launches for " + string(start_date) + " to " + string(end_date) + ".")
        num_dates = height(rrlf_ASEM_format_table_copy);
        rrlf_ASEM_format_table_copy.epoch_start_datime_matlab = start_date + (end_date - start_date) * rand(num_dates, 1);
        rrlf_ASEM_format_table_copy.lifetime = rrlf_ASEM_format_table_copy.epoch_end - rrlf_ASEM_format_table_copy.epoch_start;
        rrlf_ASEM_format_table_copy.epoch_start = juliandate(rrlf_ASEM_format_table_copy.epoch_start_datime_matlab);
        rrlf_ASEM_format_table_copy.epoch_end = rrlf_ASEM_format_table_copy.epoch_start + rrlf_ASEM_format_table_copy.lifetime;
        rrlf_ASEM_format_table_copy = removevars(rrlf_ASEM_format_table_copy, {'epoch_start_datime_matlab', 'lifetime'});
        cycle_launches{repitition + 1} = rrlf_ASEM_format_table_copy;
    end

    rrlf_ASEM_format_table = vertcat(cycle_launches{:});
    rrlf_ASEM_format_table = rrlf_ASEM_format_table(rrlf_ASEM_format_table.epoch_start >= juliandate(time0) & rrlf_ASEM_format_table.epoch_start <= juliandate(time0 + years(scenarioYrs)), :);
    
    if nargin >= 10 % Candidate network passed.
        disp("Adding candidate satellites...")
        candidate_cycle_launches = {};
        data_duration = years(days(max(candidate_table.epoch_end) - min(candidate_table.epoch_start))); % Assumes traffic replenished going forward
        repetitions = ceil(scenarioYrs/data_duration);
        for repitition = 0:repetitions
            candidate_table_copy = candidate_table;
            start_date = time0 + years(data_duration * repitition);
            end_date = start_date + years(data_duration);
            disp("Regenerating launches for " + string(start_date) + " to " + string(end_date) + ".")
            num_dates = height(candidate_table_copy);
            candidate_table_copy.epoch_start_datime_matlab = start_date + (end_date - start_date) * rand(num_dates, 1);
            candidate_table_copy.lifetime = candidate_table_copy.epoch_end - candidate_table_copy.epoch_start;
            candidate_table_copy.epoch_start = juliandate(candidate_table_copy.epoch_start_datime_matlab);
            candidate_table_copy.epoch_end = candidate_table_copy.epoch_start + candidate_table_copy.lifetime;
            candidate_table_copy = removevars(candidate_table_copy, {'epoch_start_datime_matlab', 'lifetime'});
            candidate_cycle_launches{repitition + 1} = candidate_table_copy;
        end
        rrlf_ASEM_format_table = vertcat(rrlf_ASEM_format_table, candidate_cycle_launches{:});
        rrlf_ASEM_format_table = rrlf_ASEM_format_table(rrlf_ASEM_format_table.epoch_start >= juliandate(time0) & rrlf_ASEM_format_table.epoch_start <= juliandate(time0 + years(scenarioYrs)), :);
    end

    % Assign IDs
    start_id = max(x0_ASEM_format_table.adept_id)+1;
    rrlf_ASEM_format_table.adept_id = (start_id:start_id+height(rrlf_ASEM_format_table)-1).';
    disp("rrlf: start: " + string(start_id) + " end: " + string(start_id+height(rrlf_ASEM_format_table)-1))

    %% Future Launch Constellations
    disp("Loading future constellation population...")
    if recompute_FLC
        fl = mega2matsats(constellationFile,constellationSheet, scenarioYrs, time0);
        fl_ASEM_format_table = matsats2ASEM(fl);
        start_id = max(rrlf_ASEM_format_table.adept_id)+1;
        fl_ASEM_format_table.adept_id = (start_id:start_id+height(fl_ASEM_format_table)-1).';
        disp("fl_LLC: start: " + string(start_id) + " end: " + string(start_id+height(fl_ASEM_format_table)-1))
        fl_ASEM_format_table{fl_ASEM_format_table.stationkeeping >0, "stationkeeping"} = 5 ; %Adding a 5 to indicate slotted orbits;
        save(fullfile(folder,"fl_ASEM_format_table.mat"), "fl_ASEM_format_table")
    else
        load(fullfile(folder,"fl_ASEM_format_table.mat"))
    end

    combined_objects = vertcat(fl_ASEM_format_table, x0_ASEM_format_table, rrlf_ASEM_format_table);
    combined_objects = sortrows(combined_objects, "adept_id");

    if save_table
        disp("Saving file.")
        [~, pop_name, ~] = fileparts(initialPopFile);
        [~, repeat_name, ~] = fileparts(randomizedRepeatedLaunchFile);
        [~, cons_name, ~] = fileparts(constellationFile);
        cons_name = cons_name + "_" + constellationSheet;
        
        %fileName = matlab.desktop.editor.getActiveFilename;
        %folder = fileparts(which(fileName)); 
        
        name = pop_name + "_" + repeat_name + "_" + cons_name;
        writetable(combined_objects, fullfile(folder, name + ".csv"))
    end
% 
% % Testing
% constellationFile = "megaconstellationLaunches.xlsx";
% constellationSheet = "forexport2";
% initialPopFile = "initial_population_20230101.mat";
% randomizedRepeatedLaunchFile = "launch_every_year_randomtime.mat";
% scenarioYrs = 200;
% save = true;
% objects = make_MC_population_ASEM_format(constellationFile,constellationSheet, initialPopFile, randomizedRepeatedLaunchFile, scenarioYrs, save);