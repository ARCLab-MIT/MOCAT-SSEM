
function T = synthLambdaData(start_date, end_date, rows, save_path)
    % start_date, end_date like start_date='2022-08-01'
    % rows is int
    % save_path is save location with name and format. Default
    % rand_sats.csv in launch traffic

    %% Parameters
    
    if ~exist('rows','var') || isempty('rows')
      rows=500;
    end
    
    % Dates
    if ~exist('start_date','var') || isempty('start_date')
      start_date='2022-08-01';
    end
    
    if ~exist('end_date','var') || isempty('end_date')
        end_date='2072-12-31';
    end
    
    if ~exist('save_path','var') || isempty('save_path')
      save_path = fullfile(pwd, 'Launch Traffic', 'rand_sats.csv');
    end
    
    % Masses
    mu_mass =500; %kg
    sigma_mass = 200; %kg
    
    %Inclination
    inc_min = 1; % deg
    inc_max = 90.9; % deg
    
    %Sats
    min_sats = 1;
    max_sats = 2000;

    %Sats
    min_alt = 200;
    max_alt = 1200;
    
    % Dates
    %https://www.mathworks.com/matlabcentral/answers/81241-how-to-generate-random-time-format-data-in-matlab
    SECONDS_PER_DAY = 24*60*60;
    startDateNum=datenum(start_date, 'yyyy-mm-dd');
    endDateNum=datenum(end_date, 'yyyy-mm-dd'); 
    dayRange = endDateNum - startDateNum;
    secondsRange = SECONDS_PER_DAY*dayRange;
    randomNumberOfSeconds = rand_in_range(0, secondsRange, [rows, 1]);
    randomDatenums = startDateNum + randomNumberOfSeconds/SECONDS_PER_DAY;
    randomDatenums = sort(randomDatenums);
    
    %Inclination
    randIs = rand_in_range(inc_min, inc_max, [rows, 1]);
    
    %Masses
    
    randMs = rand_in_range_normal(mu_mass, sigma_mass, [rows, 1]);
    
    % Sats
    randSats =randi([min_sats max_sats],[rows, 1]);
    
    % Alts
    alts = rand_in_range(min_alt, max_alt, [rows, 1]);

    data = [randomDatenums, randIs, randMs, randSats, alts];
    
    T = array2table(data, 'VariableNames', ["Launch_Date" "Inclination_deg" "Mass_kg" "Sats" "Alt_km"]);
    T.Launch_Date = datestr(T.Launch_Date);
    T.Type(:) = "Su";

    writetable(T,save_path);

end

function r = rand_in_range(lower, upper, shape)
    r = (upper-lower).*rand(shape) + lower;
end

function r = rand_in_range_normal(mu, sigma, shape)
    R = chol(sigma);    
    r = repmat(mu,shape) + randn(shape)*R;
    r(r<=0.1) = 0.1; %avoid negative mass.
end

