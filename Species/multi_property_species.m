classdef multi_property_species < handle
    %   Allows you to provide a struct formatted to construct a 
    %   species_properties_class instance. 
    % 
    %   This class uses that struct to create a multi_property_species 
    %   class with the species_list property with a set of species with 
    %   different properties.
    %
    %   If multiple masses are provided, but other values are provided
    %   single values, then radius, A, amr, and beta will be scaled based 
    %   on a spherical assumption. Trackability will be set based on scaled
    %   radius relative to theTRACKABLE_RADIUS_THRESHOLD constant.

    properties
        species_list
    end

    methods
        function obj = multi_property_species(launch_func, pmd_func, ...
                drag_func, species_properties, scen_properties)
            %UNTITLED Construct an instance of this class
            % launch_func: a function that specifies launch rate. Launch 
            %   functions takes a parameter t (float), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the launch quantity at the time for 
            %   this species as a function of time.
            %   TODO: figure out of this will be symbolic of quantity.
            % pmd_func: a function that specifies the PMD rate. The PMD 
            %   functions takes a parameter t (float), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the PMD rate for this species as a
            %   function of time.
            %   TODO: figure out of this will be symbolic of quantity.
            % drag_func: a function that specifies the reduction in 
            %   altitude per unit time. The drag functions takes a 
            %   parameter t (float), h (altitude in km), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the PMD rate for this species as a
            %   function of time.
            
            TRACKABLE_RADIUS_THRESHOLD = 0.05; %(m)
            
            % Char arrays have length based on string length which confuses
            % the logic to parse fields for multi-property species classes.
            % Instead, any char arrays, e.g. 'x' are converted to string
            % arrays e.g. "x" which respond to length() as expected by the
            % below code.
            
            if length(species_properties.mass) == 1
                error("Use species, not multi_property_species when " + ...
                    "creating species with a single mass value")
            end
            
            if contains(species_properties.sym_name, "_") 
                error(" _ is a reserved character and cannot be used in species names" )
            end

            species_properties_fieldnames = fieldnames(species_properties);
            for field_index=1:length(fieldnames(species_properties))
                     fieldname = species_properties_fieldnames(field_index);
                     if ischar(species_properties.(fieldname{:}))
                        species_properties.(fieldname{:}) = convertCharsToStrings(species_properties.(fieldname{:}));
                     end
            end
            
            num_species = length(species_properties.mass);
            for i = 1:num_species
                % Make copy of provided struct to adjust for specific
                % species.
                species_properties_copy = species_properties;
                species_properties_copy.mass = species_properties.mass(i);
                species_properties_copy.sym_name = species_properties.sym_name + "_" + num2str(species_properties.mass(i)) + "kg";
                
                our_fieldnames = fieldnames(species_properties);
                for field_index=1:length(our_fieldnames)
                    fieldname = our_fieldnames(field_index);
                    % If there is a value for each mass, use the provided
                    % values
                    if fieldname == "sym_name" %Don't want to overwrite the first name
                        continue
                    end
                    if length(species_properties.(fieldname{:})) == num_species || i == 1 %First species get parent properties
                        species_properties_copy.(fieldname{:}) = species_properties.(fieldname{:})(i);

                    % If a set of per-species properties is not explictly 
                    % provided, use the same value for all fields except
                    % radius, A. For those assume sphere
                    % scaling relationships.
                    elseif length(species_properties.(fieldname{:})) == 1 
                         if strcmp(fieldname,'radius') || strcmp(fieldname,'A')
                             org_volume = (4/3)*pi*species_properties.radius^3;
                             org_density = species_properties.mass(1)/org_volume;
                         end
                         if strcmp(fieldname,'radius')
                            species_properties_copy.radius = ((3*species_properties.mass(i))/(4*org_density^2 * pi))^(1/3); %Assume sphere density scaling
                         elseif strcmp(fieldname,'A')
                            species_properties_copy.A = ((3*species_properties.mass(i))/(4*org_density^2 * pi))^(1/3)^2 * pi; %Assume sphere density scaling
                         elseif ~strcmp(fieldname,'mass') && ~strcmp(fieldname,'sym_name')
                            species_properties_copy.(fieldname{:}) = species_properties.(fieldname{:});
                         end
                    elseif length(species_properties.(fieldname{:})) ~= length(species_properties.mass)
                        error("Field values must be provided for each " + ...
                            "mass value, or once. Mistmatch detected " + ...
                            "for field '" + (fieldname{:}) + "'.")
                    end %End fieldname length check
                end %end field loop
                
                % Scale derived properties if they aren't manually set.
                if length(species_properties.amr) ~= num_species
                    species_properties_copy.amr = species_properties_copy.A./species_properties_copy.mass;  % kg
                end
                if length(species_properties.beta) ~= num_species
                    species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
                end
                if length(species_properties.trackable) ~= num_species
                    if species_properties_copy.radius >= TRACKABLE_RADIUS_THRESHOLD
                        species_properties_copy.trackable = true;
                    else
                        species_properties_copy.trackable = false;
                    end
                end
                var_species_props = species_properties_class(species_properties_copy);
                species_instance = species(launch_func, pmd_func, drag_func, var_species_props, scen_properties);
                species_list(i) = species_instance;
                end %end loop over generated species
            obj.species_list = species_list;
        
            %% Set upper and lower bounds for mass bins from their default values
    
            % Sort by mass order
            masses = zeros(1, length(obj.species_list));
            for i = 1:length(obj.species_list)
                masses(i) = obj.species_list(i).species_properties.mass;
            end
            [~, ind] = sort(masses);
            obj.species_list = obj.species_list(1, ind);
            for i = 1:length(obj.species_list)
                if i == 1 %First element (lowest mass)
                    obj.species_list(i).species_properties.mass_ub = 1/2 * (obj.species_list(i).species_properties.mass + obj.species_list(i+1).species_properties.mass);
                elseif i == length(obj.species_list) %Last element (most mass)
                    obj.species_list(i).species_properties.mass_lb = 1/2 * (obj.species_list(i-1).species_properties.mass + obj.species_list(i).species_properties.mass);
                else    
                    obj.species_list(i).species_properties.mass_ub = 1/2 * (obj.species_list(i).species_properties.mass + obj.species_list(i+1).species_properties.mass);
                    obj.species_list(i).species_properties.mass_lb = 1/2 * (obj.species_list(i-1).species_properties.mass + obj.species_list(i).species_properties.mass);
                end
%                 disp([obj.species_list(i).species_properties.sym_name, ...
%                     "mass", num2str(obj.species_list(i).species_properties.mass), ...
%                     "mass lb", num2str(obj.species_list(i).species_properties.mass_lb), ...
%                     "mass ub", num2str(obj.species_list(i).species_properties.mass_ub), ...
%                     ])
            end % End lower and upper bound loop
        end %end constructor
    end % endmethods
end %end classdef