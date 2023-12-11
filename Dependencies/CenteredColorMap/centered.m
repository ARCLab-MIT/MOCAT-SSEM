%% Function to generate white-centered color maps
% Author:  Matteo Courhoud
% Website: https://matteocourthoud.github.io/
% Github:  https://github.com/matteocourthoud



function centeredmap = centered(name, m)
    %   Generate colormap centered on zero
    %
    %   Parameters
    %   ----------
    %   name: str, optional
    %       Name of the colormap, default='RdBu'
    %
    %   m: int, optional
    %       Length of the colormap
    %       
    %   Returns
    %   -------
    %   centeredmap: array
    %       Colormap centered on zero

    % Assign optional parameters
    if nargin == 0
        name = 'RdBu';
    end
    if nargin < 2
        m = size(get(gcf,'colormap'),1);
    end

    % Find limits of the colormap
    lims = get(gca, 'CLim');

    % Find ratio of negative to positive
    if (lims(1) < 0) && (lims(2) > 0)
        % If it has both negative and positive
        ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
        neglen = round(m*ratio);
        poslen = m - neglen;

        % Get two separate colormaps
        posmap = get_map(poslen, name, "positive");
        negmap = get_map(neglen, name, "negative");
        centeredmap = [negmap; posmap];

    elseif lims(1) >= 0
        % Just positive
        centeredmap = get_map(m, name, "positive");

    else
        % Just negative
        centeredmap = get_map(m, name, "negative");

    end
end



function map = get_map(m, name, option)
    %   Generate map
    %
    %   Parameters
    %   ----------
    %   m: int
    %       Length of the colormap
    %       
    %   name: str
    %       Name of the colormap, default='RdBu'
    %
    %   option: str
    %       Positive or negative colors only
    %       
    %   Returns
    %   -------
    %   map: array
    %       Colormap

    % Select map
    switch name
        case 'Spectral'
            map = [158,1,66;213,62,79;215,25,28;244,109,67;252,141,89;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;153,213,148;102,194,165;43,131,186;50,136,189;94,79,162];
        case 'RdYlGn'
            map = [165,0,38;215,48,39;215,25,28;244,109,67;252,141,89;253,174,97;254,224,139;255,255,191;217,239,139;166,217,106;145,207,96;102,189,99;26,150,65;26,152,80;0,104,55];
        case 'RdYlBu'
            map = [165,0,38;215,48,39;215,25,28;244,109,67;252,141,89;253,174,97;254,224,144;255,255,191;224,243,248;171,217,233;145,191,219;116,173,209;44,123,182;69,117,180;49,54,149];
        case 'RdGy'
            map = [103,0,31;178,24,43;202,0,32;214,96,77;239,138,98;244,165,130;253,219,199;255,255,255;224,224,224;186,186,186;153,153,153;135,135,135;64,64,64;77,77,77;26,26,26];
        case 'RdBu'
            map = [103,0,31;178,24,43;202,0,32;214,96,77;239,138,98;244,165,130;253,219,199;255,255,255;209,229,240;146,197,222;103,169,207;67,147,195;5,113,176;33,102,172;5,48,97];
        case 'PuOr'
            map = [127,59,8;179,88,6;230,97,1;224,130,20;241,163,64;253,184,99;254,224,182;255,255,255;216,218,235;178,171,210;153,142,195;128,115,172;94,60,153;84,39,136;45,0,75];
        case 'PrGn'
            map = [64,0,75;118,42,131;123,50,148;153,112,171;175,141,195;194,165,207;231,212,232;255,255,255;217,240,211;166,219,160;127,191,123;90,174,97;0,136,55;27,120,55;0,68,27];
        case 'PiYG'
            map = [142,1,82;197,27,125;208,28,139;222,119,174;233,163,201;241,182,218;253,224,239;255,255,255;230,245,208;184,225,134;161,215,106;127,188,65;77,172,38;77,146,33;39,100,25];
        case 'BrBg'
            map = [84,48,5;140,81,10;166,97,26;191,129,45;216,179,101;223,194,125;246,232,195;255,255,255;199,234,229;128,205,193;90,180,172;53,151,143;1,133,113;1,102,94;0,60,48];
    end

    % Trim for positive/negative only cases
    switch option
        case "positive"
            map = map(8:-1:1,:);
        case "negative"
            map = map(15:-1:8,:);
    end
    map = interp1(linspace(0,1,size(map,1)),map/256,linspace(0,1,m));
    
end