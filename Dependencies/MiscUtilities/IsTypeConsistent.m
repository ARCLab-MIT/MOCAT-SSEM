function [all_same] = IsTypeConsistent(myCell, type)
% If no type is provided, returns true if all elements of myCell are the same type.
% If type is provided as arg, return true if all elements of mc Cell are the same
% type and that type is the passed type arg.
% else returns false
% adapted from https://stackoverflow.com/questions/48317906/fastest-way-to-get-class-types-of-elements-of-a-cell-array
%[
    [nrows, ncols] = size(myCell);
    types = cell(nrows, ncols);
    for ci = 1:ncols
        for cr = 1:nrows
            cellclass = class(myCell{cr, ci});
            types{cr, ci} = cellclass; 
        end
    end
    type_list = unique(types);
    all_same = numel(type_list) == 1;

    if nargin==2 & all_same == true
        if type_list(1) == type
            all_same = true;
        else
            all_same = false;
        end
    end
end