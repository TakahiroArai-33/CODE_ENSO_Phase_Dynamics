function tNum = customDatemnth(startDate, nMonths)
    % customDatemnth: Add months to a given date
    % Inputs:
    %   startDate - Scalar or array, base date in datenum format
    %   nMonths   - Array of months to add
    % Outputs:
    %   tNum      - Array of resulting dates in datenum format

    % Convert the start date to datetime
    dt = datetime(startDate, 'ConvertFrom', 'datenum');
    
    % Preallocate output
    tNum = zeros(size(nMonths));
    
    for i = 1:numel(nMonths)
        % Add months and handle month overflow
        newDate = dt + calmonths(nMonths(i));
        
        % Convert back to datenum
        tNum(i) = datenum(newDate);
    end
end
