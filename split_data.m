% SPLIT_DATA separates the dataset into an identification and validation
% dataset

function [X_id, X_val, Y_id, Y_val] = split_data(X, Y)
    
    X_id   = X(1:2,1:2:end);     % odd column matrix
    X_val  = X(1:2,2:2:end);     % even column matrix
    
    Y_id   = Y(:,1:2:end);       % odd column matrix
    Y_val  = Y(:,2:2:end);       % even column matrix

end 