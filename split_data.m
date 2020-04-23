function [X_id, X_val, Y_id, Y_val] = split_data(X, Y)

%     perc_val = 0.5
%     perc_id  = 1 - perc_val;
%     rng(25);
%     idx_X   = randperm(size(X, 2));
%     idx_Y   = randperm(size(Y, 2));
% 
% 
%     indexToGroup11 = (idx_X<= perc_id * size(X,2));
%     indexToGroup12 = (idx_X>perc_val * size(X,2));
%     indexToGroup21 = (idx_Y<= perc_id * size(Y, 2));
%     indexToGroup22 = (idx_Y>perc_val * size(Y, 2));
% 
%     X_id    = X(:,indexToGroup11);
%     X_val   = X(:,indexToGroup12);
%      
%     Y_id    = Y(:,indexToGroup21);
%     Y_val   = Y(:,indexToGroup22);

    X_val = X(1:2,2:2:end);     % even column matrix
    X_id = X(1:2,1:2:end);      % odd column matrix
    Y_val = Y(:,2:2:end);       % even column matrix
    Y_id = Y(:,1:2:end);        % odd column matrix
    

end 