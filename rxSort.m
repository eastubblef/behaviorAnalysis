function [ rxSort ] = rxSort(rxMat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Use the matrix of rx times:     col 6: CL = -1; CR = 1; IL = 0; IR = 0;
%col1 = trial number
%col2 = rx times for CL trials
%col3 = rx times for CR trials
%col4 = rx times for IL trials
%col5 = rx times for IR trials

%bin the rx times:
% hist(x,xbins)
load(tb_file);
rxCondensed = [];
for m = 1:length(rxMat(:,2))
    if isfinite(rxMat(m,2))
        rxCondensed(m,1) = rxMat(m,1);
        rxCondensed(m,2) = rxMat(m,2);
        rxCondensed(m,3) = rxMat(m,6);
    else if isfinite(rxMat(m,3))
            rxCondensed(m,1) = rxMat(m,1);
            rxCondensed(m,2) = rxMat(m,3);
            rxCondensed(m,3) = rxMat(m,6);
        else if isfinite(rxMat(m,4))
                rxCondensed(m,1) = rxMat(m,1);
                rxCondensed(m,2) = rxMat(m,4);
                rxCondensed(m,3) = rxMat(m,6);
            else if isfinite(rxMat(m,5))
                    rxCondensed(m,1) = rxMat(m,1);
                    rxCondensed(m,2) = rxMat(m,5);
                    rxCondensed(m,3) = rxMat(m,6);
                else rxCondensed(m,1) = rxMat(m,1);
                     rxCondensed(m,2) = NaN;
                     rxCondensed(m,3) = rxMat(m,6);
                end
            end
        end
    end
end

% Sorted Rx times & drag trial info with ea:
rx = rxCondensed(:,2);
corrVals = rxCondensed(:,3);
[rxSort1, index] = sort(rx, 1, 'ascend');
rxSort = [];
rxSort(:,1) = index;
rxSort(:,2) = rxSort1;
rxSort(:,3) = corrVals(index);

end

