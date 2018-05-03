function [rxSortNHistory] = rxSortOneBack(rxMat);
%Summary of this function goes here
% 11.2017
% 12.8.17 called by RxAvg_performance; updated to be called by RxAvg_Performance2.m
%% Use the matrix of rx times:     col 6: CL = -1; CR = 1; IL = 0; IR = 0;
%col1 = trial number
%col2 = rx times for CL trials
%col3 = rx times for CR trials
%col4 = rx times for IL trials
%col5 = rx times for IR trials

%if col 6(i) == 0 & col 4 has a value, find when col 2 also had a value on i-1;
%if col 6(i) == 0 & col 5 has a value, find when col 3 also had a value on i-1;
%Currently, this is independent of rx time, but just a readout of incorrects
rxMat6 = rxMat(:,6);
rxMat4 = rxMat(:,4);
rxMat5 = rxMat(:,5);
rxMat3 = rxMat(:,3);
rxMat2 = rxMat(:,2);
rxPrior = [];

for i = 1:length(rxMat2)
    if rxMat6(i) == 0 && isfinite(rxMat4(i)) && isfinite(rxMat2(i-1))
       rxPrior(i) = 1;      %incorrect L trial due to previously chosen/rewarded L trial
    else if rxMat6(i) == 0 && isfinite(rxMat5(i)) && isfinite(rxMat3(i-1))
            rxPrior(i) = 2; %incorrect R trial due to previously chosen/rewarded R trial
        else rxPrior(i) = NaN;
        end
    end
end
rxPrior = rxPrior';
rxMat(:,7) = rxPrior;

%These are incorr trials in which they chose the previously rewarded side on 1 prior back
incorrNpriorL = find(rxPrior == 1);
incorrNpriorR = find(rxPrior == 2);
catIncorrPriors = vertcat(incorrNpriorL, incorrNpriorR);  %trial num where prior dictated current mvmt dir
numCatIncorrPriors = numel(catIncorrPriors);

totNumIncorr = numel(find(rxMat6 == 0));
percNumIncorrPriors = numCatIncorrPriors/totNumIncorr * 100;

rxCondensed = [];
for m = 1:length(rxMat(:,2))
    if isfinite(rxMat(m,2))
        rxCondensed(m,1) = rxMat(m,1);
        rxCondensed(m,2) = rxMat(m,2);
        rxCondensed(m,3) = rxMat(m,6);
        rxCondensed(m,4) = rxMat(m,7);
    else if isfinite(rxMat(m,3))
            rxCondensed(m,1) = rxMat(m,1);
            rxCondensed(m,2) = rxMat(m,3);
            rxCondensed(m,3) = rxMat(m,6);
            rxCondensed(m,4) = rxMat(m,7);
        else if isfinite(rxMat(m,4))
                rxCondensed(m,1) = rxMat(m,1);
                rxCondensed(m,2) = rxMat(m,4);
                rxCondensed(m,3) = rxMat(m,6);
                rxCondensed(m,4) = rxMat(m,7);
            else if isfinite(rxMat(m,5))
                    rxCondensed(m,1) = rxMat(m,1);
                    rxCondensed(m,2) = rxMat(m,5);
                    rxCondensed(m,3) = rxMat(m,6);
                    rxCondensed(m,4) = rxMat(m,7);
                else rxCondensed(m,1) = rxMat(m,1);
                     rxCondensed(m,2) = NaN;
                     rxCondensed(m,3) = rxMat(m,6);
                     rxCondensed(m,4) = rxMat(m,7);
                end
            end
        end
    end
end


% Sorted Rx times & drag trial info with ea:
rx = rxCondensed(:,2);
corrVals = rxCondensed(:,3);
priorVals = rxCondensed(:,4);

[rxSort1, index] = sort(rx, 1, 'ascend');
rxSort = [];
rxSort(:,1) = index;
rxSort(:,2) = rxSort1;
rxSort(:,3) = corrVals(index);
rxSort(:,4) = priorVals(index);

rxSortNHistory = rxSort;

end

