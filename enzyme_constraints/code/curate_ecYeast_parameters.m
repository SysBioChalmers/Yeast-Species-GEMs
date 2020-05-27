function model = curate_ecYeast_parameters(model)
enzymes = {'P09624' ...
           'P16387' 'P32473' ...
           'Q00711' 'P21801' 'P37298' 'P33421' ...
           'P06169' 'P16467' 'P26263' ...
           'P00330' 'P00331' 'P07246' 'P10127' 'P25377' 'P38113' 'Q04894'};
%PDAs: 'P09624' 'P16387' 'P32473'  41.65
%SDHs: 'Q00711' 'P21801' 'P37298' 'P33421'  219.3      
%PDCs: 'P06169' 'P16467' 'P26263' 62     
%ATPs: 'P61829''P07251' 'P00830' 120
%ADHs: 'P00330' 'P00331' 'P07246' 'P10127' 'P25377' 'P38113' 'Q04894' 900
kcats   = [41.65 41.65 41.65 ...
           219.3 219.3 219.3 219.3 ...
           62 62 62 ...
           900 900 900 900 900 900 900];
cd GECKO/geckomat/utilities
for i=1:length(enzymes)
    pIdx = strcmpi(model.metNames,['prot_' enzymes{i}]);
    [kecat,rxnIdx,rxnName,~] = getKcat(model,enzymes{i});
    kcat = kcats(i);
    model.S(pIdx,rxnIdx) = -1/(3600*kcat);
 end
cd ../../..
end