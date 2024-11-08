%% AI Clinician Identifying MIMIC-III Pneumonia Cohort

% Adapted from Matthieu Komorowski, Imperial College London 2015-2019
% Purpose: Identify ICU patients with pneumonia using MIMIC-III

%% ########################################################################
% Import All Data

tic
abx = table2array(readtable('D:/exportdir/abx.csv'));
culture = table2array(readtable('D:/exportdir/culture.csv'));
microbio = table2array(readtable('D:/exportdir/microbio.csv'));
demog = readtable('D:/exportdir/demog.csv');
ce010 = table2array(readtable('D:/exportdir/ce010000.csv'));
ce1020 = table2array(readtable('D:/exportdir/ce1000020000.csv'));
ce2030 = table2array(readtable('D:/exportdir/ce2000030000.csv'));
ce3040 = table2array(readtable('D:/exportdir/ce3000040000.csv'));
ce4050 = table2array(readtable('D:/exportdir/ce4000050000.csv'));
ce5060 = table2array(readtable('D:/exportdir/ce5000060000.csv'));
ce6070 = table2array(readtable('D:/exportdir/ce6000070000.csv'));
ce7080 = table2array(readtable('D:/exportdir/ce7000080000.csv'));
ce8090 = table2array(readtable('D:/exportdir/ce8000090000.csv'));
ce90100 = table2array(readtable('D:/exportdir/ce90000100000.csv'));
labU = [table2array(readtable('D:/exportdir/labs_ce.csv')); table2array(readtable('D:/exportdir/labs_le.csv'))];
MV = table2array(readtable('D:/exportdir/mechvent.csv'));
inputpreadm = table2array(readtable('D:/exportdir/preadm_fluid.csv'));
inputMV = table2array(readtable('D:/exportdir/fluid_mv.csv'));
inputCV = table2array(readtable('D:/exportdir/fluid_cv.csv'));
vasoMV = table2array(readtable('D:/exportdir/vaso_mv.csv'));
vasoCV = table2array(readtable('D:/exportdir/vaso_cv.csv'));
UOpreadm = table2array(readtable('D:/exportdir/preadm_uo.csv'));
UO = table2array(readtable('D:/exportdir/uo.csv'));
toc

%% ########################################################################
% Filter for Pneumonia Patients Using ICD-9 Codes

pneumonia_icd9_codes = ['01166' '01170' '01171' '01172' '01173' '01174' '01175' '01176' '00322'
 '01160' '01161' '01162' '01163' '01164' '01165' '0412' '0413' '0521'
 '0551' '0382' '0203' '0204' '0205' '11505' '11515' '11595' '1304' '1363'
 '0730' '3201' '3523' '48249' '48281' '48282' '48283' '4800' '4801' '4802'
 '4803' '4808' '4809' '481' '4820' '4821' '4822' '48230' '48231' '48232'
 '48239' '48241' '48289' '4829' '4830' '4831' '4838' '4841' '4843' '4845'
 '4846' '4847' '4848' '485' '486' '4870' '5671' '4957' '4958' '4959' '500'
 '502' '503' '504' '505' '5060' '5070' '5071' '5078' '5120' '5121' '51281'
 '51282' '51283' '51289' '51633' '5168' '5169' '5171' '7700' '8600' '8601'
 '8604' '8605' '99731' 'V0382' 'V1261'];
demog = demog(ismember(demog.icd9_code, pneumonia_icd9_codes), :);

%% ########################################################################
% Set Pneumonia Onset Based on ICU Admission Time

onset = zeros(height(demog), 3);

for i = 1:height(demog)
    onset(i, 1) = demog.subject_id(i);
    onset(i, 2) = demog.icustay_id(i);
    onset(i, 3) = demog.intime(i);
end

onset = array2table(onset, 'VariableNames', {'subject_id', 'icustay_id', 'onset_time'});

%% ########################################################################
% Initial Data Manipulations

bacterio = [microbio; culture];
demog.morta_90(isnan(demog.morta_90)) = 0;
demog.morta_hosp(isnan(demog.morta_hosp)) = 0;
demog.elixhauser(isnan(demog.elixhauser)) = 0;
inputMV(:,8) = inputMV(:,7) .* inputMV(:,6) ./ inputMV(:,5);

%% ########################################################################
% Data Reformatting for CHARTEVENTS, LABS, and MECHVENT

reformat = NaN(2000000, 68);
qstime = zeros(100000, 4);
winb4 = 49;
winaft = 25;
irow = 1;

tic
for icustayid = 1:100000
    qst = onset(icustayid, 3);
    
    if qst > 0
        d1 = table2array(demog(demog.icustay_id == icustayid, [11 5]));

        if d1(1) > 6574

            if icustayid < 10000
                temp = ce010(ce010(:,1) == icustayid, :);
            elseif icustayid < 20000
                temp = ce1020(ce1020(:,1) == icustayid, :);
            elseif icustayid < 30000
                temp = ce2030(ce2030(:,1) == icustayid, :);
            elseif icustayid < 40000
                temp = ce3040(ce3040(:,1) == icustayid, :);
            elseif icustayid < 50000
                temp = ce4050(ce4050(:,1) == icustayid, :);
            elseif icustayid < 60000
                temp = ce5060(ce5060(:,1) == icustayid, :);
            elseif icustayid < 70000
                temp = ce6070(ce6070(:,1) == icustayid, :);
            elseif icustayid < 80000
                temp = ce7080(ce7080(:,1) == icustayid, :);
            elseif icustayid < 90000
                temp = ce8090(ce8090(:,1) == icustayid, :);
            else
                temp = ce90100(ce90100(:,1) == icustayid, :);
            end

            ii = temp(:,2) >= qst - (winb4 + 4) * 3600 & temp(:,2) <= qst + (winaft + 4) * 3600;
            temp = temp(ii, :);
            
            ii = labU(:,1) == icustayid;
            temp2 = labU(ii, :);
            ii = temp2(:,2) >= qst - (winb4 + 4) * 3600 & temp2(:,2) <= qst + (winaft + 4) * 3600;
            temp2 = temp2(ii, :);
            
            ii = MV(:,1) == icustayid;
            temp3 = MV(ii, :);
            ii = temp3(:,2) >= qst - (winb4 + 4) * 3600 & temp3(:,2) <= qst + (winaft + 4) * 3600;
            temp3 = temp3(ii, :);

            t = unique([temp(:,2); temp2(:,2); temp3(:,2)]);

            if t
                for i = 1:numel(t)
                    ii = temp(:,2) == t(i);
                    col = temp(ii, 3);
                    value = temp(ii, 4);
                    reformat(irow, 1) = i;
                    reformat(irow, 2) = icustayid;
                    reformat(irow, 3) = t(i);
                    reformat(irow, 3 + col) = value;

                    ii = temp2(:,2) == t(i);
                    col = temp2(ii, 3);
                    value = temp2(ii, 4);
                    reformat(irow, 31 + col) = value;

                    ii = temp3(:,2) == t(i);
                    if nansum(ii) > 0
                        value = temp3(ii, 3:4);
                        reformat(irow, 67:68) = value;
                    else
                        reformat(irow, 67:68) = NaN;
                    end

                    irow = irow + 1;
                end

                qstime(icustayid, 1) = qst;
                qstime(icustayid, 2) = t(1);
                qstime(icustayid, 3) = t(end);
                qstime(icustayid, 4) = d1(2);
            end
        end
    end
end
toc

close(h);
reformat(irow:end, :) = [];

%% ########################################################################
% Outlier Removal and Initial Data Processing

reformat = deloutabove(reformat, 5, 300);
reformat = deloutabove(reformat, 8, 250);
reformat = deloutabove(reformat, 9, 300);
reformat = deloutbelow(reformat, 10, 0);
reformat = deloutabove(reformat, 10, 200);
reformat = deloutbelow(reformat, 11, 0);
reformat = deloutabove(reformat, 11, 200);
reformat = deloutabove(reformat, 12, 80);
reformat = deloutabove(reformat, 13, 150);
ii = reformat(:, 13) > 100;
reformat(ii, 13) = 100;
reformat = deloutabove(reformat, 14, 90);
ii = reformat(:, 15) > 25 & reformat(:, 15) < 45;
reformat(ii, 14) = reformat(ii, 15);
reformat(ii, 15) = NaN;
ii = reformat(:, 14) > 70;
reformat(ii, 15) = reformat(ii, 14);
reformat(ii, 14) = NaN;

%% ########################################################################
% Save the Final Pneumonia Cohort

writetable(onset, 'pneumonia_mimiciii.csv', 'Delimiter', ',');
reformat_table = array2table(reformat);
writetable(reformat_table, 'pneumonia_reformat.csv', 'Delimiter', ',');
