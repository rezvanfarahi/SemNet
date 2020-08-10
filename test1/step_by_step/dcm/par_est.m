clear
bms_file='/imaging/rf02/TypLexMEG/dcm/latest/150ms_32models_2ROIs/typlex/BMS.mat';
% bms_file='/imaging/rf02/TypLexMEG/dcm/latest/150ms_7models/semloc/BMS.mat';
load(bms_file)
params.type='A'; params.hier='forward';
params.to=2; params.from=1; params.xt=1;

exp(BMS.DCM.ffx.bma.mEp.A{1})
exp(BMS.DCM.ffx.bma.mEp.A{2})
exp(BMS.DCM.ffx.bma.mEp.B{1})