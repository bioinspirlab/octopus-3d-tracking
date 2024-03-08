%% Information

% Batch script to process a whole list of clips (so we don't have to do
% this later)

%%
% All used sequences:
cliparr = {"O10_6946_L1","O10_7875_R1","O14_24216_L1","O14_24216_L2",...
    "O14_24216_L3","O15_1611_7882_R3","O15_1611_11589_R3",...
    "O15_1611_13209_L3_1","O15_1611_13209_L3_2","O15_1611_13209_L4",...
    "O15_1611_15512_L1","O15_1611_19083_L2","O15_1611_19083_L3",...
    "O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4"};
% All sequences:
% cliparr = {"O10_6946_L1","O10_7875_R1","O14_24216_L1","O14_24216_L2",...
%     "O14_24216_L3","O15_1611_5596","O15_1611_7318_L1","O15_1611_7318_L2",...
%     "O15_1611_7882_R1","O15_1611_7882_R3","O15_1611_7882_R4",...
%     "O15_1611_11589_R3","O15_1611_11589_R4","O15_1611_13209_L2",...
%     "O15_1611_13209_L3_1","O15_1611_13209_L3_2","O15_1611_13209_L4",...
%     "O15_1611_15512_L1","O15_1611_19083_L1","O15_1611_19083_L2",...
%     "O15_1611_19083_L3","O15_1611_22417_L1",...
%     "O15_1611_22417_L2","O15_1611_22417_L3","O15_1611_22417_L4"};
% cliparr = {"O15_1611_15512_L1","O14_24216_L2","O15_1611_13209_L4"};

%% Loop!

for i = 1:numel(cliparr)
    disp(append('Now processing ',cliparr{i}))

    try
%         octo_PreProcess(cliparr{i});
%         octo_ExportVideo_OverlayPointsAndSpline(cliparr{i},"overwrite");
%         octo_ExportVideo_3Dviews(cliparr{i});
        octo_ExportVideo_Quad3D(wcliparr{i},"overwrite");
    catch ME
        disp(append(cliparr{i}, ' failed with error ', ME.identifier));
        disp(ME.message);
    end
end