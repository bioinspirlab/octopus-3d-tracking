%% Information

% Create visualization that looks at the peaks in curvature along each arm,
% and plots them all together

%% Inputs

clips = {"octo_1611_11589_R3",...
    "octo_1611_19083_L2",...
    "octo_1611_19083_L3",...
    "octo_1611_22417_L2",...
    "octo_1611_22417_L3",...
    "octo_1611_22417_L4"};

%% Processing

clrs = jet(numel(clips));

for i = 1:numel(clips)
    matpath = octo_InitializeData(clips{i});
    load(matpath,'armlift','armtrimloc','curvframestart',...
        'curvframeend','armtouchpoint_mm');

    % This is where the actual analysis is performed
    analysisdatapath = octo_AnalyzeSegments(clips{i});
    load(analysisdatapath,'cumdist','curvdat','curvdatbinned','timax',...
        'segdists','segdistbig','curvPeakInd','ptdistarr');


    %% Plot curvature surface

    fsplcurv = figure(8);
    fsplcurv.Position(3) = 440;
    fsplcurv.Position(4) = 480;

    if numel(curvPeakInd)~=numel(timax)
        warning('Different array lengths');
    end

    % Make x axis
%     temp = find(~isnan(curvPeakInd));
    temp = 1:numel(timax);
    temp(temp<curvframestart) = [];
    temp(temp>curvframeend) = [];

    timax_mod = timax(temp);

    figure(fsplcurv)
    plot(timax_mod,interp1(1:numel(cumdist),...
        cumdist,curvPeakInd(temp))-armtouchpoint_mm,'.','MarkerFaceColor',clrs(i,:),'MarkerSize',8);
    hold on

end

%% Modify plot appearance

hold off

ylim([-60 20])
xlim([-30 5])
xlabel('Time from substrate release (s)')
ylabel('Average position along arm (mm) [0 = location of substrate attachment]')
title('Location of peak curvature','Interpreter','none')
legend(clips,'Interpreter','none')
optimizeFig;