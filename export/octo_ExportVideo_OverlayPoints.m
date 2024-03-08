% Make video overlaying point data on original clip
function [] = octo_ExportVideo_OverlayPoints(vidpath,ptdat,outvid)

numpts = size(ptdat,3);
ff = figure(3);
ff.Name = 'Octo';


% If video exists, show frame
v = VideoReader(vidpath);
tempim = readFrame(v);
ff.Position(3) = 0.6*size(tempim,2);
ff.Position(4) = 0.6*size(tempim,1);
imm = imagesc(tempim);
hold on

% Open output video
vw = VideoWriter(outvid,'MPEG-4');
% vw = VideoWriter(outvid,'Motion JPEG AVI');
vw.Quality = 95;
open(vw);

% Loop through
for i=1:size(ptdat,1)

    if i==1
        pp = plot(reshape(ptdat(i,1,:),[numpts 1]),...
            reshape(ptdat(i,2,:),[numpts 1]),'ro','MarkerSize',8);
        hold on
        % Add numbers of parapodia
        for k = 1:numpts
            tx(k) = text(ptdat(i,1,k),ptdat(i,2,k),['  ' num2str(k)],...
                'HorizontalAlignment','left','Color','red','FontSize',11);
        end
        hold off
        currax = gca;
        axis image
        currax.YDir = 'reverse';
        figtxt = text(currax.XLim(1)+10,currax.YLim(2)-10, 'Frame: 1',...
            'HorizontalAlignment','left','VerticalAlignment','bottom',...
            'FontWeight','bold','FontSize',12,'Color','white');
    else
        tempim = readFrame(v);
    end

    %     Update movie background
    imm.CData = tempim;

    % Update point data and labels
    pp.XData = reshape(ptdat(i,1,:),[numpts 1]);
    pp.YData = reshape(ptdat(i,2,:),[numpts 1]);
    for j = 1:numpts
        tx(j).Position = [ptdat(i,1,j), ptdat(i,2,j)];
    end

    % Update frame number
    figtxt.String = ['Frame: ', num2str(i,'%3.0f')];
    drawnow;

    pause(0.001);

    % Capture frame
    frame = getframe(currax);
    frame = imresize(frame.cdata,2);
    writeVideo(vw,frame);
end
close(vw)
end