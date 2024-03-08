%% Information

% Function to copy a folder filled with 4D tiffs, but keep only the first
% (depth) channel, and put just one reference file that keeps the other
% channels

%% Function
function outfolder = convert4DtiffsTo1D(infolder,channel)
% infolder should not have a file separator at the end

if nargin<2
    channel = 1;
end

% Determine and create new outfolder
outfolder = [infolder '_stripped'];
mkdir(outfolder);

% Get folder contents
dircont = dir(infolder);
dircont = dircont(~[dircont.isdir]);

% Make sure only tif files are on our list
for i = 1:numel(dircont)
    [~,~,ext] = fileparts(dircont(i).name);
    if strfind(lower(ext),'.tif')==1
        excluded(i) = 0;
    else
        excluded(i) = 1;
    end
end
dircont = dircont(excluded==0);

% Copy a single file as a reference for other channels
ind = strfind(dircont(1).name,'_');
fname = [dircont(1).name(1:ind(end)-1) '_ref.tif'];
infile = [dircont(1).folder,filesep,dircont(1).name];
outfile = [outfolder,filesep,fname];
copyfile(infile,outfile);

% Loop through all files, maintaining their names, but only copying the
% requested channel
for i=1:numel(dircont)
    disp([num2str(i) ' of ' num2str(numel(dircont))]);
    infile = [dircont(i).folder,filesep,dircont(i).name];
    outfile = [outfolder,filesep,dircont(i).name];
    t = Tiff(infile,'r');
    imageData = read(t);
    imageData = imageData(:,:,channel);
    t_out = Tiff(outfile,'w');
    setTag(t_out,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(t_out,'Compression',getTag(t,'Compression'));
    setTag(t_out,'BitsPerSample',getTag(t,'BitsPerSample'));
    setTag(t_out,'SamplesPerPixel',1);
    setTag(t_out,'SampleFormat',getTag(t,'SampleFormat'));
    setTag(t_out,'ImageLength',getTag(t,'ImageLength'));
    setTag(t_out,'ImageWidth',getTag(t,'ImageWidth'));
    setTag(t_out,'PlanarConfiguration',getTag(t,'PlanarConfiguration'));
    write(t_out,imageData);
    close(t_out);
end

disp('Done.')
end