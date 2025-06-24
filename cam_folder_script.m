%% Script to split tiff files into CamA and CamB folders
% Remember to change campath AND spec before running
% epoch number
epochnum='1726691401451'; % 1726504201906, 1726687801658 1726507801090
% paths 
genpath='/Volumes/Elements';% path to external drive, path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
epochpath=append(genpath,'/RawImages/',epochnum); % path to folder with photos
campath=append(epochpath,'/CamA'); % '/CamB' '/CamA'

% Loop through photos in folder and sort into CamA and CamB
spec='*CameraA*'; % '*CameraB*' '*CameraA'
filePattern=fullfile(epochpath,spec); %CameraA_1727362860214.tiff
listofFiles=dir(filePattern);

for k=1:length(listofFiles)
    baseFilename=listofFiles(k).name;
    fullFilename=fullfile(listofFiles(k).folder,baseFilename);
    %fprintf(1, 'Now moving %s\n', fullFilename);
    % check for folder and save
    if isfolder(campath) == true
        movefile(fullFilename);
        %fprintf(k);
        movefile(append('/Users/bagaenzl/Desktop/MasonBEAST_GIT/',baseFilename),campath);
        %fprintf('Moved %s\n',baseFilename,'to %s\n',campath);
    else
        mkdir(campath);
        fprintf('Created new folder for %s\n',spec);
        movefile(fullFilename);
        movefile(append('/Users/bagaenzl/Desktop/MasonBEAST_GIT/',baseFilename),campath);
        %fprintf('Moved %s\n',baseFilename,'to %s\n',campath);
    end
   
end

campath=append(epochpath,'/CamB'); % '/CamB' '/CamA'

% Loop through photos in folder and sort into CamA and CamB
spec='*CameraB*'; % '*CameraB*' '*CameraA'
filePattern=fullfile(epochpath,spec); %CameraA_1727362860214.tiff
listofFiles=dir(filePattern);

for k=1:length(listofFiles)
    baseFilename=listofFiles(k).name;
    fullFilename=fullfile(listofFiles(k).folder,baseFilename);
    fprintf(1, 'Now moving %s\n', fullFilename);
    % check for folder and save
    if isfolder(campath) == true
        movefile(fullFilename);
        %frpintf(k);
        movefile(append('/Users/bagaenzl/Desktop/MasonBEAST_GIT/',baseFilename),campath);
        %fprintf('Moved %s\n',baseFilename,'to %s\n',campath);
    else
        mkdir(campath);
        %fprintf('Created new folder for %s\n',spec);
        movefile(fullFilename);
        movefile(append('/Users/bagaenzl/Desktop/MasonBEAST_GIT/',baseFilename),campath);
        %fprintf('Moved %s\n',baseFilename,'to %s\n',campath);
    end
   
end
