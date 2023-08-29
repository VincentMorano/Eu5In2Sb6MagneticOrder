function [scans] = importDataBT7(fileLoc, fileNames, headerLines)
%importDataBT7 Import .bt7 datafiles from a diffraction experiment
%   Provide array of strings indicating the filename (with ".bt7"
%   extension) file location, and number of header lines. Receive a
%   structure array with important variables. Each element of the structure
%   array is a particular datafile (i.e. usually a particular scan). May
%   need to update so get column number of variable first rather than
%   assuming the index.

scans(length(fileNames))=struct;
for i=1:length(fileNames)
    scans(i).file = readcell(strcat(fileLoc, fileNames(i)), 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', headerLines);
    scans(i).fileName = fileNames(i);
    scans(i).T = cell2mat(scans(i).file(:,8));
    scans(i).B = cell2mat(scans(i).file(:,50));
    scans(i).E = cell2mat(scans(i).file(:,40));
    scans(i).det = cell2mat(scans(i).file(:,10));
    scans(i).mon = cell2mat(scans(i).file(:,9));
    scans(i).time = cell2mat(scans(i).file(:,7));
    scans(i).H = cell2mat(scans(i).file(:,46));
    scans(i).K = cell2mat(scans(i).file(:,48));
    scans(i).L = cell2mat(scans(i).file(:,49));
    scans(i).a3 = cell2mat(scans(i).file(:,13));
    scans(i).a4 = cell2mat(scans(i).file(:,6));
    
    scans(i).HKL = [scans(i).H, scans(i).K, scans(i).L];
    
    scans(i).meanT = mean(scans(i).T);
    scans(i).meanB = mean(scans(i).B);
    scans(i).meanE = mean(scans(i).E);
    scans(i).meanTime = mean(scans(i).time);
    scans(i).meanH = round(mean(scans(i).H), 3);
    scans(i).meanK = round(mean(scans(i).K), 3);
    scans(i).meanL = round(mean(scans(i).L), 3);
    scans(i).meanHKL = round(mean(scans(i).HKL), 3);
    scans(i).meanA3 = round(mean(scans(i).a3), 3);
    scans(i).meanA4 = round(mean(scans(i).a4), 3);
    
    scans(i).detErr = sqrt(cell2mat(scans(i).file(:,10)));
    scans(i).monErr = sqrt(cell2mat(scans(i).file(:,9)));
    scans(i).intMon = scans(i).det./scans(i).mon;
    scans(i).intMonErr = (scans(i).det./scans(i).mon).*sqrt((scans(i).detErr./scans(i).det).^2 + (scans(i).monErr./scans(i).mon).^2);
    scans(i).intTime = scans(i).det./scans(i).time;
    scans(i).intTimeErr = scans(i).detErr./scans(i).time;
end
end

