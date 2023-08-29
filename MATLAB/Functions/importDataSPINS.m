function [scans] = importDataSPINS(fileLoc,fileNames,T,mon,H,K,L,E)
%importDataSPINS Import .ng5 datafiles from a diffraction experiment
%   Provide array of strings indicating the filename (with ".ng5"
%   extension) file location, and number of header lines. Receive a
%   structure array with important variables. Each element of the structure
%   array is a particular datafile (i.e. usually a particular scan). May
%   need to update so get column number of variable first rather than
%   assuming the index.

scans(length(fileNames))=struct;
for i=1:length(fileNames)
    scans(i).file=readcell(strcat(fileLoc, fileNames(i)), 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 1);
    scans(i).fileName=fileNames(i);
    scans(i).T=T(i);
    scans(i).det=cell2mat(scans(i).file(:,3));
    scans(i).mon=mon(i);
    scans(i).H=H(i);
    scans(i).K=K(i);
    scans(i).L=L(i);
    scans(i).E=E(i);
    scans(i).a3=cell2mat(scans(i).file(:,2));
    scans(i).N=length(scans(i).det);
    
    scans(i).HKL=[scans(i).H, scans(i).K, scans(i).L];
    scans(i).meanA3=round(mean(scans(i).a3), 3);
    scans(i).detErr=sqrt(scans(i).det);
    scans(i).monErr=sqrt(scans(i).mon);
    scans(i).intMon=scans(i).det./scans(i).mon;
    scans(i).intMonErr=(scans(i).det./scans(i).mon).*sqrt((scans(i).detErr./scans(i).det).^2+(scans(i).monErr./scans(i).mon).^2);
end
end

