function recordinglist = get_recordinglist(subjects, path_data)
% recordinglist = get_recordinglist(subjects, path_data)
%
% Create list with recordings
% 
% Parameters
% ----------
% subjects : cell
%     cell array with subject names
% path_data : string
%     string with the location of the datafiles
%
% Returns
% -------
% recordinglist : table
%     table with recording info
% 


recordinglist = [];

for isub = 1:length(subjects)
    recordings = dir(fullfile(path_data, subjects{isub}, '2*'));
    for irec = 1:length(recordings)
        recinfo = load(fullfile(recordings(irec).folder, recordings(irec).name, 'recinfo.mat'));
        recordinglist = [recordinglist; recinfo.recinfo];
    end
end