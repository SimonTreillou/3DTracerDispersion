function [res] = BAKER_read_data(path_directory,filetype)
    %
    % Read data from various captors in the same directory
    %
    original_files=dir([path_directory+filetype+'*']); 
    for k=1:length(original_files)
        filename=[path_directory+original_files(k).name];
        fieldname = erase(original_files(k).name,".txt");
        res.(fieldname)=BAKER_read_file(filename);
    end
return
