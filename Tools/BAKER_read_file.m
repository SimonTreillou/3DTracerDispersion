function [arr]=BAKER_read_file(fname)
    arr.data = load(fname);
    fid = fopen(fname);
    tline = fgetl(fid);
    while contains(tline,'%')
        disp(tline)
        if contains(tline,"X:")
            arr.x = str2double(tline(6:end));
        elseif contains(tline,"Y:")
            arr.y = str2double(tline(6:end));
        elseif contains(tline,"Z:")
            arr.z = str2double(tline(6:end));
        end
        tline = fgetl(fid);
    end
    fclose(fid);
return 