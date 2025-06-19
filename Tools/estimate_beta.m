function [beta] = estimate_beta(v,tmax) 
    image = double(read(v,1));
    for j=2:tmax
        image=image+double(read(v,j));
    end
    image=uint8(image./tmax);
    
    image = image(500:1500,:,:);  % remove the person
    index=compute_video_index(image);
    flat_index=mean(index,1);
    beta = mean(flat_index(500:1000));
end