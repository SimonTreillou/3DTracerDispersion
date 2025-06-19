function [index] = compute_video_index(image)
    eps=1.e-10;
    Index1=double(image(:,:,1))./(double(image(:,:,3))+eps);
    Index1(Index1>255)=0.0;
    %Index1=image(:,:,1)./image(:,:,3);
    Index2=double(image(:,:,1))./(double(image(:,:,2))+eps);
    Index2(Index2>255)=0.0;
    %Index2=image(:,:,1)./image(:,:,2);
    index=squeeze(Index1.*Index2);
end