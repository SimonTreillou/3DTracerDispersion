function [res] = filter_array(array,time,interval)
    if length(size(array))==4
        [Mt,Mz,My,Mx]=size(array);
        res=array.*0;
        for l=1:Mz
            flat=reshape(array(:,l,:,:),[Mt My*Mx]);
            count1=timeseries(flat,time); 
            idealfilter_countn = idealfilter(count1,interval,'pass'); clear count1;
            flatdp=idealfilter_countn.data;
            flatdp=reshape(flatdp,[Mt 1 My Mx]);
            clear flat idealfilter_countn
            res(:,l,:,:)=flatdp;
            clear flatdp
        end
    elseif length(size(array))==3
        [Mt,My,Mx]=size(array);
        flat=reshape(array(:,:,:),[Mt My*Mx]);
        count1=timeseries(flat,time); 
        idealfilter_countn = idealfilter(count1,interval,'pass'); 
        flatdp=idealfilter_countn.data;
        res=reshape(flatdp,[Mt 1 My Mx]);
        clear flat count1 idealfilter_countn flatdp
    end
end