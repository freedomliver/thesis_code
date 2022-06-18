function newdata = ReadSIF_all(filename)
% retrieve all frame in the SIF file
% 2020/10/10 qifengfeng
%
% filename: the path and name of the file
% newdata: the array of the pattern

atsif_setfileaccessmode(0); % read all
rc=atsif_readfromfile(filename); % open a SIF file
if (rc == 22002) % read successfully 
    signal=0; % determine if a particular data source is present in the SIF file; 0-Signal
    [~,present]=atsif_isdatasourcepresent(signal);
    if present
        [~,num_frames]=atsif_getnumberframes(signal); % retrieve the number of frames in the SIF file
        if (num_frames > 0)
            [~,size]=atsif_getframesize(signal); % retrieve the number of pixels in each frame
            [~,left,bottom,right,top,hBin,vBin]=atsif_getsubimageinfo(signal,0); % retrieve the information about the first sub-image 
            [~,data]=atsif_getallframes(signal,size); % retrieve all frame in the SIF file
            [~,pattern]=atsif_getpropertyvalue(signal,'ReadPattern'); % retrieve image information from the SIF file
            if(pattern == '4')
                width = ((right - left)+1)/hBin;
                height = ((top-bottom)+1)/vBin;
                newdata=reshape(data,width,height);
            end
        end
    end
    atsif_closefile; %  close the currently opened  SIF file
else
    error_information=['ERROR: Could not load file. Error number: ',num2str(rc)];
    error(error_information);
    
end