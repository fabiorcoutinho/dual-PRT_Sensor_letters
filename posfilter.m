%apply a median filter in the data
function [out] = posfilter(in,size)
out=in;
    if size>1     
        out(2:end-1,:) = medfilt2(in(2:end-1,:), [size size],'symmetric');
    end
end