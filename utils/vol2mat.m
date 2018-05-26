function [mat] = vol2mat(data, selected)

selected = logical(selected);

dim = size(selected);

if (length(dim) == 3)
N_samples = size(data, 4);
mat = zeros(length(find(selected)), N_samples);
vol_temp = zeros(dim(1), dim(2), dim(3));

for t = 1: N_samples
    vol_temp = data(:, :, :, t);
    mat(:, t) = vol_temp(selected);
end

else
    N_samples = size(data, 3);
    mat = zeros(length(find(selected)), N_samples);
    slice_temp = zeros(dim(1), dim(2));
    for t = 1: N_samples
        slice_temp = data(:,:,t);
       mat(:, t) = slice_temp(selected); 
    end
end

end