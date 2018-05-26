function s = getSlice(ind_k, mask)

[N_row, N_col, N_slice] = size(mask);

slice_dim = zeros(N_slice, 1);
for k = 1: N_slice
    slice_dim(k) = length(find(mask(:,:,1:k)));
end

s = zeros(size(ind_k));
for k = 1: length(ind_k)
    s(k) = find(slice_dim>=ind_k(k), 1, 'first');
end

s = reshape(s, length(s), 1);
end