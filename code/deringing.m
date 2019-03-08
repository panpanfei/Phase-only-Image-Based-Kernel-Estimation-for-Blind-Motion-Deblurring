function [result] = deringing(y, kernel, lambda_tv, lambda_l0, weight_ring)
%% Removing artifacts in non-blind deconvolution step
[h,w,~] = size(y);
y_pad = wrap_boundary_liu(y, opt_fft_size([h w]+size(kernel)-1));
Latent_tv = zeros(size(y_pad));
for c = 1:size(y,3)
    Latent_tv(:,:,c) = deblurring_adm_aniso(y_pad(:,:,c), kernel, lambda_tv, 1);
end
Latent_tv = Latent_tv(1:h, 1:w, :);

Latent_l0 = L0Restoration(y_pad, kernel, lambda_l0, 2);

Latent_l0 = Latent_l0(1:h, 1:w, :);
%%
diff = Latent_tv - Latent_l0;
bf_diff = bilateral_filter(diff, 3, 0.1);
result = Latent_tv - weight_ring*bf_diff;

result = im2double(result);
end

