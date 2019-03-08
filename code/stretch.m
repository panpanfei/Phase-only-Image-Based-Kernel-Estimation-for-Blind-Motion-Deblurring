function a = stretch (im)
    vmax = max(im(:));
    vmin = min(im(:));
    scale = 1.0 / (vmax - vmin);
    a = (im - vmin) * scale;
end
