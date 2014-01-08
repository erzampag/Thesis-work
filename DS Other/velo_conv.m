
function velocity_conversion = velo_conv(spatial, temporal)
% converts pixels and frames to degrees/sec
velocity_conversion = (spatial).*(9/31) ./ ((temporal).*120);

end
