function spatial_conversion = spat_conv(spatial)
% converts pixels/cyc (period) to cyc/mm (freq)

% spatial_conversion = 1/(double(spatial)*9); % 9um = 1 pixel, 31um = 1 degree

% can also potentially convert pixels/cyc (period) to cycles/degree (freq)
spatial_conversion = 1./(spatial.*(9/31)); % 9um = 1 pixel, 31um = 1 degree

end

