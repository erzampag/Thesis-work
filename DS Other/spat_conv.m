function spatial_conversion = spat_conv(spatial)

% spatial period is given in pixels/cyc
% projection of CRT on the array: 9um = 1 pixel
% retinal subtense for mouse: 31um = 1 degree

% converts pixels/cyc (period) to cyc/mm (freq)
% can also potentially convert pixels/cyc (period) to cycles/degree (freq)




% spatial_conversion = 1./(spatial.*9); % parameter in cycles/um

spatial_conversion = 1./(spatial.*(9/31)); %  parameter in cycles/degree

end

