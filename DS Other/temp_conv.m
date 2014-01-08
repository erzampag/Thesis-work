function temporal_conversion = temp_conv(temporal)
% converts frames/cyc (period) to cyc/sec or Hz (freq)
temporal_conversion = 1./(temporal./120); % 120 frames = 1 sec

end
