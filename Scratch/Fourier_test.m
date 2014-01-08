% % FFT of turning curve test sets
% 
% % true_spike_rate = [0.66; 1.92; 3.08; 2.56; 0.32; 0.60; 0.50; 0.96; 0.52; 2.53; 2.02; 2.74; 1.04; 0.52; 0.62; 0.46];
% 
%  true_spike_rate = [7.06; 8.64; 8.62; 8.12; 4.04; 1.28; 1.06; 1.00; 0.76; 0.73; 0.46; 0.42; 0.50; 0.58; 1.30; 3.80];
% 
% %  true_spike_rate = [1.64; 1.46; 2.86; 2.46; 2.78; 2.02; 2.30; 1.78; 1.68; 1.78; 2.74; 3.12; 4.4; 2.78;2.48;1.78];
% 
% orientation = 0:22.5:337.5;
% 
% true_spike_rate = circshift(true_spike_rate, 7);
% 
% 
% Y=fft(true_spike_rate);  %% all these index_of_mins were thebesttimulus
% n=length(Y);
% % if n>1;
% %     Y(1)=[];
% % end
% Y;
% 
% phase_1 = angle(Y(2));
% phase_2 = angle(Y(3));
% 
% % if phase_1 <0
% %     phase_1 = phase_1 +2*pi;
% % end
% % if phase_2 <0
% %     phase_2 = phase_2 +2*pi;
% % end
% 
% power = abs(Y(1:floor(n/2))).^2;
% % % % % power = abs(Y(1:ceil(n/2)+1)).^2
% nyquist = 1/2;
% % % freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);
% 
% 
% figure
% h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
% set(h,'LineWidth',2.5);
% hold on
% h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
%     [true_spike_rate(16), true_spike_rate(1)], '-k');
% set(h,'LineWidth',2.5);
% 
% p = polar([phase_1, phase_1], [0, max(true_spike_rate)], '-b');
% set(p,'LineWidth',2.5);
% 
% p = polar([phase_2, phase_2], [0, max(true_spike_rate)], '-r');
% set(p,'LineWidth',2.5);
% 
% keyboard
% 
% % figure;
% % % plot(freq, power)
% % plot(power(2:8))
% keyboard

%%

x = 0:0.1:2*pi;
% y = sin(x+3*pi/2);
y = cos(x - pi/2); 
y = cos(x - 2*pi); 


figure
plot(x,y);

Y = fft(y);

phase = angle(Y(2));
if phase <0
    phase = phase +2*pi;
end
hold on;

plot([phase phase], [0 1]);

keyboard

