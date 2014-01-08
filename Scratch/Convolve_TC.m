
ONMBT = [-0.0157   -0.0151   -0.0317   -0.0327   -0.0502   -0.0615   -0.0932   -0.0852   -0.0962   -0.1020...
    -0.1420   -0.2019   -0.2935   -0.3626   -0.3462   -0.2215...
    0.0308    0.2753    0.3933    0.3807    0.2260    0.0586    0.0059   -0.0028    0.0111];



OFFSSS = [0.0610    0.0679    0.0722    0.0978    0.0429    0.0821    0.0670    0.0748    0.0734    0.0368...
    0.0532    0.0111   -0.0007   -0.0642   -0.1759   -0.3054...
    -0.4232   -0.5135   -0.4162   -0.2898   -0.1792   -0.0492   -0.0147    0.0016    0.0139];


OFFLBT = [0.0146    0.0082    0.0062    0.0120    0.0135    0.0176    0.0234    0.0306    0.0460    0.0707...
    0.1005    0.1396    0.1878    0.2365    0.2789    0.2757...
    0.1942    0.0076   -0.2609   -0.4727   -0.4805   -0.3031   -0.0900   -0.0063    0.0064];


 step = [zeros(1,1000), ones(1,1000)];
 
%  step2 = [ones(1,12), zeros(1,13)];
 step2 = [ones(1,1000), zeros(1,1000)];
 
 
% %  conv1 = conv(step,fliplr(ONMBT));
  conv1 = conv(ONMBT, step, 'same');

 conv2 = conv(OFFSSS, step, 'same');
 
 conv3 = conv(OFFLBT, step, 'same');
 
 conv4 = conv(zeros(1,25),ONMBT);
 
 
 figure
 subplot(3,1,1)
 plot(ONMBT);
 subplot(3,1,2)
 plot(step); ylim([-0.5 1.5]);
 subplot(3,1,3)
 plot(conv1);
 
  figure
 subplot(3,1,1)
 plot(OFFSSS);
 subplot(3,1,2)
 plot(step); ylim([-0.5 1.5]);
 subplot(3,1,3)
 plot(conv2);
 
   figure
 subplot(3,1,1)
 plot(OFFLBT);
 subplot(3,1,2)
 plot(step); ylim([-0.5 1.5]);
 subplot(3,1,3)
 plot(conv3);
 
 
 
 
 