% Final set of simulation results for JAR review - Tabulate size
% for all the estimators calculated in FRQSim_FinalSize.m
% Use Gaussian critical values


% Load results
load FRQSim_FINAL_PART1 ;

% Size 
sizefirm = mean(abs(b - ones(nSim,1)*bTest3')./sefirm > 1.96);  
sizestate = mean(abs(b - ones(nSim,1)*bTest3')./sestate > 1.96);
sizesic1 = mean(abs(b - ones(nSim,1)*bTest3')./sesic1 > 1.96); 
sizesic2 = mean(abs(b - ones(nSim,1)*bTest3')./sesic2 > 1.96);  
sizesize = mean(abs(b - ones(nSim,1)*bTest3')./sesize > 1.96);
sizeT8 = mean(abs(b - ones(nSim,1)*bTest3')./seT8 > 1.96); 
sizeT6 = mean(abs(b - ones(nSim,1)*bTest3')./seT6 > 1.96); 
sizeT4 = mean(abs(b - ones(nSim,1)*bTest3')./seT4 > 1.96); 
sizeT2 = mean(abs(b - ones(nSim,1)*bTest3')./seT2 > 1.96); 
sizeSxTstate = mean(abs(bSxT - ones(nSim,1)*bTest3')./seSxTstate > 1.96);  
sizeSIC1xTsic1 = mean(abs(bSIC1xT - ones(nSim,1)*bTest3')./seSIC1xTsic1 > 1.96);  
sizeSIC2xTsic1 = mean(abs(bSIC2xT - ones(nSim,1)*bTest3')./seSIC2xTsic1 > 1.96);  
sizeSIC2xTsic2 = mean(abs(bSIC2xT - ones(nSim,1)*bTest3')./seSIC2xTsic2 > 1.96);  
sizeSIZExTsize = mean(abs(bSIZExT - ones(nSim,1)*bTest3')./seSIZExTsize > 1.96);  
sizeFxT8 = mean(abs(bFxT8 - ones(nSim,1)*bTest3')./seFxT8 > 1.96);  
sizeFxT6 = mean(abs(bFxT6 - ones(nSim,1)*bTest3')./seFxT6 > 1.96);  
sizeFxT4 = mean(abs(bFxT4 - ones(nSim,1)*bTest3')./seFxT4 > 1.96);   
sizeFxT2 = mean(abs(bFxT2 - ones(nSim,1)*bTest3')./seFxT2 > 1.96);   
sizestatetime = mean(abs(b - ones(nSim,1)*bTest3')./sestatetime > 1.96);
sizesic1time = mean(abs(b - ones(nSim,1)*bTest3')./sesic1time > 1.96); 
sizesic2time = mean(abs(b - ones(nSim,1)*bTest3')./sesic2time > 1.96);
sizeSxTstatetime = mean(abs(bSxT - ones(nSim,1)*bTest3')./seSxTstatetime > 1.96);
sizeSIC1xTsic1time = mean(abs(bSIC1xT - ones(nSim,1)*bTest3')./seSIC1xTsic1time > 1.96); 
sizeSIC2xTsic1time = mean(abs(bSIC2xT - ones(nSim,1)*bTest3')./seSIC2xTsic1time > 1.96); 
sizeSIC2xTsic2time = mean(abs(bSIC2xT - ones(nSim,1)*bTest3')./seSIC2xTsic2time > 1.96); 
sizeFMstate = mean(abs(bFMstate - ones(nSim,1)*bTest3')./seFMstate > 1.96); 
sizeFMSIC1 = mean(abs(bFMSIC1 - ones(nSim,1)*bTest3')./seFMSIC1 > 1.96); 
sizeFMsize = mean(abs(bFMsize - ones(nSim,1)*bTest3')./seFMsize > 1.96); 
sizeFMT8 = mean(abs(bFMT8 - ones(nSim,1)*bTest3')./seFMT8 > 1.96);
sizeFMT6 = mean(abs(bFMT6 - ones(nSim,1)*bTest3')./seFMT6 > 1.96);
sizeFMT4 = mean(abs(bFMT4 - ones(nSim,1)*bTest3')./seFMT4 > 1.96);
sizeFMT2 = mean(abs(bFMT2 - ones(nSim,1)*bTest3')./seFMT2 > 1.96);
sizeadapt = mean(abs(badapt - ones(nSim,1)*bTest3')./seadapt > 1.96);

SIZE = [sizefirm ; sizestate ; sizesic1 ; sizesic2 ; sizesize ; sizeT8 ; ...
    sizeT6 ; sizeT4 ; sizeT2 ; sizeSxTstate ; sizeSIC1xTsic1 ; sizeSIC2xTsic1 ; ...
    sizeSIC2xTsic2 ; sizeSIZExTsize ; sizeFxT8 ; sizeFxT6 ; sizeFxT4 ; ...
    sizeFxT2 ; sizestatetime ; sizesic1time ; sizesic2time ; sizeSxTstatetime ; ...
    sizeSIC1xTsic1time ; sizeSIC2xTsic1time ; sizeSIC2xTsic2time ; ...
    sizeFMstate ; sizeFMSIC1 ; sizeFMsize ; sizeFMT8 ; sizeFMT6 ; sizeFMT4 ; ...
    sizeFMT2 ; sizeadapt];


save FRQSim_FINAL_PART2_GAUSS ;
