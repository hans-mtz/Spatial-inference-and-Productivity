
clc 
clear

nSim=1000;

load('SizePt1_1000X2000_bsize3_17Obs_sim_1000.mat')

Bsize_firm    = mean(abs(b-ones(nSim,1)*bTest3')./se_firm > bcv_firm)
Bsize_state   = mean(abs(b-ones(nSim,1)*bTest3')./se_state > bcv_state)
Bsize_sic1    = mean(abs(b-ones(nSim,1)*bTest3')./se_sic1 > bcv_sic1)
Bsize_sic2    = mean(abs(b-ones(nSim,1)*bTest3')./se_sic2 > bcv_sic2)
Bsize_size    = mean(abs(b-ones(nSim,1)*bTest3')./se_size > bcv_size)

disp([ 'firm & firm, year                     & ' num2str(Bsize_firm(1,1),'%12.3f') ' & '  num2str(Bsize_firm(1,2),'%12.3f') ' & ' ...
     num2str(Bsize_firm(1,3),'%12.3f') ' & '  num2str(Bsize_firm(1,4),'%12.3f') ' & '  num2str(Bsize_firm(1,5),'%12.3f') ' & ' ...
     num2str(Bsize_firm(1,6),'%12.3f') ' & '  num2str(Bsize_firm(1,7),'%12.3f') ' & '  num2str(Bsize_firm(1,8),'%12.3f') ' & '  num2str(Bsize_firm(1,9),'%12.3f') '\\ '])
 
 disp([ 'state & firm, year                   & ' num2str(Bsize_state(1,1),'%12.3f') ' & '  num2str(Bsize_state(1,2),'%12.3f') ' & ' ...
     num2str(Bsize_state(1,3),'%12.3f') ' & '  num2str(Bsize_state(1,4),'%12.3f') ' & '  num2str(Bsize_state(1,5),'%12.3f') ' & ' ...
     num2str(Bsize_state(1,6),'%12.3f') ' & '  num2str(Bsize_state(1,7),'%12.3f') ' & '  num2str(Bsize_state(1,8),'%12.3f') ' & '  num2str(Bsize_state(1,9),'%12.3f') '\\ '])
 
 disp([ 'one-digit SIC & firm, year           & ' num2str(Bsize_sic1(1,1),'%12.3f') ' & '  num2str(Bsize_sic1(1,2),'%12.3f') ' & ' ...
     num2str(Bsize_sic1(1,3),'%12.3f') ' & '  num2str(Bsize_sic1(1,4),'%12.3f') ' & '  num2str(Bsize_sic1(1,5),'%12.3f') ' & ' ...
     num2str(Bsize_sic1(1,6),'%12.3f') ' & '  num2str(Bsize_sic1(1,7),'%12.3f') ' & '  num2str(Bsize_sic1(1,8),'%12.3f') ' & '  num2str(Bsize_sic1(1,9),'%12.3f') '\\ '])
 
 disp([ 'two-digit SIC & firm, year           & ' num2str(Bsize_sic2(1,1),'%12.3f') ' & '  num2str(Bsize_sic2(1,2),'%12.3f') ' & ' ...
     num2str(Bsize_sic2(1,3),'%12.3f') ' & '  num2str(Bsize_sic2(1,4),'%12.3f') ' & '  num2str(Bsize_sic2(1,5),'%12.3f') ' & ' ...
     num2str(Bsize_sic2(1,6),'%12.3f') ' & '  num2str(Bsize_sic2(1,7),'%12.3f') ' & '  num2str(Bsize_sic2(1,8),'%12.3f') ' & '  num2str(Bsize_sic2(1,9),'%12.3f') '\\ '])
 
  disp([ 'size category & firm, year          & ' num2str(Bsize_size(1,1),'%12.3f') ' & '  num2str(Bsize_size(1,2),'%12.3f') ' & ' ...
     num2str(Bsize_size(1,3),'%12.3f') ' & '  num2str(Bsize_size(1,4),'%12.3f') ' & '  num2str(Bsize_size(1,5),'%12.3f') ' & ' ...
     num2str(Bsize_size(1,6),'%12.3f') ' & '  num2str(Bsize_size(1,7),'%12.3f') ' & '  num2str(Bsize_size(1,8),'%12.3f') ' & '  num2str(Bsize_size(1,9),'%12.3f') '\\ '])