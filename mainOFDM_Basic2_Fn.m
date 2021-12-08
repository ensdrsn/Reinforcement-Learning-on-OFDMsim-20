function [BER,throughput] = mainOFDM_Basic2_Fn(scs,cyclic,delay,modulation)

numberofSC = 2^scs;
cyclicPrefix = (2^cyclic);
% delay = [0 3 5 6 8]./delayS; 
% delaySpread = {[0 0 0 0 0],[0 0 0 10 0],[0 0 20 10 0],[0 0 30 10 0],[0 0 40 10 0]};
% delayNum = randi(length(delaySpread));
% delay = cell2mat(delaySpread(delayNum));

% delayPower = [0 -8 -17 -21 -25].*delayS;
ebno = 20;
modulation = 2*modulation;

[BER,throughput] = fn_OFDM_Basic2(numberofSC, cyclicPrefix, ebno, delay, modulation);
fprintf('Delay Spread:  \n ----------------------\n');
end

