% function [] = plotter(data,snr,delay,scs,cyc,modorder)
% opengl software
data = [];
snr = 20;
% delay = {[0 -1 -2 -3 -4],[0 -2 -4 -6 -8],[0 -4 -5 -9 -12],[0 -3 -7 -11 -15],[0 -5 -10 -15 -20]};
delay = {[0 -5 -10 -15 -20],[0 -3 -7 -10 -14],[0 -4 -5 -10 -12],[0 -2 -4 -5 -8],[0 0 0 0 0]};
scs = [9 10 11 12];
cyc = [1 2 3 4 5 6 7 8 9];
modorder = [1 2 3];

for i=1:length(delay)

    % high values for numerology
    [floorBER(i),floorThroughput(i)] = mainOFDM_Basic2_Fn(scs(1),cyc(end),...
        (delay{i}), modorder(1));

    % low values for numerology
    [ceilBER(i),ceilThroughput(i)] = mainOFDM_Basic2_Fn(scs(end),cyc(1),...
        (delay{i}), modorder(end));
    
%     % adaptive graph
%     line = find(data(:,5) == delay(i));
%     adapBER(i) = data(line(1),1);
%     adapThroughput(i) = data(line(1),2); 
    
    
    
end

delay1 = 1:1:5;
figure;
set(gcf,'color','w')
% set(findall(gcf,'-property','FontSize'),'FontSize',20)
plot(delay1,floorBER,'r');
hold on;
plot(delay1,ceilBER,'g');
plot(delay1,adapBERQ,'b--o','LineWidth',2);
plot(delay1,adapBERdeepQ,'--','LineWidth',2);
plot(delay1,adapBERPPO,'-.','LineWidth',2);
set(gca,'FontSize',14,'FontName','Times New Roman')
hold off;
grid on;
title('BER vs Delay Spread with adaptive learning algorithms')
xlabel('Delay Spread')
ylabel('BER')
legend('Floor','Ceiling','Q-learning','Deep Q-learning','PPO learning','Location','best')

figure;
set(gcf,'color','w')
% set(findall(gcf,'-property','FontSize'),'FontSize',20)
plot(delay1,floorThroughput,'r');
hold on;
plot(delay1,ceilThroughput,'g');
plot(delay1,adapThroughputQ,'b--o','LineWidth',2);
plot(delay1,adapThroughputdeepQ,'--','LineWidth',2);
plot(delay1,adapThroughputPPO,'-.','LineWidth',2);
set(gca,'FontSize',14, 'FontName', 'Times New Roman')
hold off;
grid on;
ylim([0 800])
title('Data rate vs Delay Spread with adaptive learning algorithms')
xlabel('Delay Spread')
ylabel('Data rate (Mbps)')
legend('Floor','Ceiling','Q-learning','Deep Q-learning','PPO learning','Location','best')

% end