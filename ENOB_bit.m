fileId = fopen('E:\CAN\CAN_disturbance_4\TH7011_V02\ADC_sample\ADC_matlab\ENOB\data\ENOB_B.txt', 'r');

data = textscan(fileId, '%s');

fclose(fileId);

% get ADC data
ADC_data = data{1};

num_values = length(ADC_data);

% Transform ADC data from hex to decimal considering sign bit
ADC_data_d = arrayfun(@(x) convertHexToSignedDec(x{1}), ADC_data);

% delet offset error
signal_value = ADC_data_d -7;

% ADC sample frequency
Fs = 40e6;

% ADC bits number
num_bit = 8;

% set time_vector
time_vector = (0:25:25*(length(signal_value)-1)) * 1e-3; 

subplot(2, 1, 1);
% Plot with DisplayName property for legend
plot(time_vector(1:60), signal_value(1:60), 'Color', [0.85, 0.9, 0.3], 'Linewidth', 1.5, 'DisplayName', 'Channel B'); 
xlabel('Time (us)');
ylabel('Digital Value');
title('Time Domain');
% set background black
set(gca, 'Color', 'k');
ax = gca;
ax.FontSize = 22;
% Add legend and set its color
lgd = legend;  % This will show legend with names defined in 'DisplayName'
lgd.TextColor = 'white';  % Setting the legend text color to white
% limit x coordinate
xlim([time_vector(1) time_vector(60)]);
ylim([-128,128]);

% get 3dB band width by empirical formula
span = max(round(num_values/200),5);

% Approximate search span for harmonics on each side
spanh = 5;

% get normalization value
data33 = signal_value(1:num_values)/(2^(num_bit-1)-1);

hammingWindow = zeros(num_values, 1);
for i = 1:num_values
    hammingWindow(i) = 0.54 - 0.46 * cos(2 * pi * (i-1) / (num_values-1));
end
windowedSignal = data33 .* hammingWindow;

%% calculate fft
Y = fft(windowedSignal(1:65536));

V1 = abs(Y); 
P1_dB = 20*log10(V1);
P2 = V1 .* V1;
P2_dB = 10*log10(P2);

max_voltage = max(abs(data33));
max_voltagedb = 20*log10(max_voltage/1);

bili_factor_1=1/max_voltage;
bili_factor=bili_factor_1^2;

maxdb=max(P1_dB(1:num_values/2));
y1 = P1_dB(1:num_values/2) - maxdb;

f = (0:(length(Y)/2)) * (Fs / length(Y)) / 1e3;

subplot(2, 1, 2);
plot(f(1:32768), y1);
title('Frequency domain');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
ax = gca;
ax.FontSize = 22;
grid on;
xlim([f(1) f(32768)]);
ylim([-100,0]);

fin = find(P1_dB(1:num_values/2) == maxdb);     %find the signal bin number

Fh= zeros(1:8);     % frequency of harmonics
Ph= zeros(1:8);     % power of harmonics

for har_num=1:8      %1到15次谐波分量
    tone = rem((har_num * fin + 1) / num_values, 1);  
    %各次谐波对应的能量
    if tone>0.5
            tone=1-tone;      %对称
    end
Fh(har_num) = tone;
%谐波，在Fh中加入tone
if round(tone*num_values)-spanh>0
har_peak=max(P2(round(tone*num_values)-spanh:round(tone*num_values)+spanh));%谐波附近+-spanh最大值
har_bin=find(P2(round(tone*num_values)-spanh:round(tone*num_values)+spanh)==har_peak);%最大值对应频率
    har_bin=har_bin+round(tone*num_values)-spanh-1;
    Ph(har_num)=sum(P2(round(tone*num_values)-spanh:round(tone*num_values)+spanh)); 
%  har_bin-spanh:har_bin+spanh;
else 
har_peak=max(P2(1:round(tone*num_values)+spanh));%谐波附近+-spanh最大值
har_bin=find(P2(1:round(tone*num_values)+spanh)==har_peak); %最大值对应频率
    har_bin=har_bin+round(tone*num_values)-spanh-1;
    Ph(har_num)=sum(P2(1:har_bin+spanh)); %旧版本，每次迭代的时候append 数据，速度慢，现改称预设内存大小，根据数组下标赋值
end
end

Pdc=sum(P2(1:span));                                            % 直流能量
Ps=bili_factor*sum(P2(fin-span:fin+span));                      % 信号能量，前后span个都算作信号能量，
Pd=bili_factor*sum(Ph(2:8));                                    % 谐波能量
Pn=sum(P2(1:num_values/2))-Pdc-Pd/bili_factor-Ps/bili_factor;   % 噪声能量%

fin_MHz=fin/num_values*Fs;                      % 信号频率
SNR=10*log10(Ps/Pn);                            % 信噪比
SINAD=10*log10(Ps/(Pn+Pd));                     % 信纳比
SFDR=10*log10(Ph(1)/max(Ph(2:8)));              % 无杂散动态范围
ENOB=(SINAD-1.76+20*log10(1/max_voltage))/6.02; % 幅度校正后有效位

function decValue = convertHexToSignedDec(hexValue)
    % Convert hex to a 16-bit binary string
    binaryStr = dec2bin(hex2dec(hexValue), 8);

    % Check if the number is negative (highest bit is 1)
    if binaryStr(1) == '1'
        % Convert from two's complement and negate
        decValue = -(2^8 - bin2dec(binaryStr));
    else
        % Convert directly to decimal
        decValue = bin2dec(binaryStr);
    end
end