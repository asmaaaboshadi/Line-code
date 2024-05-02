clc;
clear;

% Number of bits to generate
num_bits = 100;

% Generate a random binary sequence
rand_bits = randi([0, 1], 1, num_bits);

% Display the original data
figure;
stairs(1:num_bits, rand_bits, 'LineWidth', 2.5); % Ensure x-axis is correct
title('Original Data');
xlabel('Bit Index');
ylabel('Value (0 or 1)');
%% Create the NRZ signal
% Initialize NRZ signal
% The signal has the same length as the number of bits
nrz_signal = zeros(1, num_bits); 
for i = 1:num_bits
    if rand_bits(i) == 1
        nrz_signal(i) = 1; % Polar NRZ: "1" is positive
    else
        nrz_signal(i) = -1; % Polar NRZ: "0" is negative
    end
end

% Plot the NRZ signal
figure;
stairs(1:num_bits, nrz_signal, 'LineWidth', 2.5);
title('Non-Return-to-Zero (NRZ) Polar');
xlabel('Bit Index');
ylabel('Amplitude');

%% NRZ-I
% Initialize the NRZ-I signal
nrz_i_signal = zeros(1, num_bits); % The signal has the same length as the number of bits

% Define the initial level (start with a known state, e.g., low level)
current_level = -1;

% Create the NRZ-I signal
for i = 1:num_bits
    if rand_bits(i) == 1
        % Transition for binary "1"
        current_level = -current_level; % Toggle the level
    end
    % Maintain the same level for binary "0"
    nrz_i_signal(i) = current_level;
end

% Plot the NRZ-I signal
figure;
stairs(1:num_bits, nrz_i_signal, 'LineWidth', 2.5);
title('Non-Return-to-Zero Inverted (NRZ-I)');
xlabel('Bit Index');
ylabel('Amplitude');

%% Initialize Return-to-Zero (RZ) signal
RZ_signal = zeros(1, num_bits * 2); % Double the length for RZ
time = linspace(0, num_bits, num_bits * 2); % Adjust time axis to match RZ_signal length

% Create the RZ signal
for i = 1:num_bits
    % Index into the RZ_signal array
    index = (i - 1) * 2 + 1;  % Find the position in the RZ array
    if rand_bits(i) == 1
        RZ_signal(index) = 1;  % Set the first half of the bit duration
        RZ_signal(index + 1) = 0;  % Return to zero for the second half
    else
        RZ_signal(index) = 0;  % No signal for a "0"
        RZ_signal(index + 1) = 0;  % Remains zero
    end
end

% Plot the Return-to-Zero (RZ) signal
figure;
stairs(time, RZ_signal, 'LineWidth', 2.5);  % Ensure arrays are same length
title('Return-to-Zero (RZ) Signal');
xlabel('Time');
ylabel('Amplitude');


%% Create the AMI signal

% Initialize AMI signal and a variable to track the sign of "1"s
AMI_signal = zeros(1, num_bits); % Starts as all zeros
z = 1; % This variable will determine the alternating sign
for i = 1:num_bits
    if rand_bits(i) == 1 % Check if the current bit is 1
        AMI_signal(i) = z; % Assign the signal with the current sign
        z = -z; % Alternate the sign for the next "1"
    else
        AMI_signal(i) = 0; % If it's a zero, keep it zero
    end
end

% Plot the AMI signal
figure;
stairs(1:num_bits, AMI_signal, 'LineWidth', 2.5);
title('Alternate Mark Inversion (AMI)');
xlabel('Bit Index');
ylabel('Amplitude');


%% create manchester 
% Initialize Manchester signal
% Each bit is represented by two transitions, so the signal length is doubled
manchester_signal = zeros(1, num_bits * 2); % Each bit has two transitions
time = linspace(0, num_bits, num_bits * 2); % Time axis

% Create the Manchester signal
for i = 1:num_bits
    index = (i - 1) * 2 + 1; % Position in the Manchester array
    
    if rand_bits(i) == 0
        % Encode "0" with low-to-high transition
        manchester_signal(index) = 0;     % First half low
        manchester_signal(index + 1) = 1; % Second half high
    else
        % Encode "1" with high-to-low transition
        manchester_signal(index) = 1;     % First half high
        manchester_signal(index + 1) = 0; % Second half low
    end
end

% Plot the Manchester signal
figure;
stairs(time, manchester_signal, 'LineWidth', 2.5);
title('Manchester Code');
xlabel('Time');
ylabel('Amplitude');

%% multi level transmission 3 
% Initialize 3-PAM signal
% Group every two bits to determine the 3-PAM value
num_symbols = floor(num_bits / 2); % Each symbol represents two bits
pam3_signal = zeros(1, num_symbols);
time = linspace(0, num_symbols, num_symbols); % Time axis

% Create the 3-PAM signal
for i = 1:num_symbols
    % Get the bits for this symbol
    bit1 = rand_bits(2*i - 1);
    bit2 = rand_bits(2*i);

    % Determine the PAM-3 level
    if bit1 == 0 && bit2 == 0
        pam3_signal(i) = -1; % Lowest level
    elseif bit1 == 0 && bit2 == 1
        pam3_signal(i) = 0; % Middle level
    elseif bit1 == 1 && bit2 == 0
        pam3_signal(i) = 1; % Highest level
    elseif bit1 == 1 && bit2 == 1
        pam3_signal(i) = 2; % Highest level (if you need this extra level)
    end
end

% Plot the 3-PAM signal
figure;
stairs(time, pam3_signal, 'LineWidth', 2.5);
title('3-Level Pulse Amplitude Modulation (3-PAM)');
xlabel('Time');
ylabel('Amplitude');



%% PSD CALCULATIONS 
A=1;
f = 0:1:100;
Bit_Duration= 1/20;
PSD_NRZ = (A^2)*Bit_Duration * (sinc(f*Bit_Duration)).^2;
PSD_NRZI = ((A^2)*Bit_Duration/2) * ((sinc(f*Bit_Duration)).^2).* (1+cos(pi*f*Bit_Duration));
PSD_RZ= ((A^2)/2) * ((sinc(f*(Bit_Duration/2))).^2);
PSD_AMI=  ((A^2)*Bit_Duration/4) * ((sinc(f*(Bit_Duration/2))).^2) .* (sin (pi *f*Bit_Duration)).^2;
PSD_MANCHESTER = (A^2)*Bit_Duration  * ((sinc(f*(Bit_Duration/2))).^2) .* (sin (pi *f*Bit_Duration/2)).^2;
PSD_MLT3= ((A^2)*Bit_Duration/4) * ((sinc(f*(Bit_Duration/2))).^2) + ((1./ (pi*f *Bit_Duration)).^2).*((sin(2*pi*f*(Bit_Duration))/2).^2)-((A^2)*Bit_Duration/4).*((sinc(f.*Bit_Duration)).^2);


%% plot PSD

figure ;
plot(f, PSD_NRZ, 'b', 'LineWidth', 2.5);
title('Power Spectral Density of non Return-to-Zero (NRZ)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


figure ;
plot(f, PSD_NRZI, 'b', 'LineWidth', 2.5);
title('Power Spectral Density of NRZI');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

figure ;
plot(f, PSD_RZ, 'b', 'LineWidth', 2.5);
title('Power Spectral Density of Return-to-Zero (RZ)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

figure;
plot(f, PSD_AMI, 'r', 'LineWidth', 2.5);
title('Power Spectral Density of Alternate Mark Inversion (AMI)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

figure ;
plot(f, PSD_MANCHESTER, 'g', 'LineWidth', 2.5);
title('Power Spectral Density of Manchester Code');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


figure;
plot(f, PSD_MLT3, 'c', 'LineWidth', 2.5);
title('Power Spectral Density of 3-Level PAM');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

