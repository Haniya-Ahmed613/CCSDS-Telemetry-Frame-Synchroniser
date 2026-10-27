
clear;
clc;
close all;

rng(7, 'twister');  % Repeatable random results

%% 1. PROJECT PARAMETERS

SYNC_HEX = '1ACFFC1D';
markerBits = double(dec2bin(hex2dec(SYNC_HEX), 32) - '0');

MARKER_LEN = 32;

% Frame length excludes the 32-bit sync marker.
FRAME_LEN_BITS = 1784;

NUM_FRAMES = 500;

% Desired theoretical BER values for uncoded BPSK over AWGN.
targetBERs = [0, 1e-4, 1e-3, 1e-2, 0.05, 0.10, 0.15];

% Maximum number of bit mismatches allowed in a marker window.
thresholds = [0, 2, 4];

fprintf('CCSDS Telemetry Frame Synchronizer - Level 1\n');
fprintf('Sync marker: %s\n', SYNC_HEX);
fprintf('Frame length excluding ASM: %d bits\n', FRAME_LEN_BITS);
fprintf('Number of frames: %d\n\n', NUM_FRAMES);

%% 2. GENERATE THE TRANSMITTED BITSTREAM

% Each frame is preceded by a 32-bit ASM.
% The remaining FRAME_LEN_BITS are random demonstration data.

framePeriod = MARKER_LEN + FRAME_LEN_BITS;
streamLength = NUM_FRAMES * framePeriod;

txBits = zeros(1, streamLength);

% Starting index of each true sync marker.
trueStarts = 1 + (0:NUM_FRAMES-1) * framePeriod;

for n = 1:NUM_FRAMES

    pos = trueStarts(n);

    % Insert the marker, MSB first.
    txBits(pos:pos+MARKER_LEN-1) = markerBits;

    % Fill the Transfer Frame with random demonstration bits.
    dataStart = pos + MARKER_LEN;
    dataEnd = dataStart + FRAME_LEN_BITS - 1;

    txBits(dataStart:dataEnd) = randi([0 1], 1, FRAME_LEN_BITS);

end

fprintf('Generated %d transmitted bits.\n\n', streamLength);

%% 3. BPSK MODULATION

% Bit 0 becomes +1.
% Bit 1 becomes -1.

txSymbols = 1 - 2 * txBits;

%% 4. PREPARE THE RESULT ARRAYS

numBER = numel(targetBERs);
numThresholds = numel(thresholds);

measuredBER = zeros(numBER, 1);

detectionRate = zeros(numBER, numThresholds);
missRate = zeros(numBER, numThresholds);

falseCandidateCount = zeros(numBER, numThresholds);
falseAlarmRate = zeros(numBER, numThresholds);

% Every possible 32-bit window is tested.
numWindows = streamLength - MARKER_LEN + 1;

% Exclude the known correct marker starting positions.
numNonMarkerWindows = numWindows - NUM_FRAMES;

%% 5. TEST DIFFERENT AWGN CONDITIONS

for ib = 1:numBER

    p = targetBERs(ib);

    % Convert the desired theoretical BPSK BER into noise strength.
    % Eb = 1 for our unit-amplitude BPSK symbols.
    
    if p == 0

        noiseSigma = 0;

    else

        EbN0 = (erfcinv(2*p))^2;
        noiseSigma = sqrt(1/(2*EbN0));

    end

    % Add real Gaussian noise to the transmitted BPSK samples.
    rxSymbols = txSymbols + ...
        noiseSigma * randn(size(txSymbols));

    % Hard-decision receiver: recover bits from noisy amplitudes.
    rxBits = double(rxSymbols < 0);

    % Measure the actual BER in this trial.
    measuredBER(ib) = mean(rxBits ~= txBits);

    %% 6. SLIDING-WINDOW CORRELATION

    % A marker bit equal to 0 has weight +1.
    % A marker bit equal to 1 has weight -1.
    %
    % This convolution calculates the mismatch count for
    % every possible 32-bit window efficiently.

    weights = 1 - 2*markerBits;

    mismatchCount = sum(markerBits) + ...
        conv(rxBits, fliplr(weights), 'valid');

    mismatchCount = round(mismatchCount);

    %% 7. APPLY DIFFERENT DETECTION THRESHOLDS

    for it = 1:numThresholds

        T = thresholds(it);

        % Accept a candidate when mismatches are <= T.
        candidateStarts = find(mismatchCount <= T);

        % Did each true marker get detected at its exact position?
        detectedFrames = ismember(trueStarts, candidateStarts);

        detectionRate(ib,it) = mean(detectedFrames);
        missRate(ib,it) = 1 - detectionRate(ib,it);

        % Candidate detections at incorrect bit positions.
        isFalse = ~ismember(candidateStarts, trueStarts);

        falseCandidateCount(ib,it) = sum(isFalse);

        % Fraction of non-marker windows falsely accepted.
        falseAlarmRate(ib,it) = ...
            falseCandidateCount(ib,it) / numNonMarkerWindows;

    end

    fprintf('Target BER: %.4g | Measured BER: %.6f\n', ...
        p, measuredBER(ib));

end

%% 8. DISPLAY RESULTS TABLE

numRows = numBER * numThresholds;

TargetBER = zeros(numRows,1);
MeasuredBER = zeros(numRows,1);
ThresholdBits = zeros(numRows,1);
DetectionProbability = zeros(numRows,1);
MissProbability = zeros(numRows,1);
FalseCandidates = zeros(numRows,1);
FalseAlarmProbability = zeros(numRows,1);

row = 0;

for ib = 1:numBER
    for it = 1:numThresholds

        row = row + 1;

        TargetBER(row) = targetBERs(ib);
        MeasuredBER(row) = measuredBER(ib);
        ThresholdBits(row) = thresholds(it);

        DetectionProbability(row) = detectionRate(ib,it);
        MissProbability(row) = missRate(ib,it);

        FalseCandidates(row) = falseCandidateCount(ib,it);
        FalseAlarmProbability(row) = falseAlarmRate(ib,it);

    end
end

results = table(TargetBER, MeasuredBER, ThresholdBits, ...
    DetectionProbability, MissProbability, ...
    FalseCandidates, FalseAlarmProbability);

disp(' ');
disp('========== FINAL RESULTS ==========');
disp(results);

%% 9. PLOT DETECTION PROBABILITY

% BER = 0 is omitted from logarithmic x-axis plots.
% Its noise-free result is still included in the table.

figure('Name', 'Detection Probability');

hold on;

for it = 1:numThresholds

    semilogx(targetBERs(2:end), ...
        detectionRate(2:end,it), '-o', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Threshold T = %d', thresholds(it)));

end

grid on;
xlabel('Target BER of BPSK/AWGN Channel');
ylabel('Detection Probability');
title('Sync Marker Detection vs BER');
ylim([0 1.02]);
legend('Location', 'best');

%% 10. PLOT FALSE ALARM PROBABILITY

figure('Name', 'False Alarm Probability');

hold on;

for it = 1:numThresholds

    semilogx(targetBERs(2:end), ...
        falseAlarmRate(2:end,it), '-o', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Threshold T = %d', thresholds(it)));

end

grid on;
xlabel('Target BER of BPSK/AWGN Channel');
ylabel('False Alarm Probability per Non-Marker Window');
title('False Alarm Probability vs BER');
legend('Location', 'best');

fprintf('\nSimulation completed.\n');
