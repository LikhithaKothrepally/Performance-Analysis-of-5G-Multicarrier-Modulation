symbolsIn = qammod(inpData(:), 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Process all sub-bands together
offset = subbandOffset; 
symbolsInOFDM = [zeros(offset, 1); symbolsIn; ...
                 zeros(numFFT-offset-subbandSize*numSubbands, 1)];
ifftOut = sqrt(numFFT).*ifft(ifftshift(symbolsInOFDM));

% Plot power spectral density (PSD) over all subcarriers
[psd,f] = periodogram(ifftOut, rectwin(length(ifftOut)), numFFT*2, ...
                      1, 'centered'); 
hFig1 = figure; 
plot(f,10*log10(psd)); 
grid on
axis([-0.5 0.5 -100 20]);
xlabel('Normalized frequency'); 
ylabel('PSD (dBW/Hz)')
title(['OFDM, ' num2str(numSubbands*subbandSize) ' Subcarriers'])
set(hFig1, 'Position', figposition([46 50 25 30]));

% Compute peak-to-average-power ratio (PAPR)
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = PAPR2(ifftOut);
disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);