

function peakVec = plotExperimentalData(filename,type)
    % This function takes in the accelerometer data from the shaker table 
    % Plots the signal responece of each of the accelerometrs in both the time
    % and frequency domain
    
    % Load experimental data file
%       cd("Data/")
        data = readmatrix(filename);
%        cd("../")

    % Asign data to variables (only looking at the test_2min_all_8)
    totalTest = data;

    time = totalTest(:,1); 

   % Plot accelerometer signal in time and frequency domain
   name = ["Accelerometer 01" "Accelerometer 02" "Accelerometer 03"];
    j = 1; % name indexing variable
    for i = 3:5
        accel = totalTest(:,i);
        Fs = 1/(sum(diff(time))/length(diff(time))); % sampling frequency
        t = time; % time vector
        N = length(accel); % length of signal
        xdft = fft(accel); % calculate FFT of signal
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/N:Fs/2;

%         % Plot in frequency domain
%         figure
%         hold on
%         grid on
%         % since this psdx already squared the fft (xdft), it is in power units; to
%         % plot in dB, take 10*log10:
%         %plot(freq,10*log10(psdx))
%         plot(freq,psdx)
%         title(name(j)+' Periodogram Using FFT ' + "("+type+")")
%         xlabel('Frequency (Hz)')
%         ylabel("Amplitude") %'Power/Frequency dB(Vrms^2/Hz)'
%         xlim([0 55])
%         hold off
        peak = max(psdx);
        peakIndex = find(psdx == peak);
        peakVec(j) = freq(peakIndex);
        j = j+1;

    end

end