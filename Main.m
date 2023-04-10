

recieve(channel(tansmitt()));

%% Transmitter Side
function s = tansmitt()
    prompt = {'Station 1 Frequency (Hz)(Total Available Bandwith 1Khz-10)','Station 1 sound content (file name)','Station 1 bandwith (Hz)','Station 2 Frequency (Hz)(Total Available Bandwith 1Khz-10)','Station 2 sound content (file name)','Station 2 bandwith (Hz)(Do not overlap with station 1!!!)'};

    dlgtitle = 'Radio Transmission Station';
    dims = [2 250];
    definput = {'3000','1.wav','4000','8000','2.wav','4000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    %answer{1}
    % Load the audio signals
    [y1, Fs1] = audioread(answer{2});
    [y2, Fs2] = audioread(answer{5});
    
    % Resample if necessary to a common sampling frequency
    if Fs1 ~= Fs2
        if Fs1 > Fs2
            y2 = resample(y2, Fs1, Fs2);
            Fs2 = Fs1;
        else
            y1 = resample(y1, Fs2, Fs1);
            Fs1 = Fs2;
        end
    end

    % pad the smaller array with zeros-
    if length(y1) > length(y2)
        y2 = [y2.', zeros(1, length(y1) - length(y2))];
        y2=y2.';
    else
        y1 = [y1.', zeros(1, length(y2) - length(y1))];
        y1=y1.';
    end
    
  
   
    
    y1=lowpass(y1,str2double(answer{3})/4,Fs1);
    y2=lowpass(y2,str2double(answer{6})/4,Fs1);
    
    %zeros(1, length(y1) - length(y2))
    % Normalize the audio signals to have maximum amplitude of 1
    y1 = y1 / max(abs(y1));
    y2 = y2 / max(abs(y2));
    
    % Set the carrier frequency and sampling frequency
    fc1 = str2double(answer{1}); % Carrier frequency
    fc2 = str2double(answer{4});
    Fcs1 =  fc1; % Sampling frequency
    Fcs2 = fc2; 
    % Create the carrier signal
    t1 = linspace(0, length(y1)/Fs1, length(y1));
    t2 = linspace(0, length(y2)/Fs1, length(y2));
    carrier1 = cos(2*pi*fc1*t1);
    carrier2 = cos(2*pi*fc2*t2);
    
    % Modulate the audio signals onto the carrier signal
    s1 = ((y1) .* carrier1');
    s2 = ((y2) .* carrier2');
    
    % Combine the modulated signals
    %s = s1 + s2;
    s = s1 + s2;
    
    % Plot the time-domain signals
    figure;
    subplot(3,2,1);
    plot(y1);
    title('Station 1 Time Domain');
    subplot(3,2,3);
    plot(y2);
    title('Station 2 Time Domain');
    subplot(3,2,5);
    plot(s);
    title('Combained Time Domain');
    
    % Plot the frequency-domain signals
    N = length(s);
    N1= length(y1)
    f1 = (-N1/2:N1/2-1) * Fs1/N1;
    f = (-N/2:N/2-1) * Fs1/N;
    S = fftshift(abs(fft(s))/N);
  
    subplot(3,2,2);
    plot(f1, abs(fftshift(fft(y1))/N1));
    title('Station 1 Frequency Domain');
    subplot(3,2,4);
    plot(f1, abs(fftshift(fft(y2))/N1));
    title('Station 2 Frequency Domain');
    subplot(3,2,6);
    plot(f, S);
    ylim([0 0.003])
    xlim([0 11025])
    title('Modulated Signal (Frequency Domain)');
    prompt = {'Modulated signal is ready! Please verify all the plots before proceeding to channel.'};
    dlgtitle = 'Complete!';
    dims = [0 250];
    definput = {'4000',''};
    r = inputdlg(prompt,dlgtitle,dims,definput);
    %f = msgbox("Modulated signal is ready! Please verify all the plots before proceeding to channel.");


end



%% channel 
function output= channel(s)
    prompt = {'Enter SNR value'};

    dlgtitle = 'Welcome to channel service provider';
    dims = [3 250];
    definput = {'34'};
    an2 = inputdlg(prompt,dlgtitle,dims,definput);
    output=awgn(s,str2double(an2{1}),'measured');

    prompt = {'Signal transmitter through channel, tune using radio to listrn'};
    dlgtitle = 'Complete!';
    dims = [0 250];
    definput = {'4000',''};
    r = inputdlg(prompt,dlgtitle,dims,definput);
        
end
%% Receiver Side
function recieve(s)
    prompt = {'Bandwith Of Radio','Frequency to Tune to..'};
    dlgtitle = 'Welcome to Radio!';
    dims = [2 250];
    definput = {'4000','3000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Fs1=22050;
    t1 = linspace(0, length(s)/Fs1, length(s));

    carrier1 = cos(2*pi*3000*t1);
    carrier2 = cos(2*pi*8000*t1);
    %yr1=s;
    figure;
    yr1=s.*cos(2*pi*str2double(answer{2})*t1)';
    subplot(3,1,1);
    plot(yr1);
    %str2double(answer{1})/2
    gauss_filter=fspecial('gaussian',[1 length(s)],10);
    yr1=lowpass(yr1,str2double(answer{1})/2,Fs1);
    yr1= imfilter(yr1',gauss_filter);
    subplot(3,1,2);
    plot(yr1);
    %yr1=yr1-1;
    subplot(3,1,3);
    plot(yr1);
    title('Time Domain Output');
    f = msgbox("Demodulation complete! Sound playing!");
    sound(yr1*2,Fs1)        
end











