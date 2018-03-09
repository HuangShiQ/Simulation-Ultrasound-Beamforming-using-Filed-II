function [sys,probe,image]= Para()
sys.channel = 128;
sys.c = 1540;
sys.fs = 40e6;
%% probe parameter
fs  = sys.fs;                             %sample frequency
dt = 1/fs;
probe.f0 = 5e6;                     %  Transducer center frequency [Hz]
probe.height = 6.15e-4;              %  Height of element [m]

probe.kerf = 7.7e-6;                       % gap between elements [m]
probe.width     = 1.46e-4;       % Width of element [m]
probe.pitch = probe.kerf+probe.width;                     %  Kerf [m]
probe.focus = [0 0 40]/1000;                       %  Fixed focal point [m]
probe.element = 256;                         %  Number of physical elements
lambda = sys.c/probe.f0;
probe.no_sub_x = round(probe.width/(lambda/8));
probe.no_sub_y = round(probe.height/(lambda/8));

%% impulse response signal
bandwidth = 0.7;        % probe bandwidth
pulse_duration  = 2.5;             % pulse duration [cycles]
t0 = (-1/bandwidth/probe.f0): dt : (1/bandwidth/probe.f0);
impulse_response = gauspuls(t0, probe.f0, bandwidth);
impulse_response = impulse_response-mean(impulse_response); % get rid of DC
probe.impulse_response =  impulse_response;
te = (-pulse_duration/2/probe.f0): dt : (pulse_duration/2/probe.f0);
excitation = square(2*pi*probe.f0*te+pi/2);
probe.excitation = excitation;
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;
figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% Trans Rev
sys.tx_dynamic_focus = 0;
sys.tx_apo = rectwin(sys.channel)';
sys.rx_dynamic_focus = 1;
sys.rx_apo = hamming(sys.channel)';
sys.tx_focus_loc = 0.04;
sys.tx_focus_times = 0;
sys.emit_line = 128;
%% image 
image.width = 15/1000; 

