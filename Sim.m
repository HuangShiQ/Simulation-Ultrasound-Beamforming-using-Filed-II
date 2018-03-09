%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field II Simulation beamforming signal
% Element->256
% Channel->128
% A linear scan of the phantom was done with a 256 element transducer,
% using 128 active elements with a Rect apodization in transmit and hamming receive


% support receive dynamic focus ,fix focus ,transmit dynamic focus ,fix focus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all
clc;
[sys,probe img]= Para();

addpath('Field_II');

field_init(-1);

Th = xdc_linear_array (probe.element, probe.width, probe.height, probe.kerf,...
    probe.no_sub_x,probe.no_sub_y, probe.focus);
Rh = xdc_linear_array (probe.element, probe.width, probe.height, probe.kerf,...
    probe.no_sub_x,probe.no_sub_y, probe.focus);

xdc_impulse (Th, probe.impulse_response);
xdc_impulse (Rh, probe.impulse_response);
xdc_excitation (Th, probe.excitation);

d_x = img.width/(sys.emit_line -1);
%Create phantom
x_pts = [-4.5 0 4.5]*1e-3;
y_pts = zeros(1,length(x_pts));
z_pts = ones(1,length(x_pts))*probe.focus(3);
phantom_amplitudes = ones(length(x_pts),1);


%-- Loop for each line
tic;
for i=1:sys.emit_line
    
    i
    teta = -img.width/2+(i-1)*d_x;
    
    N_pre  = round(teta/(probe.pitch) + probe.element/2 - sys.channel /2);
    N_post = probe.element - N_pre - sys.channel;
    tx.apod_x=[zeros(1,N_pre) sys.tx_apo zeros(1,N_post)];
    rx.apod_x=[zeros(1,N_pre) sys.rx_apo zeros(1,N_post)];
    %         %-- transmit apodization
    xdc_apodization (Th, 0, tx.apod_x);
    
    %         %-- receive apodization
    xdc_apodization (Rh, 0,rx.apod_x);
    %-- position  focus in tx
    
    xdc_center_focus(Th,[teta 0 0]);
    %         -- transmit focus
    if sys.tx_dynamic_focus ==0
        xdc_focus (Th,sys.tx_focus_times,[teta 0 sys.tx_focus_loc]);
    else
        xdc_dynamic_focus (Th, 0, 0,0);
    end
    
    %-- position of coordinate system for calculation of focus in rx
    xdc_center_focus(Rh,[teta 0 0]);
    if sys.rx_dynamic_focus ==1
        %dynamic focus
        xdc_dynamic_focus (Rh, 0, 0,0);
    else
        %fix focus
        xdc_focus(Rh,0,[teta 0 sys.tx_focus_loc]);
    end
    clf
    subplot(211)
    plot(xdc_get(Th,'apo').*xdc_get(Rh,'focus'))
    
    subplot(212)
    plot(xdc_get(Rh,'apo'))
    pause(0.1)
    
    phantom_positions = [x_pts;y_pts;z_pts]';
    [rf_data,t_start] = calc_scat_multi(Th,Rh,phantom_positions,phantom_amplitudes);
    
    times(i) = t_start;
    
    rf{i} = rf_data(:,find(tx.apod_x ~= 0));
    
    
    
end

%%  Read raw data and sum it
for k=1:sys.emit_line
    rft = rf{k};
    
    for l=1:sys.channel
        y_win(l,:) = sys.rx_apo(l)*rft(:,l);
    end
    %   DAS
    final(1:max(size(y_win)),k) = sum(y_win,1);
    t1_array(k) = times(k);
    
    clear y_win
end
%%  Normalization of the image
final = final / max(max(final));

for i = 1:sys.emit_line
    final_t1(fix(t1_array(i)*sys.fs):(fix(t1_array(i)*sys.fs)+size(final,1)-1),i) = ...
        final(:,i);
    
end

%%
%  Do logarithmic compression
dB_range=60;  % Dynamic range for display in dB

disp('Finding the envelope')
log_env=20*log10(abs(hilbert(final_t1)));
log_env=log_env+dB_range;

% Interpolation
ID=20;
[n,m]=size(log_env);
new_env=zeros(n,m*ID);
for i=1:n
    new_env(i,:)=interp(log_env(i,:),ID);
end
[n,m]=size(new_env);

fn=sys.fs/10;

figure;
image(((1:(ID*sys.emit_line -1))*d_x/ID-sys.emit_line *probe.pitch/2),((1:n)/sys.fs)*1540/2,log_env)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray(60))







