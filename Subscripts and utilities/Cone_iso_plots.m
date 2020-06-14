%% Clear variables and close all windows

%The idea here is to get the spectra of the cones and then cross compare
%them with the spectra of the LEDs to generate a multidimensional array of
%excitations to be able to then calculate cone isolating stimuli (being
%able to vary the excitation to only one cone type while leaving the
%stimulation to the other cone types flat (but non-zero)

%% clean up
clearvars
close all force
load('paths.mat')
addpath(genpath(paths(1).main_path))

% define the figure path
fig_path = strcat(paths(1).fig_path,'Setup\');

%define figure counter
f_c = 0;
%define the save variable (if ==1 it saves)
save_val = 0;
%define the cone excitation save variable
saveconeexc_var = 0;
%% Get the spectra of the zebrafish cones

%absorption maxima for the 4 zebrafish cones
l_max = [565 482 411 361];
% l_max = [552 530 426 361];
% l_max = [510 450 400 360];

%define number of cones
cone_num = 4;

%absorption maxima for the rods
% rod  = 502;

%color vector for the plots
colors = [1 0 0;0 1 0;0 0 1;1 0 1];

%Fit values taken from Govardovskii et al. 2000
A = 69.7;
B = 28;
C = -14.9;
D = 0.674;

b = 0.922;
c = 1.104;

%define vector with the wavelengths used for evaluation
wav1 = (200:800)';
%allocate memory for the cone alpha bands
S = zeros(size(wav1,1),cone_num);
%allocate memory for the cone beta bands
Sb = S;

%define the constant for the beta band fits (Govardovskii et al. 2000)
Ab = 0.26;

%for all the cones
for r = 1:cone_num
    
    %Alpha band calculation
    a = 0.8795 + 0.0459*exp((-(l_max(r)-300)^2)/11940);
     
    x = l_max(r)./wav1;
    
    S(:,r) = 1./(exp(A.*(a-x)) + exp(B.*(b-x)) + exp(C.*(c-x)) + D);
    
%     figure(1)
%     hold('on')
%     plot(wav1,S(:,r),strcat(colors(r),'-'))
%     hold('off')
    
    %Beta band calculation (be careful with the low wavelength ones, check Govardovskii et al)
    bb = -40.5 + 0.195*l_max(r);
    lb = 189+ 0.315*l_max(r);
    
    Sb(:,r) = Ab.*exp(-((wav1-lb)./b).^2);
    
%     figure(1)
%     hold('on')
%     plot(wav1,Sb(:,r),strcat(colors(r),'-'))
%     hold('off')
       
end

%combine the alpha and beta bands
S_all = S;% + Sb;

%include the Risner et al. spectral sensitivity functions

%allocate memory for the weights in each spectral sensitivity function
%(taken from Risner et al. 2006). The first dimension is for cones, the
%second one is for excitation vs inhibition and the third one for retinal
%(ERG), tectal or behavioral sensitivity. Keep in mind these correspond to
%the adult, but I'll try them.
spec_weights = zeros(cone_num,2,3); 

%Values taken from Risner 2006
spec_weights(:,1,1) = [0.07 0.15 0.43 1.11];
spec_weights(:,2,1) = [-0.11 -0.28 0 0];

spec_weights(:,1,2) = [0.12 0.12 0.25 1.01];
spec_weights(:,2,2) = [-0.31 -0.91 0 0];

spec_weights(:,1,3) = [0.33 0.93 0.2 0.27];
spec_weights(:,2,3) = [-11.53 -2.65 0 0];

%allocate memory for the recalculated cone spectra
S_recalc = zeros(size(S_all));

%also allocate memory for a dummy spectra variable to calculate the
%inhibiting weights
S_inh = S_all;
S_inh(:,1:3) = S_all(:,2:4);

%define which spectral sensitivity to use. 1 is retinal, 2 is tectal and 3
%is behavioral
spec_sel = 1;

%for all the cones
for cone = 1:cone_num
    %calculate the weighted spectrum
    S_recalc(:,cone) = spec_weights(cone,1,spec_sel).*S_all(:,cone) + ...
        spec_weights(cone,2,spec_sel).*S_inh(:,cone);
    %make all the negative values 0
    S_recalc(S_recalc(:,cone)<0,cone) = 0;
end

% %replace the original spectra with the weighted ones
% S_all = S_recalc;

% normalize the spectra to integrate to 1

%for all the cones
for cone = 1:cone_num
    %divide the cone spectrum by the integrated power to make it add up to
    %1 (so I can modify it arbitrarily later)
    S_all(:,cone) = S_all(:,cone)./sum(S_all(:,cone));
end

%plot the combined cone spectra
f_c = f_c + 1;
figure(f_c)
hold('on')

%for all the cones
for cone = 1:cone_num
%     %get rid of the zeros in the trace (for figure purposes only)
%     non_zero = S_all(:,cone)>0.001;
    non_zero = true(length(S_all(:,cone)),1);
    plot(wav1(non_zero),S_all(non_zero,cone)./max(S_all(non_zero,cone)),'Color',colors(cone,:)...
        ,'Linewidth',3)
end
title('Absorption spectra of the Zebrafish cones','Fontsize',30)
xlabel('Wavelength (nm)','Fontsize',30)
ylabel('Absorption (A.U.)','Fontsize',30)
set(gca,'FontSize',50)
% set(gca,'FontSize',20,'YLim',[0 1.1],'XLim',[200 800])
%% OFF Calculate cone contribution from an arbitrary stimulus (for paper)

% %turn the measurement from umols of photons into mW at each wavelength
% %define Planck's constant, the speed of light and Avogadro's number
% h = 6.64e-31; %in J*s
% c = 3e8;%in m/s
% Na = 6.02e23;%in particles per mol
% %define the wavelength range
% wav_range = 310:370;%in nm, based on Saul's LED spectrum
% %define the stimulus (assuming a more or less square region of the spectrum
% %in Saul's supplemental for their rig)
% arb_stim = zeros(601,1);
% arb_stim(wav_range-199) = (0.25e-3)/1e4; %from umol/m2/s to mmol/cm2/s
% 
% %allocate memory to store the power at each wavelength
% wav_power = zeros(601,1);
% %for all the wavelength ranges
% for wav = 1:size(wav_power,1)
%     %if the point is within the wavelength range (accounting for the fact
%     %that point 1 is wavelength 200 nm)
%     if any(wav == (wav_range-199))
%         %calculate the power density
%         arb_stim(wav) = arb_stim(wav) * Na*h*c;
%     end
% end
% 
% %calculate the contribution to each cone using the cone excitation matrix
% %calculated above
% 
% %allocate memory for the contributions
% cone_cont = zeros(cone_num,1);
% 
% %for all cones
% for cone = 1:cone_num
%     cone_cont(cone) = sum(arb_stim.*S_all(:,cone));
% end
% 
% %show the calculated values
% cone_cont
%% Load raw projector spectra

close all

%define the number of LEDs
led_num  = 4;
% %select the spectra files
% spectra_files = uipickfiles;

%define the path to search
spec_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\20160729_spectra';

%get the files in the selected path
spec_list = dir(spec_path);
%eliminate the . and ..
spec_list = spec_list(3:end);
%get the names from the cell array
spec_names = {spec_list(:).name}';
%for each color, extract the corresponding file name
spec_red = spec_names{~cellfun(@isempty,strfind(spec_names,'red'))};
spec_green = spec_names{~cellfun(@isempty,strfind(spec_names,'green'))};
spec_blue = spec_names{~cellfun(@isempty,strfind(spec_names,'blue'))};
spec_uv = spec_names{~cellfun(@isempty,strfind(spec_names,'UV'))};

%define a cell array with the spectra file names
spec_files = {spec_red,spec_green,spec_blue,spec_uv};

%allocate memory for the spectra (1 is raw, 2 is low passed)
spec_cell = cell(led_num,3);

%allocate a cell structure for the different fits
spectra_m = zeros(length(wav1),led_num);

figure
%for all the spectra
for led = 1:led_num
    subplot(2,2,led)
    %open the target power file
    FID = fopen(fullfile(spec_path,spec_files{led}),'r');
    %read the text from the file
    spec_raw = textscan(FID,'%f %f',...
        'HeaderLines',1,'MultipleDelimsAsOne',1);
    %close the file
    fclose(FID);
    %combine the info into a single wavelength/intensity matrix to put on the cell
    spec_cell{led,1} = [spec_raw{1},spec_raw{2}];
    %plot the data
    plot(spec_cell{led,1}(:,1),spec_cell{led,1}(:,2))
    hold('on')
    %also low pass the data with a moving average
    wind_size = 20;
    spec_cell{led,2} = filtfilt(ones(1,wind_size)/wind_size,1,spec_cell{led,1}(:,2));
    
    plot(spec_cell{led,1}(:,1),spec_cell{led,2})
    
    %now interpolate the desired range for each LED
    spec_cell{led,3} = interp1(spec_cell{led,1}(:,1),spec_cell{led,2},wav1);
    
    plot(wav1,spec_cell{led,3})
    
%     %normalize the spectra to integrate to 1 and save
%     spectra_m(:,led) = spec_cell{led,3}./abs(sum(spec_cell{led,3}));

    temp_model = fit(wav1,spec_cell{led,3},'gauss1');
    
    %evaluate the data using the model
    temp_data = feval(temp_model,wav1);
    
    %save the normalized version
    spectra_m(:,led) = temp_data./sum(temp_data);
    
end



figure
for led = 1:led_num
    plot(wav1,spectra_m(:,led)./max(spectra_m(:,led)))
    hold('on')
end
xlabel('Wavelength (nm)','FontSize',30)
ylabel('Intensity (a.u.)','FontSize',30)
set(gca,'FontSize',50)
%% Create a combined spectra figure

close all

% define the font size
fontsize = 12;
% define the level threshold
level_threshold = 1e-2;
% define the tint factor
tint_factor = 0.5;
% plot the cone spectra
figure

% define the cone colors
cone_colors = [1 0 0;0 1 0;0 0 1;1 0 1];
% define the projector colors
projector_colors = tint_colormap(cone_colors,tint_factor);

% subplot(2,1,1)
%for all the cones
for cone = 1:cone_num
%     %get rid of the zeros in the trace (for figure purposes only)
%     non_zero = S_all(:,cone)>0.001;
    non_zero = true(length(S_all(:,cone)),1);
    % normalize the spetrum
    norm_spectrum = S_all(non_zero,cone)./max(S_all(non_zero,cone));
    % remove values below the threshold
    norm_spectrum(norm_spectrum<level_threshold) = NaN;
    plot(wav1(non_zero),norm_spectrum,'Color',cone_colors(cone,:)...
        ,'Linewidth',3)
    hold on
end
set(gca,'TickLength',[0 0],'FontSize',fontsize,'LineWidth',2)
set(gca,'XLim',[250 750])
box off
ylabel('Absorption/Emission  ','Fontsize',fontsize)
% pbaspect([7,1,1])

% subplot(2,1,2)

% plot the projector spectra
for led = 1:led_num
    % normalize the spectrum
    norm_led = spectra_m(:,led)./max(spectra_m(:,led));
    % remove sub threshold values
    norm_led(norm_led<level_threshold) = NaN;
    plot(wav1(non_zero),norm_led(non_zero),'Color',projector_colors(led,:),'LineWidth',3,'LineStyle','--')
    hold on
    plot(wav1(non_zero),norm_led(non_zero),'Color',projector_colors(led,:),'LineWidth',1,'LineStyle','-')
end
% title('Absorption spectra of the Zebrafish cones','Fontsize',30)
xlabel('Wavelength (nm)','Fontsize',fontsize)
ylabel('Absorption/Emission  ','Fontsize',fontsize)
set(gca,'TickLength',[0 0],'FontSize',fontsize,'LineWidth',2)
set(gca,'XLim',[250 750])
box off
set(gcf,'Color','w')
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5])
pbaspect([2,1,1])

file_path = fullfile(fig_path,'combinedSpectra.png');
export_fig(file_path,'-r600')
%% OFF Get the spectra from the projector (also fits)

% %here I create the Gaussian-based spectra of the projector LEDs based on
% %their maximum and standard deviation so I can use them in later
% %calculations. Since there are only 3 channels I either used the spectrum
% %for the deep blue LED or the actual UV one
% 
% %define the wavelengths to calculate (same as above)
% wav1 = wav1;
% 
% %define the number of LEDs
% led_num = 4;
% 
% %allocate a cell structure for the different fits
% spectra_m = zeros(length(wav1),led_num);
% 
% %define the colors and symbols for the projector spectra
% c_spec = 'kkkk';
% s_spec = 'os*d';
% 
% % %create a new figure (COMMENT TO PUT THE SPECTRA ALONG WITH THE FISH CONE ONES)
% f_c = f_c + 1;
% figure(f_c)
% hold('on')
% 
% %for all three LEDs
% for led = 1:led_num
%     
%     
%     %Parameters for the Dell M109s projector LEDs
% %         switch t
% %             case 1
% %                 s_0 = 15.118;
% %                 x_0 = 625.65;
% %             case 2
% %                 s_0 = 21.7861;
% %                 x_0 = 523.0325;
% %             otherwise
% %                 s_0 = 12.4151;
% %                 x_0 = 455.1870;
% %         end
%     
%     %Parameters for the custom LED source with the deep UV LEDs
%     switch led
%         
%         case 1 %Amber LED
%             s_0 = 33.495;
%             x_0 = 594.113;
%         
%         case 2 %Blue LED
%             s_0 = 12.423;
%             x_0 = 462.43;
%         case 3 %UV1 LED
%             s_0 = 7.6121;
%             x_0 = 396.93;
% 
%         otherwise %HPUV LED
%             s_0 = 5.613;
%             x_0 = 380.86;
%     end
%     
%     %calculate the actual fit
%     spectra_m(:,led) = exp(-(wav1-x_0).^2./(2*s_0^2));
%     
%     %normalize the LED spectra to integrate up to 1
%     spectra_m(:,led) = spectra_m(:,led)./sum(spectra_m(:,led));
% %     %normalize the LED spectra to max 1
% %     spectra_m(:,led) = spectra_m(:,led)./max(spectra_m(:,led));
%     %get a vector to get rid of the zero values
%     non_zero = spectra_m(:,led)>0.0005;
%     %plot the fit
%     plot(wav1(non_zero),spectra_m(non_zero,led),strcat(c_spec(led),s_spec(led)),'LineWidth',2)
%     
% end
% title('Spectra of the projector LEDs','FontSize',30)
% xlabel('Wavelength (nm)','FontSize',30)
% ylabel('Power (A.U.)','FontSize',30)
% % set(gca,'FontSize',20,'YLim',[0 1],'XLim',[200 800])
% 
% % % Plot an external spectrum (like from a filter)
% % %get the path to the spectrum
% % % [spec_n,spec_p] = uigetfile('*.*');
% % spec_p = 'E:\Behavioral data\Matlab\';
% % % spec_n = 'W5089.csv';
% % spec_n = '383DRLP.csv';
% % 
% % %load the spectrum
% % spec_data = csvread(fullfile(spec_p,spec_n),1,0);
% % %normalize the spectrum data
% % spec_data(:,2) = spec_data(:,2)./max(spec_data(:,2));
% % 
% % %now plot the loaded spectrum
% % plot(spec_data(:,1),spec_data(:,2))
%% Calibration data from the power measurements of the projector

%define the number of power levels
pow_num = 256;

%define the LED colors for plots
led_colors = 'rbmk';

%define the path to the power files
% p_power = 'C:\Users\engert lab 13\Documents\Thorlabs\Powermeter';
% p_power = 'E:\Behavioral data\Thorlabs\Powermeter';
% p_power = 'E:\Behavioral data\Matlab\20160512_powermeter';
p_power = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\20160729_powermeter';

%define the power file names
% amber_file = 'PM100_23-Feb-15_11-29_Amber.txt';
% blue_file = 'PM100_23-Feb-15_11-15_Blue.txt';
% uv_file = 'PM100_23-Feb-15_10-34_UV.txt';
% uv365_file = 'PM100_23-Feb-15_10-34_UV.txt';

% amber_file = 'PM100_13-May-16_02-14_red_mod.txt';
% blue_file = 'PM100_13-May-16_02-27_green_mod.txt';
% uv_file = 'PM100_13-May-16_02-20_blue_mod2.txt';
% uv365_file = 'PM100_13-May-16_02-29_uv_mod.txt';

%get the files in the selected path
p_list = dir(p_power);
%eliminate the . and ..
p_list = p_list(3:end);
%get the names from the cell array
p_names = {p_list(:).name}';
%for each color, extract the corresponding file name
p_red = p_names{~cellfun(@isempty,strfind(p_names,'red'))};
p_green = p_names{~cellfun(@isempty,strfind(p_names,'green'))};
p_blue = p_names{~cellfun(@isempty,strfind(p_names,'blue'))};
p_uv = p_names{~cellfun(@isempty,strfind(p_names,'UV'))};

%define a cell array with the spectra file names
power_files = {p_red,p_green,p_blue,p_uv};

%allocate memory for the power measurements
power_read = cell(led_num,1);

%allocate memory for the power fits
power_fits = cell(led_num,1);

%allocate memory to store the values per intensity level
power_vals = zeros(256,led_num);

%open a new figure and set hold to on
f_c = f_c + 1;
figure(f_c)
hold('on')
%for all the spectra
for led = 1:led_num
    %open the target power file
    FID = fopen(fullfile(p_power,power_files{led}),'r');
    %read the text from the file
    power_raw = textscan(FID,'%s %s %s %f %s',...
        'HeaderLines',1,'MultipleDelimsAsOne',1);
    %close the file
    fclose(FID);
    %load the actual measurements on the allocated cell
    power_read{led} = power_raw{4};
    %low pass the curve
    wind_size = 10;
    power_read{led} = filtfilt(ones(1,wind_size)/wind_size,1,power_read{led});
    
%     %only for mod files
%     pow_temp = cellfun(@str2double,power_raw{2});
%     power_read{led} = pow_temp;
    
    %fit the data to a polynomial
%     fit_custom = fittype('(K)/(1 + C*exp(-x*B))^1/n');
%     power_fits{led} = fit((0:length(power_read{led})-1)',power_read{led}...
%         ,'smoothingspline');%,'exclude',[1:5,end-3:end]);
    %define the vector with the x values to fit
    x_vec = (255:-256/length(power_read{led}):0)';
    power_fits{led} = fit(x_vec,power_read{led}...
    ,'smoothingspline');%,'exclude',[1:5,end-3:end]);
    
    %plot the power curve
    plot(x_vec,power_read{led},strcat(led_colors(led),'*'))
    %plot the fits
    plot(power_fits{led},led_colors(led))
%     %show fit values
%     coeffvalues(power_fits{led})

    %for all 256 intensity values
    for pow = 1:pow_num
        %evaluate the fit at the desired values
        power_vals(pow,led) = feval(power_fits{led},pow-1);
    end
end
%% Load the response function of the power meter

%define the filename and path for the responsitivity function file (units
%of nm and A/W)
res_name = 'NOFILTER_S130VC_responsivity.txt';
res_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\PM100D_response';

%open the file
FID = fopen(fullfile(res_path,res_name),'r');
%extract the data
res_raw = textscan(FID,'%f %f','MultipleDelimsAsOne',1);
%close the file
fclose(FID);

%load the responsitivity function into a matrix
res_fun = [res_raw{1} res_raw{2}];

%define target wavelengths for interpolation
tar_wav = wav1;

%interpolate linearly to get all the points of the curve
res_final = interp1q(res_fun(:,1),res_fun(:,2),tar_wav);
%% OFF Calculate the matrix of excitations per cone and LED

% %define the LED colors
% led_col = [255 190 0;0 0 255;255 0 255;0 0 0]./255;
% 
% %allocate memory for the matrix
% cone_exc = zeros(cone_num,led_num,pow_num);
% 
% %allocate memory for a projector power density matrix
% proj_power = zeros(256,led_num);
% 
% 
% %for all the LEDs
% for led = 1:led_num
%     %calculate the "combined", normalized spectrum of the led and sensor
%     %A = f(wav1) function
%     led_sensor = res_final.*spectra_m(:,led);
%     %normalize so it integrates to 1
%     led_sensor = led_sensor./sum(led_sensor);
%     %for all the cones
%     for cone = 1:cone_num
%         %for all the power levels
%         for pow = 1:256
%             
%             %scale the combined spectrum by the detected current
%             led_sensor_scale = power_vals(pow,led).*led_sensor;
%             
%             %calculate the spectrum of the led with power units
%             led_watts = led_sensor_scale./res_final;
%             
%             %multiply by the normalized cone spectrum
%             cone_watts = led_watts.*S_all(:,cone);
%             
%             %calculate the corresponding cone excitation and store
%             cone_exc(cone,led,pow) = sum(cone_watts);
%             %calculate power, but only for the first cone since it's the
%             %same for the other ones.
%             if cone == 1
%                 %calculate the projector power to I value curve per LED
%                 proj_power(pow,led) = sum(led_watts);
%             end
%         end
%     end
%     
% end
% 
% %normalize the excitation matrix
% cone_exc = cone_exc.*255./max(cone_exc(:));
% 
% %plot the projector power for each LED and power level
% f_c = f_c + 1;
% figure(f_c)
% hold('on')
% 
% %for all LEDs
% for led = 1:led_num
%     %plot the integrated power density of the LED at that wavelength
%     plot(proj_power(:,led).*1e6,'LineStyle','none','MarkerFaceColor',led_col(led,:),...
%         'MarkerEdgeColor',[0 0 0],'Marker','o')
% end
% 
% title('LED I value to Power curves','FontSize',20)
% xlabel('I value (A.U.)','FontSize',20)
% ylabel('Power density (uW/cm2)','FontSize',20)
% 
% % % UNCOMMENT to save the projector power file
% % %define the path to save the file and the format as .mat
% % [f_name,p_name] = uiputfile('*.mat');
% % %save the file to the target folder
% % save(fullfile(p_name,f_name),'proj_power')
%% Calculate cone excitations (new version 20160512)

close all

%define the LED colors (for when plotting them)
% led_col = [255 190 0;0 0 255;255 0 255;0 0 0]./255;
led_col = [1 0 0;0 1 0;0 0 1;1 0 1];

%allocate memory for the matrix
cone_exc = zeros(cone_num,led_num,pow_num);

%allocate memory for a projector power density matrix
proj_power = zeros(256,led_num);

figure

%for all the LEDs
for led = 1:led_num 
    
    %for all the power levels
    for pow = 1:256
        
        %spread the current detected across the detector
        led_curr = spectra_m(:,led).*power_vals(pow,led);
        %now calculate the power at each wavelength
        led_pow = led_curr./res_final;
        
        %calculate the corrected projector power to I value curve per LED
        %(since the power_vals vector includes the values measured by the
        %meter, so not the real power)
        proj_power(pow,led) = sum(led_pow);

        %for all the cones
        for cone = 1:cone_num
            
            %calculate the power from each wavelength absorbed by a cone
            power_cone = led_pow.*S_all(:,cone);
            
            %finally calculate the integrated excitation at that cone
            cone_exc(cone,led,pow) = sum(power_cone);
            
            %at max power
            if pow == 256
                %plot the spectrum of each LED in power units
                plot(wav1,led_pow.*1e6,'LineStyle','none','MarkerFaceColor',led_col(led,:),...
                    'MarkerEdgeColor',[0 0 0],'Marker','o')
                ylabel('Power density (uW/cm2)','FontSize',20)
                xlabel('Wavelength (nm)','FontSize',20)
                hold('on')
            end
        end
    end
    
end

%normalize the excitation matrix to 8 bits
cone_exc = cone_exc.*255./max(cone_exc(:));

%plot the projector power for each LED and power level
% f_c = f_c + 1;
% figure(f_c)
figure
hold('on')

%for all LEDs
for led = 1:led_num
    %plot the integrated power density of the LED at that wavelength
    plot(proj_power(:,led).*1e6,'LineStyle','none','MarkerFaceColor',led_col(led,:),...
        'MarkerEdgeColor',[0 0 0],'Marker','o')
    hold('on')
%     plot(power_vals(:,led).*1e6)
end
set(gca,'FontSize',15,'TickLength',[0 0])

% title('LED I value to Power curves','FontSize',20)
xlabel('Projector pixel value (A.U.)','FontSize',20)
ylabel('Power density (uW/cm2)','FontSize',20)

box off
set(gcf,'Color','w')
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5])
% pbaspect([2,1,1])

file_path = fullfile(fig_path,'projectorPower.png');
export_fig(file_path,'-r600')
%% Stacked bar plot with the excitations for each cone from each LED

close all
figure

%for max power
t_power = 256;

%allocate memory to store the excitations
cone_bar = zeros(cone_num,led_num);
%for all the cones
for cone = 1:cone_num
    %store the corresponding excitation values
    cone_bar(cone,:) = cone_exc(cone,:,t_power);
end

%normalize the matrix to the sum of the highest activated cone
cone_bar = cone_bar./max(sum(cone_bar,2));

if led_num == 4
    cone_labels = {'Red cone','Green cone','Blue cone','UV cone'};
else
    cone_labels = {'Red cone','Blue cone','Green cone'};
end
%plot the excitations in a stacked bar plot
b = bar(cone_bar,'stacked');
set(gca,'XTickLabel',cone_labels,'FontSize',20,'YLim',[0 1.2],'XTickLabelRotation',45)
set(b(1),'FaceColor',[1 0 0])
set(b(2),'FaceColor',[0 1 0])
set(b(3),'FaceColor',[0 0 1])
if led_num == 4
    set(b(4),'FaceColor',[1 0 1])
end
ylabel('Relative cone excitation (a.u.)','Fontsize',20)
error('Stop here')
%% Calculate cone isolating stimuli (based on Estevez and Spekreijse)

%the variables I need are S_all, containing all the cone spectra normalized
%to integrate to 1 and also spectra_m, which contains all the LED spectra,
%also normalized to integrate to 1. The idea is for this code to give me
%the ratios in total power between the LEDs to get cone isolating stimuli
%and then I can input the actual powers to get those total ratios

close all
%calculate the cone excitation matrix to be inverted

%define the starting I levels for each one of the LEDs
start_I = [127 127 127 127];
%calculate the power corresponding to the specified I values
start_power = zeros(size(start_I));
%for all the LEDs
for leds = 1:led_num
    %get the power corresponding to each I value
    start_power(leds) = proj_power(start_I(leds),leds).*1e9;
end
%allocate memory for the cone excitation matrix (where the dimensions
%should be 4 cones by 4 LEDs)
cone_excm = zeros(cone_num,led_num);

%for all the cones
for cone = 1:cone_num
    %for all the LEDs
    for led = 1:led_num
        %calculate the excitation at that point (in units of power in nW)
        cone_excm(cone,led) = sum(S_all(:,cone).*spectra_m(:,led).*proj_power(start_I(led),led).*1e9);
    end
end

% %calculate the inverse of the excitation matrix
% cone_inv = inv(cone_excm);
figure
imagesc(cone_excm)
xlabel('Cone')
ylabel('LED')
title(strcat('Excitation matrix for LED profile:',num2str(start_I)))
set(gca,'XTick',1:cone_num,'YTick',1:led_num)
%calculate the ratios of the LEDs for each cone

%define the number of excitation levels to be calculated
num_lev = 40;
%define the number of cycles
wave_cyc = 5;

%get the vector to calculate the sine with
wave_vec = 2*pi()*wave_cyc*linspace(0,1,num_lev);
%finally, calculate the sine
wave_sine = sin(wave_vec);


% exc_lev = [linspace(-6,6,num_lev);linspace(-5,5,num_lev);...
%     linspace(-1,1,num_lev);linspace(-0.5,0.5,num_lev)];

%allocate memory for the LED ratios (one dimension will contain the different
%levels of excitation the next one the cones and the last one the LEDs,)
cone_iso = zeros(num_lev,cone_num,led_num);

%allocate memory for the I level output matrix
I_mat = zeros(size(cone_iso));

%allocate memory to store the dynamic range values
range_val = zeros(cone_num,1);

%for each cone
for cone = 1:cone_num
    %for a range of values
    for vals = 5:-0.01:0.01
        %initialize a counter for the levels
        lev_count = 1;
        %for each level of excitation
%         for levs = linspace(-vals,vals,num_lev)
        for levs = vals.*wave_sine

            %     for levs = exc_lev(cone,:)
            %define the vector of desired modulation
            mod_vec = zeros(cone_num,1);
            mod_vec(cone) = levs;
            %calculate the ratios of the LEDs
            cone_iso(lev_count,cone,:) = cone_excm\mod_vec;
            %update the value of the counter
            lev_count = lev_count + 1;
        end
        
        
        
        %Calculate the I levels corresponding to the modified power levels
        
        
        %for all the LEDs
        for leds = 1:led_num
            %load the part of the matrix corresponding to this LED and add it to
            %the starting power value
            temp_mat = cone_iso(:,cone,leds) + start_power(leds);
            %interpolate and round the I values from the power
            I_mat(:,cone,leds) = interp1(proj_power(:,leds).*1e9,(0:255)',temp_mat,'nearest');
        end
        %when the NaNs disappear
        if all(~isnan(I_mat(:)))
            %store the range value
            range_val(cone) = vals;
            %stop the loop and go to the next cone
            break
        end
    end
end
figure
for cone = 1:cone_num
    subplot(2,2,cone)
    imagesc(squeeze(cone_iso(:,cone,:)))
    xlabel('Cone')
    ylabel('Frames')
    title('Cone isolation matrix')
end
figure
for cone = 1:4
    subplot(2,2,cone)
    imagesc(squeeze(I_mat(:,cone,:)))
    xlabel('Cone')
    ylabel('Frames')
    title('LED isolation I values')
end
%% Add the second side

% %Phototaxis (1 side modulated)
% dark_mat = zeros(size(I_mat));
% %dark on side 2
% Out_mat1 = cat(3,I_mat,dark_mat);
% 
% %dark on side 1
% Out_mat2 = cat(3,dark_mat,I_mat);
% 
% Out_mat = cat(2,Out_mat1,Out_mat2);

%Full field (modulated)
Out_mat = cat(3,I_mat,I_mat);
%% OFF Save the file with the I values

% %if the save variable is on
% if save_val == 1
%     %define the main saving path
%     I_path = 'E:\Behavioral data\Matlab\Ivals\20160817_Ivals';
%     %for all the cones
%     for cone = 1:size(Out_mat,2)
%         %define the file name including the time stamp
%         I_name = strcat(datestr(now,30),'_Ivals_',num2str(cone),'.txt');
%         %open the file for writing
%         FID = fopen(fullfile(I_path,I_name),'w');
%         %write the text file with the plane corresponding to the cone
%         fprintf(FID,'%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\r\n',squeeze(Out_mat(:,cone,:))');
%         %close the file
%         fclose(FID);
%     end
% end
%% Save the cone excitation matrix

%if the saveconeexc_var is active
if saveconeexc_var == 1
    %define the save path
    cexc_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\201609019_coneexc';
    %define the file name
    cexc_name = strcat(datestr(now,30),'_ConeExc.mat');
    %save the file
    save(fullfile(cexc_path,cexc_name),'cone_exc')
end
%% OFF Plot the cone isolation levels?

% f_c = f_c + 1;
% figure(f_c)
% hold('on')
% %for all the cones
% for cone = 1:cone_num
%     
%     
%     %for all the excitation levels
%     for levs = 1:exc_lev
%         uv_color = [cone_iso(levs,cone,4)./255,cone_iso(levs,cone,4)./255,...
%             cone_iso(levs,cone,4)./255];
%         plot3(cone_iso(levs,cone,1),cone_iso(levs,cone,2),cone_iso(levs,cone,3)...
%             ,'MarkerFaceColor',uv_color,'Marker','o',...
%             'MarkerEdgeColor',uv_color)
%     end
% end
%% OFF Rearrange the excitation matrix using the excitations as coordinates

% %define the max value of the range for excitation normalization (which will
% %in turn define the resolution of the normalization)
% norm_range =  255;
% 
% %take the log of the cone excitation matrix to let the values be within the
% %same order of magnitude. Also normalize it to 0-norm_range
% cone_log = log(cone_exc);
% cone_log = round(norm_range.*(cone_log-min(cone_log(:)))./max((cone_log(:))-min(cone_log(:))));
% 
% %allocate memory for the power levels
% cone_power = zeros(size(cone_log));
% 
% %allocate memory to store the repeating values when present
% cone_reps = cell(size(cone_log));
% 
% %for all the cones
% for cone = 1:cone_num
%     %for all the leds
%     for led = 1:led_num
%         %for all the excitation levels
%         for exc = 1:norm_range+1
%             %get the power level associated with the given condition
%             pow_level = find(cone_log(cone,led,:)==exc-1);
%             %if the power level is empty
%             if isempty(pow_level)
%                 %store a zero in the cell
%                 cone_reps{cone,led,exc} = 0;
%             else %if not
%                 %store the value(s) in the cell
%                 cone_reps{cone,led,exc} = pow_level;
%             end
%             %get the average of the cell value and store it in the matrix
%             cone_power(cone,led,exc) = mean(cone_reps{cone,led,exc});
%         end
%     end
% end
% %% Plot the results of the cone isolation
% 
% % cell_size = cellfun(@length,cone_reps,'UniformOutput',0);
% % sum([cell_size{:}]~=0)
% 
% %for all cones
% for cone = 1:cone_num
%     f_c = f_c + 1;
%     figure(f_c)
% %     subplot(2,2,cone)
%     mesh(squeeze(cone_power(cone,:,:)))
%     title(['Cone ',num2str(cone),' excitation'],'FontSize',20)
%     xlabel('Excitation','FontSize',20)
%     ylabel('LED','FontSize',20)
%     zlabel('Power level','FontSize',20)
% end
% 
% %for all LEDs
% for led = 1:led_num
%     f_c = f_c + 1;
%     figure(f_c)
% %     subplot(2,2,cone)
%     mesh(squeeze(cone_power(:,led,:)))
%     title(['LED ',num2str(led),' excitation'],'FontSize',20)
%     xlabel('Excitation','FontSize',20)
%     ylabel('Cone','FontSize',20)
%     zlabel('Power level','FontSize',20)
% end
% %% Attempts to plot the data
% 
% %plot the LED power level as a function of excitation
% f_c = f_c + 1;
% figure(f_c)
% hold('on')
% 
% %define the cone colors
% cone_color = 'rgbm';
% 
% %define the led symbols
% led_sym = '*so';
% 
% %for all the leds
% for led = 1:led_num
%     %for all the cones
%     for cones = 1:cone_num
%         plot_tar = squeeze(cone_power(cones,led,:));
%         plot_tar(plot_tar==0) = NaN;
%         plot(plot_tar,strcat(led_sym(led),cone_color(cones)))
%     end
% end
% 
% title('Power level for each LED and cone given an excitation','FontSize',20)
% xlabel('Excitation (A.U.)','FontSize',20)
% ylabel('LED power level (A.U.)','FontSize',20)
% set(gca,'Xgrid','on')
% 
% %plot the excitation as a function of LED power
% f_c = f_c + 1;
% figure(f_c)
% hold('on')
% 
% %define the cone colors
% cone_color = 'rgbm';
% 
% %define the led symbols
% led_sym = '*so';
% 
% %for all the leds
% for led = 1:led_num
%     %for all the cones
%     for cones = 1:cone_num
%         plot(squeeze(cone_exc(cones,led,:)),strcat(led_sym(led),cone_color(cones)))
%     end
% end
% 
% title('Excitation for each LED and cone given a power level','FontSize',20)
% ylabel('Excitation (A.U.)','FontSize',20)
% xlabel('LED power level (A.U.)','FontSize',20)
% set(gca,'Xgrid','on')
%% OFF Save the excitation file as a .mat file
% %define the path to save the file and the format as .mat
% [f_name,p_name] = uiputfile('*.mat');
% %save the file to the target folder
% save(fullfile(p_name,f_name),'cone_exc')