% MASS LOSS RATE CONTROL OF THE FPA HEATING - RORY HADDEN August 2022

%Use at your own risk
clear
close all
clc
instrreset

load ramp_calibrations.mat

% desired_MLR = inputdlg('What mass loss rate do you want (g/s/m^2)?');
% sample_area = 0.96*0.96;

prompt = {'Desired MLR, g/s/m^2', 'Sample area, m^2', 'Material'};
dlgtitle = 'Some details of the test';
dims = [1 100];
definput = {num2str(2), num2str(0.096*0.096), 'Test'};
test_info = inputdlg(prompt,dlgtitle,dims,definput);
desired_MLR = str2double(cell2mat(test_info(1)));
sample_area = str2double(cell2mat(test_info(2)));
test_ID = test_info(3);

DL = visa('Agilent','GPIB0::9::INSTR');
loadcell = tcpip('192.168.127.254', 4001);

% connect to the data logger and the load cell
fopen(DL);
fopen(loadcell);

%get some initial data on teh sample mass (note that that split function
%separates the text string
%zero_loadcell = query(loadcell, 'Z');
mass_query = query(loadcell, 'SI');
mass_query = split(mass_query);
initial_mass = str2double(cell2mat(mass_query(3)));

%do some house keeping on the data logger
fprintf(DL,'ABOR');     % fprintf:  write data to instrument DL
fprintf(DL,'*RST');     % Reset the instrument to the Factory configuration.

%tell it what channels you care about
DL_CHANNELS = '110';%'101:116';
DL_CHANNELS_VDC_OUTPUT = '304'; % THIS IS THE CHANNEL THAT 'TAKES' THE INPUT VOLTAGE

%the calibration factor for the HFG - not needed for constant MLR?
HFG_conversion = 9.833*1000; %kW/m2 per V;

% define the format of the time
TimeStampFormat = 'yyyymmddHHMMSS.FFF';

%make a filename that avoids overwritign data
FileName = [datestr(now, TimeStampFormat) '_' num2str(desired_MLR) '_gsm2_' cell2mat(test_ID)];

% define the PID constants
PID_kp = 0.2 %0.2;
PID_ki = 0.2%0.04;
PID_kd = 0.2%0.2;

% PROTECT THE FPA! (lamps)
max_lamp_voltage = 4.5;
min_lamp_voltage = 0.25;

%define a waiting period to allow calculation of the MLR
pretest_period = 60;
averaging_window = 30;
ave_frame = 10;

%a counting variable
i=1;
% define initial value for integral_term - maybe not necessary
integral_term = 0;
%define initial output as 0
output = 0;
VDC_out = 0;
%define an irradiation rate
irradiation_rate = 0.025;
VDC_init = 2;

%define the variable setpoint
setpoint = desired_MLR;

% reset the voltage to the lamps to 0
fprintf(DL,['SOUR:VOLT 0, (@' DL_CHANNELS_VDC_OUTPUT ')']);
%tell the DL the channels to scan
fprintf(DL,['ROUT:SCAN (@' DL_CHANNELS ')']);

%make a figure
test_fig = figure;
mass_axes = axes('Parent', test_fig, 'Position', [0.13 0.61 0.775 0.35], 'box', 'on');
HFG_axes = axes('Parent', test_fig, 'Position', [0.13 0.11 0.775 0.35], 'box', 'on');

HFG_axes.YLabel.String = 'VDC out, V';
HFG_axes.XLabel.String = 'Test time, s';
HFG_axes.YLim = [0 5];
HFG_axes.XLim = [0 1000];

mass_axes.YLabel.String = 'Mass loss rate, g/s/m2';
mass_axes.XLabel.String = 'Test time, s';
mass_axes.XLim = [0 1000];
mass_axes.YLim = [-setpoint*2 setpoint*2];

hold(HFG_axes)
hold(mass_axes)

HF_plot = plot(0, 0, 'x', 'parent', HFG_axes, 'color', 'r');
mass_plot = plot(0, 0, 'o', 'parent', mass_axes, 'color', 'b');
mass_plot = plot([0 1000], [setpoint setpoint], '-', 'color', 'k');
drawnow

stop_button = uicontrol('Style', 'Pushbutton', 'String', 'Stop', 'fontsize', 12, ...
    'BackgroundColor',[0.8 0.8 0.8], 'units', 'normalized', 'Position', [0.9 0 0.1 0.07], ...
    'Callback', 'delete(gcbf)');

% VDC = [-fit_data(2)+sqrt(-fit_data(2)^2 - 4*fit_data(1)*(fit_data(3)-heatflux), ]

while(ishandle(stop_button))

    if i <= pretest_period

        tic

        ScanDataAsText = query(DL,'READ?');
        ScanDataAsNum_Raw = str2num(ScanDataAsText);

        HF_VDC(i) = -ScanDataAsNum_Raw(1);
        HF_eng(i) = -ScanDataAsNum_Raw(1).*HFG_conversion;

        mass_query = query(loadcell, 'SI');
        mass_query = split(mass_query);
        mass(i) = str2double(cell2mat(mass_query(3)));
        time_stamp(i) = str2double(datestr(now, TimeStampFormat));
        test_time(i) = (datenum(num2str(time_stamp(i), '%.3f'), TimeStampFormat)-datenum(num2str(time_stamp(1), '%.3f'), TimeStampFormat)).*(24*60*60);

        MLR = [0 (-diff(mass)./diff(test_time))];
        MLR = MLR./sample_area;

        if i > averaging_window
            smooth_MLR(i) = mean(MLR(i-averaging_window:i));
        else
            smooth_MLR(i) = 0;
        end

        VDC_out(i) = i*irradiation_rate; %VDC_init
        proportional_term(i) = 0;
        integral_term(i) = 0;
        derivative_term(i) = 0;

        fprintf(DL,['SOUR:VOLT ' num2str(VDC_out(i)), ', (@' DL_CHANNELS_VDC_OUTPUT ')']);         % Set the output voltage level on the specified DAC channel. You can set the output voltage to any value between +12 Vdc and -12 Vdc, in 1 mV steps. Each DAC channel is capable of 10 mA maximum output current. The DAC channels are numbered “s04” and “s05”, where s represents the slot number. The :VOLT? query returns the output voltage level on the specified DAC channel. Returns a number in the form “+8.00000000E+00”.

        disp(['VDC to lamps = ' num2str(VDC_out(i)) ' V']);
        disp(['Mass loss rate = ' num2str(MLR(i)) ' g/s/m2']);

        HF_plot = plot(test_time(i), VDC_out(i), 'x', 'parent', HFG_axes, 'color', 'r');
        mass_plot = plot(test_time(i), MLR(i), 'o', 'parent', mass_axes, 'color', [0.5 0.5 0.5]);

        i=i+1;

        toc
        pause(1-toc)

    end

    if i > pretest_period

        tic

        % this bit of the loop gets the data
        ScanDataAsText = query(DL,'READ?');
        ScanDataAsNum_Raw = str2num(ScanDataAsText);

        HF_VDC(i) = -ScanDataAsNum_Raw(1);
        HF_eng(i) = -ScanDataAsNum_Raw(1).*HFG_conversion;

        mass_query = query(loadcell, 'SI');
        mass_query = split(mass_query);
        mass(i) = str2double(cell2mat(mass_query(3)));
        time_stamp(i) = str2double(datestr(now, TimeStampFormat));
        test_time(i) = (datenum(num2str(time_stamp(i), '%.3f'), TimeStampFormat)-datenum(num2str(time_stamp(1), '%.3f'), TimeStampFormat)).*(24*60*60);

        % calcualte the mass loss rate and the smooth mass loss rate
        MLR = [0 (-diff(mass)./diff(test_time))];
        MLR = MLR./sample_area;
        %         if MLR(i) < 0
        %             MLR(i) = 0;
        %         end

        time_change(i) = test_time(i)-test_time(i-1);
        smooth_MLR(i) = mean(MLR(i-averaging_window:i));

        VDC_out(i) = i*irradiation_rate; %VDC_init
        proportional_term(i) = 0;
        integral_term(i) = 0;
        derivative_term(i) = 0;


        if exist('PID_status', 'var')==0
            if mean(smooth_MLR(i-ave_frame:i)) > 0.8*desired_MLR %isempty(find(smooth_MLR) > 0.8*desired_MLR, 1))==0
                PID_status = 1;
            end
        end

        if exist('PID_status', 'var') == 1

            disp('PID active')

            current_input = smooth_MLR(i);
            last_input = smooth_MLR(i-1);

            %calculate the PID terms
            error(i) = setpoint - current_input;
            error_sum(i) = error(i) * time_change(i);
            d_input(i) = (current_input - last_input) / time_change(i);

            % compute all the terms
            proportional_term(i) = PID_kp*error(i);
            integral_term(i) = integral_term(i-1) + PID_ki*error_sum(i);
            derivative_term(i) = PID_kd * d_input(i);

            % constrain the integral term so that the PID understand the limits of the lamps
            if integral_term(i) > max_lamp_voltage
                integral_term(i) = max_lamp_voltage;
            end
            if integral_term(i) < min_lamp_voltage
                integral_term(i) = min_lamp_voltage;
            end

            % Compute the Output
            output(i) = proportional_term(i) + integral_term(i) + derivative_term(i);

            % protect the FPA!
            if output(i) > max_lamp_voltage
                output(i) = max_lamp_voltage;
            end
            if output(i) < min_lamp_voltage
                output(i) = min_lamp_voltage;
            end

            VDC_out(i) = output(i);

        end

        disp(['VDC to lamps = ' num2str(VDC_out(i)) ' V']);
        disp(['Mass loss rate = ' num2str(smooth_MLR(i)) ' g/s/m2']);
        disp(['The test time is ' num2str(test_time(i)) ' s.'])

        HF_plot = plot(test_time(i), VDC_out(i), 'x', 'parent', HFG_axes, 'color', 'r');
        mass_plot = plot(test_time(i), MLR(i), 'o', 'parent', mass_axes, 'color', [0.5 0.5 0.5]);
        mass_plot = plot(test_time(i), smooth_MLR(i), '.', 'parent', mass_axes, 'color', 'b');
            
        fprintf(DL,['SOUR:VOLT ' num2str(VDC_out(i)), ', (@' DL_CHANNELS_VDC_OUTPUT ')']);         % Set the output voltage level on the specified DAC channel. You can set the output voltage to any value between +12 Vdc and -12 Vdc, in 1 mV steps. Each DAC channel is capable of 10 mA maximum output current. The DAC channels are numbered “s04” and “s05”, where s represents the slot number. The :VOLT? query returns the output voltage level on the specified DAC channel. Returns a number in the form “+8.00000000E+00”.

        i=i+1;
        toc
        pause(1-toc)

    end

    i

end


fprintf(DL,['SOUR:VOLT 0, (@' DL_CHANNELS_VDC_OUTPUT ')']);
fclose(DL)
fclose(loadcell)

figure
hold on
plot(test_time, proportional_term)
plot(test_time, derivative_term)
plot(test_time, integral_term)

save([FileName '.mat'])
% mass_loss_rate = -diff(mass)./diff(test_time);
% smooth_MLR =  fit(test_time(1:end-1)',mass_loss_rate','smoothingspline', 'SmoothingParam',0.05);
% mlr_fig = figure;
% plot(smooth_MLR, test_time(1:end-1)', mass_loss_rate')
% xlabel('Test time, s')
% ylabel('Mass loss rate, g/s')
%
% SMOOTH_MLR = smooth(mass_loss_rate, 21, 'sgolay', 2);
% SMOOTH_MLR = [0
%     SMOOTH_MLR];
% mass_loss_rate = [0 mass_loss_rate];
% test_data = [test_time' mass' mass_loss_rate' Ramp_HF' SMOOTH_MLR];
% writematrix(test_data, [FileName '_test.csv'])
%
% save([FileName '.mat'])
% saveas(test_fig, [FileName '_test.png'])
% saveas(mlr_fig, [FileName '_mlr.png'])



