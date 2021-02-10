%  SPDX-License-Identifier: BSD-3-Clause
%
%  Copyright(c) 2020 Intel Corporation. All rights reserved.
%
%  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
%---------------------------------------------------
%---------------------------------------
%   History
%---------------------------------------
%   2020/12/24 Sriram Shastry       - initial version
%
clearvars;close all;clc;
WorkingDir = uigetdir();
Files      = subdir(WorkingDir);
clf();
figNo      = 100;

% figure(figNo+1);
for k = 1:length(Files)
    [~,foldername]   = fileparts(Files(k).folder);
    [~,filename,ext] = fileparts(Files(k).name);
    disp([foldername,    filename,    ext]);
    pattern = ["Results",];
    extn    = ["",'.txt'];
    % prepare for data plots
    
    
    if contains(foldername,pattern) && contains(filename,'drc_asin_fixed') && contains(ext,extn)  && Files(k).bytes > 1
        %  C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_asin_fixed\drc_asin_fixpt\Results
        %  * Input is Q2.30: (-2.0, 2.0)
        %  * Output range: [-1.0, 1.0]; regulated to Q2.30: (-2.0, 2.0)
        
        figure(figNo+1)
        drcasinfixed = get_drc_asin( Files(k).name);
        
        Index     = drcasinfixed.idx;
        Inptasine = drcasinfixed.inasine/2^30;
        Outasine  = drcasinfixed.outasine/2^30;
        
        x = [-1:0.1:1];
        y= 2/pi*asin(x);
        rms1 = rms(abs(y)) - rms(abs(Outasine));
        
        subplot(3,2,1);
        plot(1:numel(x), x,'DisplayName','input-fltpt-range');grid on;
        xlabel('Numpts');ylabel('fix(-2\pi):0.1:fix( 2\pi)');legend({'y = asin(x)'},'Location','best')
        title('Range-asin[Q2.30]');
        
        subplot(3,2,2);
        plot(1:numel(Index), Inptasine,'DisplayName','fixpt-range');grid on;
        xlabel('Numpts');ylabel('fix(-2\pi):0.1:fix( 2\pi)');legend({'y = asin(x)'},'Location','best')
        title('Range-asin[Q2.30]');
        
        
        subplot(3,2,3);
        plot(1:numel(x), y,'DisplayName','fltpt-Output');grid on;
        xlabel('Numpts');ylabel('fix(-2\pi):0.1:fix( 2\pi)');legend({'y = asin(x)'},'Location','best')
        title('Asin(floatpoint)');
        
        subplot(3,2,4);
        plot(1:numel(Index), Outasine,'DisplayName','fixpt-output');grid on;
        xlabel('Numpts');ylabel('fix(-2\pi):0.1:fix( 2\pi)');legend({'y = asin(x)'},'Location','best')
        title('Asin(fixpoint)');
        
        % Plot output 4
        subplot(3,2,5);
        polarplot(((-pi*1/pi:0.1:pi*1/pi)),y,'r--'); hold on;
        polarplot(((-pi*1/pi:0.1:pi*1/pi)),Outasine,'g-x'); hold on;legend('fltpt','fixpt')
        %         polarplot(deg2rad((-pi*1/pi:0.1:pi*1/pi)),y,'m-+'); hold on;
        %         polarplot(deg2rad((-pi*1/pi:0.1:pi*1/pi)),Outasine,'k--'); hold on;legend('fltpt','fixpt')
        
        pax = gca;
        pax.ThetaAxisUnits = 'Degree';
        pax = gca;
        pax.ThetaColor = 'magenta';
        pax.RColor = [0 .5 0];
        pax.GridColor = 'black';
        title('Degree-2/pi * asin(x)-floatpt');
        
        % Plot output 5
        subplot(3,2,6);
        
        plot(1:numel(x),20*log10(sqrt(mean(rms1.^2))),'k-+'); hold on; grid on; legend('THD+N');title('drc-sin-fixpt-SNR');
        xlabel('fix(-3\pi) < x <fix( 3\pi)');ylabel('SNR[dB]');
    end
    
    if contains(foldername,pattern) && contains(filename,'drc_sin_fixed') && contains(ext,extn)  && Files(k).bytes > 1
        % C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_sin_fixed\drc_sin_fixed\Results
        figure(figNo+2)
        drcsinfixed = get_drc_sine_fixed( Files(k).name);
        idx         = drcsinfixed.idx;
        fixptvector = drcsinfixed.insine/2.^30;    % in Radian
        fixptsine   = drcsinfixed.outsine/2.^31;
        
        %         x = fi((-1:0.1:1),1,32,30);
        x = (-1:0.1:1);
        y= sin(x*pi/2);
        rms1 = rms((y)) - rms((fixptsine));
        % Plot output 1
        subplot(3,2,1);
        plot(idx, x,'DisplayName','input-range');grid on;
        xlabel('fix(-3\pi) < x <fix( 3\pi)');ylabel('InputRange[AbsNum]');legend({'y = sin(x)'},'Location','best')
        title('Range-Sin[Q2.30]');
        
        % Plot output 2
        subplot(3,2,2);
        plot(x, sin(x*pi/2),'DisplayName','input-range');grid on;
        xlabel('fix(-3\pi) < x <fix( 3\pi)');ylabel('y = sin(x)');legend({'y = sin(x)'},'Location','best')
        title('Reference-Sin(x)-Floatpt');
        
        % Plot output 3
        subplot(3,2,3);
        plot(idx, fixptsine,'r-+','DisplayName','Johny-sine');grid on;
        xlabel('fix(-3\pi) < x <fix( 3\pi)');ylabel('y = sin(x)');legend({'y = sin(x)'},'Location','best')
        title( 'Johny-Fixpt-drc-Sin-fixed[Q1.31]');
        % Plot output 4
        subplot(3,2,4);
        polarplot(fixptvector,max(0,fixptsine),'+-r')
        pax= gca;
        pax.FontSize = 12;
        title('polar-sin(x*pi/2)-Sin-fixed[Q1.31]');
        % Plot output 5
        subplot(3,2,5:6);
        plot(idx,20*log10(sqrt(mean(rms1))),'k-+'); hold on; grid on; legend('THD+N');title('drc-sin-fixpt-SNR');
        xlabel('fix(-3\pi) < x <fix( 3\pi)');ylabel('SNR[dB]');
    end
    
    if contains(foldername,pattern) && contains(filename,'drc_pow_fixed') && contains(ext,extn)  && Files(k).bytes > 1
        % C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_pow_fixed\drc_power_fixed\Results
        figure(figNo+3)
        drcpowfixed = get_drc_pow( Files(k).name);
        
        Idxi   = drcpowfixed.idxi;  % Base
        Idxj   = drcpowfixed.idxj;  % Exponent
        InputX = drcpowfixed.inxQ626/2^26;
        InputY = drcpowfixed.inyQ230/2^30;
        Output = drcpowfixed.outpowQ1220/2^20;
        
        subplot(3,2,1);plot(1:numel(Idxi),Idxi,'r--','linewidth',1.5); hold on;  grid on
        xlabel('Base[Numpts] ');ylabel('Mag(AbsVal)');legend('Base- exponent', 'location', 'best');
        title('pow(Base,Exp)=Output:Base-Q6.26 [>0To32.0], Exponent-Q2.30 [-2.0, 2.0]',...
            'Output[Q12.20-max 2048.0]-WithOutScaling')
        
        
        subplot(3,2,2);plot(1:numel(Idxj),Idxj,'g--','linewidth',1.5); hold on;  grid on
        xlabel('Exponent[Numpts] ');ylabel('Mag(AbsVal)');legend('Base- exponent', 'location', 'best');
        title('pow(Base,Exp)=Output:Base-Q6.26 [>0To32.0], Exponent-Q2.30 [-2.0, 2.0]',...
            'Output[Q12.20-max 2048.0]-WithOutScaling')
        
        subplot(3,2,3);plot(1:numel(InputX),InputX,'b--','linewidth',1.5); hold on;  grid on
        xlabel('Base[Numpts] ');ylabel('Mag[Q626]');legend('Base- exponent', 'location', 'best');
        title('pow(Base,Exp)=Output:Base-Q6.26 [>0To32.0], Exponent-Q2.30 [-2.0, 2.0]',...
            'Output[Q12.20-max 2048.0]-WithScaling')
        
        
        subplot(3,2,4);plot(1:numel(InputY),InputY,'c--','linewidth',1.5); hold on;  grid on
        xlabel('Exponent[Numpts] ');ylabel('Mag(Q230]');legend('Base- exponent', 'location', 'best');
        title('pow(Base,Exp)=Output:Base-Q6.26 [>0To32.0], Exponent-Q2.30 [-2.0, 2.0]',...
            'Output[Q12.20-max 2048.0]-WithScaling')
        
        
        subplot(3,2,5:6);plot(1:numel(Output),Output,'m--','linewidth',2.5); hold on;  grid on
        xlabel('Output[Numpts] ');ylabel('Power[Q12.20]');legend('Base- Output', 'location', 'best');
        title('pow(Base,Exp)=Output:Base-Q6.26 [>0To32.0], Exponent-Q2.30 [-2.0, 2.0]',...
            'Output[Q12.20-max 2048.0]-WithScaling')
        
    end
    if contains(foldername,pattern) && contains(filename,'drc_inv_fixed') && contains(ext,extn)  && Files(k).bytes > 1
        % C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_inv_fixed\drc_inv_fixed\Results
        figure(figNo+4)
        drcinvfixed = get_drc_inv( Files(k).name);
        
        Index = drcinvfixed.idx;
        InputInverse = drcinvfixed.ininv/2^12;
        OuputInverse = drcinvfixed.outinv/2^30;
        
        x = [[2^-16, 2^16-1], [-2^-16, -(2^16-1)]];
        y = pinv(x);
        
        subplot(2,2,1);
        x = 1:length(x);
        plot(x,(x),'b-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-input');title('fixpt-inverse');
        
        subplot(2,2,2);
        plot(x,(y),'r-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-Output');title('floating-inverse');
        
        subplot(2,2,3);
        Index = 1:length(Index);
        plot(Index,(InputInverse),'m-o','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-input');title('fixpt-inverse');
        
        subplot(2,2,4);
        plot(Index,(OuputInverse),'k-o','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-Output');title('fixpt-inverse');
        
    end
    
    if contains(foldername,pattern) && contains(filename,'drc_db2mag') && contains(ext,extn)  && Files(k).bytes > 1
        % C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_db2mag_fixed\drc_db2mag\Results
        %         ATTENTION - this covered
        figure(figNo+4)
        drcdb2mag = get_drc_db2mag( Files(k).name);
        
        drcdb2mag.Index = dataArray{:, 1};
        drcdb2mag.Invaldb = dataArray{:, 2};
        drcdb2mag.Outvalmag = dataArray{:, 3};
        
        x = [[2^-16, 2^16-1], [-2^-16, -(2^16-1)]];
        y = pinv(x);
        
        subplot(2,2,1);
        x = 1:length(x);
        plot(x,(x),'b-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-input');title('fixpt-inverse');
        
        subplot(2,2,2);
        plot(x,(y),'r-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-Output');title('floating-inverse');
        
        subplot(2,2,3);
        Index = 1:length(Index);
        plot(Index,(InputInverse),'m-o','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-input');title('fixpt-inverse');
        
        subplot(2,2,4);
        plot(Index,(OuputInverse),'k-o','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('Inv-Output');title('fixpt-inverse');
    end
    if contains(foldername,pattern) && contains(filename,'mag2dB') && contains(ext,extn)  && Files(k).bytes > 1
        %  * Input is Q6.26: max 32.0
        %  * Output range ~ (-inf, 30.1030); regulated to Q11.21: (-1024.0, 1024.0)
        
        % C:\Users\shastry\source\Work\Audio\SourceCode\a_v_03\Models\drc_Existing_CodeWrapper\drc_lin2db_fixed\drc_lin2dB_fixed\drc_lin2dB_fixed\Results
        figure(figNo+5)
        mag2dB = get_drc_mag2db( Files(k).name);
        
        Index     = mag2dB.idx;
        Inpmag2db = mag2dB.testvector/2^26;
        Outmag2db = mag2dB.Fixlog10linear/2^21;
        
        x = (-31.9:0.1:31.9);
        x(x<0) = NaN;
        y = mag2db(x);
        
        subplot(2,2,1);
        plot(x,y,'b-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('Mag[dB]');
        legend('mag2dB');title('fltpt-mag2db');
        
        subplot(2,2,2);
        plot(Inpmag2db,Outmag2db,'r-x','linewidth',2.0); hold on; grid on;
        xlabel('Numpts');ylabel('InvMag[x]');
        legend('mag2db');title('fixpt-mag2db');
        
        subplot(2,2,3:4);
        Error = rms(y(321:end)) - rms(Outmag2db(321:end));
        plot(Index,20*log10(sqrt(mean(Error.^2))),'r-x'); hold on; grid on;
        xlabel('Numpts');ylabel('THD+N');legend('THD+N')
        title('fltpt-fixpt-mag2db');
    end
    
end
clearvars
