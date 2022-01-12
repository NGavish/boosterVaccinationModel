%% This script produces the graphs
close all;
defineColors;
pBoostColor=green;
presentCalibratedModel=false; % Figure 1: Calibrated model outcomes vs actual data.
presentBoosterSchedule=false; % Figure 2: Effect of waning and boosting.
presentEffectOfTiming=false;  % Figure 3: Timing of boost start date.
presentEffectOfBoost=false;   % Figure 6 (SM): Booster uptake in Israel.
presentVaccineWaning=true;    % Figure 7: Estimates of vaccine waning profile.
presentVEprofile=false;       % Figure 8: Vaccine efficacy profile from administration of first dose.

displaystartTime=datetime('1/7/2021','InputFormat','dd/MM/uuuu');
displayendTime=datetime('15/11/2021','InputFormat','dd/MM/uuuu');
close all;
%% presentCalibratedModel
if presentCalibratedModel

    startTimeCalibration=datetime('1/7/2021','InputFormat','dd/MM/uuuu');displaystartTimeCalibration=startTimeCalibration;
    endTimeCalibration=datetime('14/11/2021','InputFormat','dd/MM/uuuu');
    T = readtable('calibration.csv');
    TdataByVacStatus = readtable('calibrationData.csv');

    t_layout = tiledlayout(2,1,'TileSpacing','Compact');
    t_detected = tiledlayout(1,3,'TileSpacing','Compact');

    ax1=nexttile(t_detected);
    h=plot(T.dates,T.detected_nv,LineWidth=2,color=pBoostColor);hold on;
    hold on;
    p2=plot(TdataByVacStatus.dates,TdataByVacStatus.detected_nv,'k.','LineWidth',1.2);
    p3=plot(TdataByVacStatus.dates,movmean(TdataByVacStatus.detected_nv,7),'k','LineWidth',1)
    xlim([displaystartTimeCalibration displayendTime]);grid on;
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 8500]);
    ylabel('New detected cases (# of people)',FontWeight='bold')

    ax2=nexttile(t_detected);
    h=plot(T.dates,T.detected_v,LineWidth=2,color=green);hold on;
    hold on;
    plot(TdataByVacStatus.dates,TdataByVacStatus.detected_v,'k.','LineWidth',1.2);
    plot(TdataByVacStatus.dates,movmean(TdataByVacStatus.detected_v,7),'k','LineWidth',1)
    xlim([displaystartTimeCalibration displayendTime]);grid on;
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 7500]);

    ax3=nexttile(t_detected);
    h=plot(T.dates,T.detected_b,LineWidth=2,color=pBoostColor);hold on;
    hold on;
    p2=plot(TdataByVacStatus.dates,TdataByVacStatus.detected_b,'k.','LineWidth',1.2);
    p3=plot(TdataByVacStatus.dates,movmean(TdataByVacStatus.detected_b,7),'k','LineWidth',1)
    xlim([displaystartTimeCalibration displayendTime]);grid on;
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 7500]);

    sgtitle(t_detected,'Calibration - PCR confirmed detection data',fontWeight='bold')

    set(gcf,'Position',[520 143 873 579])

    lg=legend(ax3,'Model outcome','Actual data','Actual Data - moving average','Orientation','horizontal','fontsize',12)
    lg.Layout.Tile='South'
    text(ax1,0.05,0.9,'A','FontSize',14,'Units','normalized')
    text(ax2,0.05,0.9,'B','FontSize',14,'Units','normalized')
    text(ax3,0.05,0.9,'C','FontSize',14,'Units','normalized')

    title(ax1,'Non-vaccinated','FontWeight','normal','fontsize',12);%title(ax4,'Non-vaccinated','FontWeight','normal','fontsize',12)
    title(ax2,'Vaccinated','FontWeight','normal','fontsize',12);%title(ax5,'Vaccinated','FontWeight','normal','fontsize',12)
    title(ax3,'Vaccinated with booster','FontWeight','normal','fontsize',12);%title(ax6,'Vaccinated with booster','FontWeight','normal','fontsize',12)

    set(gcf,'Position',[520 402 873 320])
    print(['CalibratedModel',datestr(endTimeCalibration)],'-depsc');
    print(['CalibratedModel',datestr(endTimeCalibration)],'-djpeg');

end
close all;
%% presentBoosterSchedule
if presentBoosterSchedule
    colororder('default')

    subplot(3,3,[1 2 4 5])
    T = readtable('EffectOfBoost.csv');
    plot(T.dates,100*cumsum(T.BoostData60)/1287420,LineWidth=1.5);hold on;
    h=plot(T.dates,100*cumsum(T.BoostData4059)/1471104,'--',LineWidth=1.5);hold on;
    plot(T.dates,100*cumsum(T.BoostData039)/1977213,'-.',LineWidth=1.5)

    displaystartTime=datetime('30/7/2021','InputFormat','dd/MM/uuuu');
    displayendTime=datetime('1/12/2021','InputFormat','dd/MM/uuuu');
    xlim([displaystartTime displayendTime]);
    set(gca,'xtick',displaystartTime+[0,14,19,25,30,47,63,78,93,108],'xticklabelrotation',45);%,'xticklabel',{'1-Aug, Ages 60+ elibigle','13-Aug, Ages 50+ elibigle','20-Aug, Ages 40+ elibigle','24-Aug, Ages 30+ elibigle','29-Aug, Ages 16+ elibigle','15-Sep','1-Oct','15-Oct'})
    ytickformat('percentage')
    ylim([0,100])
    text(displaystartTime,100,{' Ages 60+',' eligible for Boost'},'VerticalAlignment','bottom','Rotation',-90)
    text(displaystartTime+14,100,' Ages 50+','VerticalAlignment','bottom','Rotation',-90)
    text(displaystartTime+19,100,' Ages 40+','VerticalAlignment','bottom','Rotation',-90)
    text(displaystartTime+25,100,' Ages 30+','VerticalAlignment','bottom','Rotation',-90)
    text(displaystartTime+30,100,' Ages 16+','VerticalAlignment','bottom','Rotation',-90)
    legend('Ages 60 and older','Age group 40-59','Age group 16-39','location','southeast','autoupdate','off')
    plot(displaystartTime+[0,0],[0,100],'k','LineWidth',1)
    for ix=[0,14,19,25,30]
        plot(displaystartTime+[0,0]+ix,[0,100],'--','LineWidth',0.7,'Color',0.7*[1 1 1])
    end
    grid on
    ylabel('Portion of vaccinated who received a booster shot')
    set(gcf,'Position',[520 318 768 479])

    print -depsc BoostSchedule
end
close all;
%% presentEffectOfBoost
if presentEffectOfBoost
    close all;

    displaystartTime=datetime('15/7/2021','InputFormat','dd/MM/uuuu');
    displayendTime=datetime('1/12/2021','InputFormat','dd/MM/uuuu');
    defineColors;
    newcolors =[orange; yellow; blue; green];

    colororder(newcolors)


    T = readtable('EffectOfBoost.csv');

    tiledlayout(3,1,'TileSpacing','tight')
    ax1=nexttile;

    hNoBoost=plot(T.dates,T.detected_NoBoost,'--',LineWidth=2);hold on;
    h60=plot(T.dates,T.detected_Boost60,'-.',LineWidth=2);
    h40=plot(T.dates,T.detected_Boost40,':',LineWidth=2);hold on;
    h16=plot(T.dates,T.detected_Boost,'-',LineWidth=2);hold on;

    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    MaxY=40000;ytickformat('%,4.0g');ylim([0 MaxY]);
    boosterStartTime=datetime('30/7/2021','InputFormat','dd/MM/uuuu');
    text(boosterStartTime,MaxY,{' Ages 60+',' eligible for Boost'},'VerticalAlignment','bottom','Rotation',-90)
    text(boosterStartTime+14,MaxY,' Ages 50+','VerticalAlignment','bottom','Rotation',-90)
    text(boosterStartTime+21,MaxY,' Ages 40+','VerticalAlignment','bottom','Rotation',-90)
    text(boosterStartTime+19,MaxY,' Green pass','VerticalAlignment','bottom','Rotation',-90)
    text(boosterStartTime+25,MaxY,' Ages 30+','VerticalAlignment','bottom','Rotation',-90)
    text(boosterStartTime+30,MaxY,' Ages 16+','VerticalAlignment','bottom','Rotation',-90)

    for ix=[0 14 19 21 25 30]
        h=plot(boosterStartTime+[0,0]+ix,[0,45000],'--','LineWidth',0.7,'Color',0.7*[1 1 1])
    end
    h.Color=0.4*[1 1 1]
    ylabel('Number of cases')
    title('Daily number of detected cases','fontsize',13)
    text(0.025,0.92,'A','FontSize',14,'Units','normalized')
    grid on;

    ax2=nexttile;
    plot(T.dates,100+0*T.severe_NoBoost,'Color',0.6*[1 1 1],LineWidth=1);hold on;
    text(T.dates(20),105,{'Lockdown threshold','in previous surges'},'VerticalAlignment','bottom','HorizontalAlignment','left')

    plot(T.dates,T.severe_NoBoost,'--',LineWidth=2,Color=hNoBoost.Color);hold on;
    h60=plot(T.dates,T.severe_Boost60,'-.',LineWidth=2,Color=h60.Color);hold on;
    h40=plot(T.dates,T.severe_Boost40,':',LineWidth=2,Color=h40.Color);hold on;
    h12=plot(T.dates,T.severe_Boost,LineWidth=2,Color=h16.Color);hold on;

    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 800]);
    ylabel('Number of cases');title('Daily number of new severe cases','fontsize',13)
    text(0.025,0.92,'B','FontSize',14,'Units','normalized')

    grid on;
    ax3=nexttile;

    plot(T.dates,T.R_NoBoost,'--',LineWidth=2);hold on;
    plot(T.dates,T.R_Boost60,'-.',LineWidth=2)
    plot(T.dates,T.R_Boost40,':',LineWidth=2);hold on;
    plot(T.dates,T.R_Boost,LineWidth=2);hold on;

    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0.65 1.4]);

    grid on
    ylabel('R_{eff}','Rotation',0)
    title('Effective reproduction number','fontsize',13)

    text(0.025,0.92,'C','FontSize',14,'Units','normalized')

    lg=legend('No boosters','Boosters - ages 60 and older','Boosters - ages 40 and older','Boosters - ages 16 and older','Orientation','horizontal');%,'Fontsize',12);
    lg.Layout.Tile='north'
    linkaxes([ax1,ax2,ax3],'x')
    set([ax1,ax2,ax3],'xtick',boosterStartTime+[0,14,19,21,25,30,47,63,77,94,108,124],'xticklabelrotation',45);

    set(gcf,'Position',[520 182 762 615])
    shg

    ax1=axes('Position',[0.60 0.47479674796748 0.278215223097113 0.12520325203252]);
    plot(T.dates,100+0*T.severe_NoBoost,'Color',0.6*[1 1 1],LineWidth=1);hold on;

    plot(T.dates,T.severe_NoBoost,'--',LineWidth=2,Color=hNoBoost.Color);hold on;
    plot(T.dates,T.severe_Boost60,'-.',LineWidth=2,Color=h60.Color);hold on;
    plot(T.dates,T.severe_Boost40,':',LineWidth=2,Color=h40.Color);hold on;
    plot(T.dates,T.severe_Boost,LineWidth=2,Color=h12.Color);hold on;

    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    ytickformat('%,4.0g');ylim([0 150]);
    set(gca,'xtick',boosterStartTime+[0,1,17,19,23,28,31,45,61,75,92,106,122],'XTickLabel', []);

    grid on
    shg

    printGraph('EffectOfBoost');

    close
    colororder('default')
end
close all;
if presentEffectOfTiming
    displaystartTime=datetime('15/7/2021','InputFormat','dd/MM/uuuu');
    displayendTime=datetime('15/11/2021','InputFormat','dd/MM/uuuu');    %% Effect of Timing
    T = readtable('EffectOfTiming.csv');
    tiledlayout(3,1,'TileSpacing','tight')
    ax1=nexttile;
    plot(T.dates,T.detected,LineWidth=2);hold on;
    plot(T.dates,T.detected_Early,'-.',LineWidth=2)
    plot(T.dates,T.detected_Late,'--',LineWidth=2)
    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 21000]);
    ylabel('Number of cases')
    title('Daily number of new cases','fontsize',13)
    grid on
    text(0.025,0.92,'A','FontSize',14,'Units','normalized')
    ax2=nexttile;
    p1=plot(T.dates,T.severe,LineWidth=2);hold on;
    p2=plot(T.dates,T.severe_Early,'-.',LineWidth=2)
    p3=plot(T.dates,T.severe_Late,'--',LineWidth=2)
    xlim([displaystartTime displayendTime]);
    ax = gca;ax.YRuler.Exponent = 0;
    xtickformat('dd-MM')
    ytickformat('%,4.0g');ylim([0 215]);

    ylabel('Number of cases')
    title('Daily number of new severe cases','fontsize',13)
    text(0.025,0.92,'B','FontSize',14,'Units','normalized')
    grid on
    ax3=nexttile;
    plot(T.dates,T.Reff,LineWidth=2);hold on;
    plot(T.dates,T.R_early,'-.',LineWidth=2)
    plot(T.dates,T.R_late,'--',LineWidth=2)

    xlim([displaystartTime displayendTime]);
    ylim([0.55 1.4])
    ylabel('R_{eff}','Rotation',0)
    title('Effective reproduction number','fontsize',13)
    xtickformat('dd-MM')
    text(0.8,0.9,'C','FontSize',14,'Units','normalized')
    grid on
    set(gcf,'Position',[520 182 675 615])
    lg=legend('Actual boost schedule','Early schedule','Late schedule','Orientation','horizontal','Fontsize',12);
    lg.Layout.Tile='north'
    linkaxes([ax1,ax2,ax3],'x')
    set([ax1,ax2,ax3],'xtick',displaystartTime+[0 17 31 48 62 78 92 109],'xticklabelrotation',45);

    printGraph('EffectOfTiming')
    
end
close all;
if presentVaccineWaning
    set(0, 'DefaultLineLineWidth', 1.5);
    subplot(2,1,1)
    tt=35+[240 210 180 150 120 90  60  30  0]+15;

    yairNewJul=[nan nan 1.2 1.5 2 2.8 3.1 4.5 8.2];
    yairNewAug=[nan 1.2 1.4 1.7 2.2 2.5 3.8 6.5 6.8];
    yairNewSep=[1.5 1.5 1.7 2.2 2.5 3.8 5.4 6 8.7];

    %     eff=linspace(93,53,5);

    plot(35+[0 30 30 60 60 90 90 120 120 150],[93 93 87 87 77 77 59 59 53 53],'k-','MarkerSize',10);hold on;

    plot(tt,100*(1-1./yairNewJul),'+');
    plot(tt,100*(1-1./yairNewAug),'+');
    plot(tt,100*(1-1./yairNewSep),'+');

    plot(35+(0:300),linspace(97,0,301),'-.')

    legend('Tartof et al.','Israel - July','Israel - August','Israel - September','Best fitted VE profile')
    ytickformat('percentage');ylabel('Vaccine efficacy');xlabel('s');
    title('Vaccine efficacy profile in preventing infections')
    text(0.05,0.9,'A','units','normalized','FontSize',12)
    set(gca,'xtick',[0 35 100:50:300])

    subplot(2,1,2)
    ts=35+[120 150 180]
    yair4059severe=[98 98 94];
    yair60severe=[91 88 86];

    h=plot(ts,yair4059severe,'s','linewidth',1.5,color=blue);hold on;
    t=(0:300);
    hold on;plot(35+t,100-0.025*t,'--',color=gray)

    h=plot(ts,yair60severe,'d','linewidth',1.5,color=red);hold on;
    hold on;plot(35+t,100-0.078*t,color='k')

    legend('40-59','Best fitted curve to 40-59','60+','Best fitted curve to 60+','autoupdate','off','location','best')
    ylim([80 100]);ylabel('Vaccine efficacy');xlabel('s')
    xlim([0 210+35])
    set(gca,'xtick',[0 35 100:50:300]);ytickformat('percentage');ylim([75 100])
    shg
    text(0.05,0.9,'B','units','normalized','FontSize',12)
    title('Vaccine efficacy profile in preventing severe outcomes')
    printGraph('VE');
end
close all;
if presentVEprofile
    load('sigmaVaccine.mat');
    t=1:365;
    subplot(2,2,1);

    plot(t,100*VEv)
    xlabel('s');ylabel('VE_v(s)       ','rotation',0);xlim([0 365]);ylim([0 100])
    title({'Vaccine efficacy profile','after first dose vaccination'})
    grid on
    set(gca,'xtick',[0 35 100:50:300])
    text(0.9,0.9,'A','units','normalized','FontSize',12)
    ytickformat('percentage');
    subplot(2,2,2);

    plot(t,100*VEb)
    xlabel('s');ylabel('VE_b(s)       ','rotation',0);xlim([0 365]);ylim([0 100])
    title({'Vaccine efficacy profile','after booster vaccination'})
    ytickformat('percentage');
    grid on
    text(0.9,0.9,'B','units','normalized','FontSize',12)
    set(gca,'xtick',[0 35 100:50:300])
    xlim([0 120])

    printGraph('RelativeSusceptibilityProfile')
end
return
