%% Read NEX5
addpath(genpath('D:\Desktop\fengrui_solve'));
clear
clc
[file,path]=uigetfile('*.nex5');
cd(path)
UnitsData= readNex5File(fullfile(path,file));
TimeVector = 0:0.001:UnitsData.tend;
OptoMarks = UnitsData.events{1,1}.timestamps;
OptoMarks([1:100])=[];
EventMarks =[];
clear path
%% Opto-tagging
addpath(genpath('D:\Desktop\CellBase-master\CellBase_R2013a'));
InputParameters = inputdlg({'Pre (ms)';'Post (ms)';'Testing Window (ms)'},...
    'Peri-opto Range',[1 60],{'20';'20';'5'});
PreEvent = abs(str2num(InputParameters{1})/1000);
PostEvent =str2num(InputParameters{2})/1000;
WN = str2num(InputParameters{3})/1000;
SpikeTrain=zeros(numel(UnitsData.neurons),length(TimeVector));
Peri_Light = zeros(length(OptoMarks),(PreEvent+PostEvent)*1000+1,numel(UnitsData.neurons));
clear OptoTagged
close all
for u = 1:numel(UnitsData.neurons)
    idx = round(UnitsData.neurons{u,1}.timestamps/0.001)+1;
    SpikeTrain(u,idx)=1;
    for i=1:length(OptoMarks)
        event_idx = round(OptoMarks(i)/0.001) + 1;
        idx_range = (event_idx - PreEvent*1000):(event_idx + PostEvent*1000);
        Peri_Light(i,:,u) = SpikeTrain(u,idx_range);
        clear idx_range event_idx
    end
    TempColor = [0.5 0.5 0.5];
    BL = Peri_Light(:,1:PreEvent*1000,u);
    TS = Peri_Light(:,end-(PostEvent*1000)+1:end,u);
    [P I]=salt(BL,TS,0.001);
    OptoTagged{u,1}=UnitsData.neurons{u,1}.name;
    OptoTagged{u,2}=0;
    if P<0.05 && mean(sum(TS(:,1:WN*1000),2))>mean(sum(BL(:,end-WN*1000+1:end),2)) && I>0.05
        TempColor = [0.9 0.1 0.2];
        OptoTagged{u,2}=1;
    end
    figure
    sgtitle(['NO. ',num2str(u),', ',UnitsData.neurons{u,1}.name,', I = ',num2str(I)],'color',TempColor)
    subplot(2,1,1)
    imagesc(-PreEvent*1000:PostEvent*1000,1:length(OptoMarks),Peri_Light(:,:,u))
    colormap([1 1 1;0 0 0])
    clim([0 1])
    subplot(2,1,2)
    stairs(-PreEvent*1000:PostEvent*1000,mean(Peri_Light(:,:,u),1))
    xlabel('Time (ms)')
    clear idx i TempColor P I BL TS
    pause(0.1)
end
clear u InputParameters WN PreEvent PostEvent
%% Peri-behavior
addpath('D:\Desktop\fengrui_solve');
InputParameters = inputdlg({'Pre (s)';'Post (s)';'Bin (ms)';'Baseline from (s)';'Baseline to (s)'},...
    'Peri-event Range',[1 60],{'-2';'3';'10';'-1';'0'});
PreEvent = abs(str2num(InputParameters{1}));
PostEvent =str2num(InputParameters{2});
Bin = str2num(InputParameters{3});
BL(1) = -abs(str2num(InputParameters{4}));
BL(2) = str2num(InputParameters{5});
Peri_Event_Time = -PreEvent:1/1000:PostEvent;
%SpikeTrain = zeros(numel(UnitsData.neurons),length(TimeVector)); 
Peri_Event = zeros(length(EventMarks),(PreEvent+PostEvent)*1000+1,numel(UnitsData.neurons));
clear FiringRate FiringRate_Z
close all
for u = 1:numel(UnitsData.neurons)
    idx = round(UnitsData.neurons{u,1}.timestamps/0.001)+1;
    SpikeTrain(u,idx)=1;
    for i=1:length(EventMarks)
        event_idx = round(EventMarks(i)/0.001) + 1;
        idx_range = (event_idx - PreEvent*1000):(event_idx + PostEvent*1000);
        Peri_Event(i,:,u) = SpikeTrain(u,idx_range);
        clear idx_range event_idx
    end

    [tempFR,~,TimePoints] = binSpikePerTrial(Peri_Event(:,:,u), Bin, -PreEvent*1000,1,2);
    FiringRate(u,:)=mean(tempFR,1)./10;
    FiringRate_Z(u,:) = (FiringRate(u,:)-mean(FiringRate(u,TimePoints>=BL(1)*1000 & TimePoints<BL(2)*1000)))./std(FiringRate(u,TimePoints>=BL(1)*1000 & TimePoints<BL(2)*1000));
    TempColor = [0.1 0.1 0.1];
    figure
    sgtitle(['NO. ',num2str(u),', ',UnitsData.neurons{u,1}.name],'color',TempColor)
    subplot(2,1,1)
    imagesc(-PreEvent*1000:PostEvent*1000,1:length(EventMarks),Peri_Event(:,:,u))
    colormap([1 1 1;0 0 0])
    clim([0 1])
    subplot(2,1,2)
    plot((-PreEvent*1000+Bin):Bin:PostEvent*1000,FiringRate(u,:))
    xlabel('Time (ms)')
    ylabel('Firing rate (Hz)')
    clear idx i TempColor tempFR n
    pause(0.1)
end
clear u InputParameters SpikeTrain
%% Excited or Inhibited or no response (1 -1 0)
addpath('D:\桌面\多通道\fengrui_solve');
InputParameters = inputdlg({'Baseline from (s)';'Baseline to (s)';'Event from (s)';'Event to'},...
    'ROI select',[1 60],{num2str(BL(1));num2str(BL(2));' ';' '});
BL(1) = str2num(InputParameters{1});
BL(2) = str2num(InputParameters{2});
EventTime(1)=str2num(InputParameters{3});
EventTime(2)=str2num(InputParameters{4});
E_or_I = zeros(numel(UnitsData.neurons),1);
clc
for u = 1:numel(UnitsData.neurons)
    BL_firing = Peri_Event(:,Peri_Event_Time>=BL(1) & Peri_Event_Time<BL(2),u);
    TS_firing = Peri_Event(:,Peri_Event_Time>=EventTime(1) & Peri_Event_Time<EventTime(2),u);
    if length(BL_firing)==length(TS_firing)
        [P,RR] = poisson_rate_test(sum(BL_firing,1),sum(TS_firing,1),length(BL_firing)/1000);
        All_P(u,1)=P;
        All_RR(u,1)=RR;
        disp(['Unit #',num2str(u)])
        if P<0.05
            if RR>1
                disp('Increase')
                E_or_I(u,1)=1;
            end
            if RR<1
                disp('Decrease')
                E_or_I(u,1)=-1;
            end
        end
        clear BL_firing TS_firing P RR
    else
        clc
        disp('Length is different !')
    end
end
clear u
%% Plor ROI and heatmap
IsOT=[];
for i = 1:size(OptoTagged,1)
   if OptoTagged{i,2}==1
   IsOT=[IsOT,i];
   end
end
Opto_Tagged_FR_Z = FiringRate_Z(IsOT,:);
Opto_Tagged_FR_O = FiringRate(IsOT,:);
Opto_Tagged_RR = All_RR(IsOT,1);
figure;
imagesc(TimePoints/1000,1:size(IsOT),FiringRate_Z(IsOT,:));
colorbar
colormap jet
clim([-5 5])
xticks([-5:1:5]) 