% Gengxi Lu
% NRD: raw data
% NEV: event data
% Sampling frequency 30 kHz
% -ADMaxValue 8388608
% -ADBitVolts 0.000000015624999999999999 || 0.015625 uV
% testing data source: https://figshare.com/articles/dataset/UltrasoundRP/25438252

parpool(4)

clear

FilePath=['H:\neuralynx processing tool\CheetahData\'];
cd(FilePath)
FileNum=['_0007'];  % _0001

%%
fs = 30000;
dV = 0.015625; % uV
low_B = fir1(256, [1 300]./(fs/2),'bandpass'); % low frequency range for LFP
high_B = fir1(256, [500 9000]./(fs/2),'bandpass'); % high frequency range for spikes
Notch60 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);   % 60Hz notch for LFP
 
%% read event data
EvFilename= [FilePath,'Events', FileNum, '.nev'];   % Events || Events_00XX  .nev
EvFieldSelection = [1 0 0 0 0];  
% FieldSelectionFlags(1): Timestamps
% FieldSelectionFlags(2): Event IDs
% FieldSelectionFlags(3): TTLs
% FieldSelectionFlags(4): Extras
% FieldSelectionFlags(5): Event Strings
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];

% [TimeStamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV( EvFilename, EvFieldSelection, ExtractHeader, ExtractMode, ModeArray );
[EvTimestamps] = Nlx2MatEV( EvFilename, EvFieldSelection, ExtractHeader, ExtractMode, ModeArray );

Ev_T = (EvTimestamps-EvTimestamps(1))./1e6;   % unit: s. TimeStamp's unit is 'us'
Ev_P = round(Ev_T.*fs);   % point number
% Ev_P=[0 77*30 77*30 2500*30 2500*30 0];
t_range(:,1) = [Ev_P(4)-0.1*fs:Ev_P(4)+1*fs];
for i=2:length(Ev_P)/2-2
    t_range(:,i) = [Ev_P((i+1)*2)-0.1*fs:Ev_P((i+1)*2)+1*fs];
end


%% read NRD data

% initialize processed data
SpikeLine =zeros([size(t_range),64]);

FileName=[FilePath,'RawData', FileNum,'.nrd'];      % RawData || RawData_00XX   .nrd

FieldSelectionFlags = [0 1];  % FieldSelectionFlags(1): Timestamps; FieldSelectionFlags(2): Samples
HeaderExtractionFlag = 0;
ExtractMode = 1;
ExtractionModeVector = [];
% [Samples] = Nlx2MatNRD( FileName, 0,FieldSelectionFlags, HeaderExtractionFlag, ExtractMode, ExtractionModeVector);  
% figure, plot(filtfilt(high_B, 1, Samples));
% [Timestamps, Samples, Header] = Nlx2MatNRD( FileName, ChannelNumber,FieldSelectionFlags, HeaderExtractionFlag, ExtractMode, ExtractionModeVector);
for ChannelNumber = 0:31
    [Samples] = Nlx2MatNRD( FileName, ChannelNumber,FieldSelectionFlags, HeaderExtractionFlag, ExtractMode, ExtractionModeVector);  
    % spike and LFP
    S_high = filtfilt(high_B, 1, Samples);
    for i=1:size(t_range,2)
        % spike
        SpikeLine(:,i,ChannelNumber+1) = dV.*S_high(t_range(:,i));
        
    end
end 

for ChannelNumber = [0:23]+32
    [Samples] = Nlx2MatNRD( FileName, ChannelNumber+32,FieldSelectionFlags, HeaderExtractionFlag, ExtractMode, ExtractionModeVector);  
    % spike and LFP
    S_high = filtfilt(high_B, 1, Samples);
    for i=1:size(t_range,2)
        % spike
        SpikeLine(:,i,ChannelNumber+1) = dV.*S_high(t_range(:,i));
        
    end
end 

% save([FilePath, '_Processed data', FileNum,'.mat']);

%% plot figures

% load([FilePath, '_Processed data', FileNum,'.mat']);

% LFP
% figure, 
% hold on
% for i=1:32
%     plot([-0.1:1/fs:1], LFPLine(:,2,i)+(i-1)*500, 'LineWidth', 2); 
% %     plot(LFPLine(:,2,i)+(i-1)*300, 'LineWidth', 2); 
% end
% ylim([-500 500*33])
% xlabel('Time (s)');
% ylabel('Volt (uV)');
% set(gca,'LooseInset',get(gca,'TightInset'));
% set(gca,'box','on','FontSize',14,'fontweight','bold','Linewidth',2);
% set (gcf,'Position',[100,100,500,900]);
% saveas(gcf,[FilePath, FileNum, 'LFP.fig']);


% spike
% remove initial artifact in spikes
SpikeLine([0:0.012*fs]+0.1*fs,:,:) = SpikeLine([-0.012*fs:0]+0.1*fs,:,:);
SpikeLine([0.9248*fs:0.9319*fs]+0.1*fs,:,:) = SpikeLine([0.932*fs:0.9391*fs]+0.1*fs,:,:);


% plot raw signals
RNo = 1;

figure,  
hold on
for i=1:56
    plot([-0.1:1/fs:1], SpikeLine(:,RNo,i)./80+i,'LineWidth', 0.8); 
end
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'box','on','FontSize',14,'fontweight','bold','Linewidth',2);
set (gcf,'Position',[100,100,700,800]);
xlim([-0.1 1])
ylim([0 57])
% title([FileNum,'SpikeRawData NormNo.4 80uV Repeat',num2str(RNo)]);
% xlabel('Time (s)');
% ylabel('Channels');
% saveas(gcf,[FilePath, FileNum, 'Repeat',num2str(RNo),'Spike.fig']);

%% plot spike raster or MUA (multi-unit activity) of all channels
NoiseLvl = squeeze(sqrt(mean((mean(SpikeLine([1:0.098*fs],:,:).^2,2)),1)));  % calculate threshold for each channel

spike = zeros(size(SpikeLine));
spike_p = spike;
for i=1:56
    spike(:,:,i) = abs(SpikeLine(:,:,i))./NoiseLvl(i);
end
threshold = 3;
spike(spike<threshold)=threshold;  
for m = 1:56
    for n = 1:size(t_range,2)
        [temp,loc] = findpeaks(spike(:,n,m));
        spike_p(loc,n,m) = 1;
        clear temp
    end
end

figure,
hold on
for i=1:32
    plot([-0.1:1/fs:1],spike_p(:,RNo,i)*i,'sk','MarkerFaceColor','k','MarkerSize',1.8); %'MarkerEdgeColor',[0.5,0.5,0.5]
end
hold off
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'box','on','FontSize',14,'fontweight','bold','Linewidth',2);
set (gcf,'Position',[100,100,700,800]);
xlim([-0.1 1]);
ylim([0.5 32.5]);
xlabel('Time (s)');
ylabel('Channels');
% title('Spike Raster Representitive Normal');
title(['SpikeRaster NormNo.4 Threshold',num2str(threshold),'Repeat',num2str(RNo)]);
% saveas(gcf,[FilePath, FileNum, 'SpikeRaster', 'R',num2str(RNo), '.fig']);

%% plot spike count every 40 ms; need to run after spike raster
bins = 5*30; %40ms, 30 points per ms
counts = floor(size(t_range,1)/bins);
spike_c = zeros(counts,size(spike_p,2), size(spike_p,3));
for i=1:counts
    spike_c(i,:,:)=sum(spike_p([1:bins]+(i-1)*bins,:,:),1);
end
spike_count=squeeze(mean(spike_c(:,[3,4,6:8],:),2)); % for FileNum18 no noise repeats: [1:6,9:11,14,15,17,20]
std_count=squeeze(std(spike_c(:,[3,4,6:8],:),0,2));    

% response trials [4,5,7:11,14,15, 37,64]
% all channels
figure,
imagesc([-0.1 1] , [1 56], spike_count(:,1:56)',[0 10]);
colormap('hot')
% xlabel('Time (s)');
% ylabel('Channels');
set(gca,'YDir','normal')
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'box','on','FontSize',14,'fontweight','bold','Linewidth',2);
set (gcf,'Position',[100,100,700,800]);
saveas(gcf,[FilePath, FileNum, 'Avg SpikeCount All Channels_5ms.fig']);
saveas(gcf,[FilePath, FileNum, 'Avg SpikeCount All Channels_5ms.tif']);

% one channel
ChNo = 36;
figure,
hold on
bar(linspace(-0.1,1,counts),spike_count(:,ChNo)+0.7*std_count(:,ChNo),1,'Facecolor','c');
bar(linspace(-0.1,1,counts),spike_count(:,ChNo),1,'Facecolor','b','Facecolor','b');
ylim([0 20])
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'box','on','FontSize',14,'fontweight','bold','Linewidth',2);
set (gcf,'Position',[200,200,800,500]);
saveas(gcf,[FilePath, FileNum, 'SpikeCount_5ms Ch.',num2str(ChNo),'.fig']);


%% mapping
Interval = 0.35; % mm
% TrialNum=[1,5,6,7,8,9,10,11];
TrialNum=[1:11];
% ResponseAmp = sum(spike_count(21:60,1:56),1);
for i = [1,6:11]
%     ResponseAmp = mean(squeeze(sum(spike_c(21:60,1:TrialNum(i),1:56),1)),1);
    ResponseAmp = squeeze(sum(spike_c(21:60,TrialNum(i),1:56),1));
    ResponseAmp([23,49,53])=0;
    ResMap = zeros(8,9); 

    % correct spatial relationship 
    ResMap(1,[5,4,3]) = ResponseAmp(1:3); 
    ResMap(2,[4:7,3,2]) = ResponseAmp(4:9); 
    ResMap(3,[8:-1:2]) = ResponseAmp(10:16); 
    ResMap(4,[8:-1:2]) = ResponseAmp(17:23); 
    ResMap(5,[9:-1:2]) = ResponseAmp(24:31); 
    ResMap(6,[9:-1:1]) = ResponseAmp(32:40);
    ResMap(7,[8:-1:1]) = ResponseAmp(41:48);
    ResMap(8,[8:-1:1]) = ResponseAmp(49:56);

    ResMap_interp = interp2(ResMap,4,'makima');
    % figure, imagesc(ResMap);
    figure, 
    imagesc(linspace(0,2.8,size(ResMap_interp,2)), linspace(0,2.45,size(ResMap_interp,1)),ResMap_interp);
    set(gca,'YDir','normal')
    colormap('hot');
    caxis([0 150])
    xticks([0 2.8])
    yticks([0 2.45])
%     savefig(['_Mapping 56Ch normal rat',FileNum,'_Trial',num2str(TrialNum(i)),'.fig'])
end

% save([FilePath, '_Processed results', FileNum,'.mat'],'spike_p');
