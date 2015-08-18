%% Matlab function to interface the Hills python algorithm with the portal
function [sens,spec,acc] = Hills_Detector_pystudy(dataset,sz_layer,tnsz,tsz,offset)
%% Nested Functions

%% Clipping
%% clip_ictal_segments
% Function to clip ictal segments into 1-second clips
    function  clip_ictal_segments(ictalSeg,freq,channels)
        
        count = 0;
        
        % For each segment
        for seg = 1:length(ictalSeg);
            ictalSeg{seg} = single(ictalSeg{seg});
            % number of clips in segment
            nclips = floor(length(ictalSeg{seg})/freq);
            
            % clipping
            for clip = 1:nclips
                count = count + 1;
                data = ictalSeg{seg}(freq*(clip-1)+1:freq*clip,:)';
                for chan = 1:size(data,1)
                    if std(data(chan,:)) == 0
                        data(chan,:) = rand(1,size(data,2));
                    end
                end
                % save clip
                latency = clip-1;
                str = sprintf('Hills/seizure-detection-master/seizure-data/Patient_1/Patient_1_ictal_segment_%d.mat',count);
                save(str,'channels','data','freq','latency')
            end
        end
    end

%% clip_inter_segments
% Function to clip interictal segments into 1-second clips
    function clip_inter_segments(interSeg,freq,channels)
        count = 0;
        
        % for each segment
        for seg = 1:length(interSeg);
            interSeg{seg} = single(interSeg{seg});
            % number of clips in segment
            nclips(seg) = floor(length(interSeg{seg})/freq);
            
            % randomly clip for half the number of clips
            ran = randperm(nclips(seg));
            
            for clip = 1:round(nclips(seg)/2)
                count = count +1;
                data = interSeg{seg}(freq*(ran(clip)-1)+1:freq*ran(clip),:)';
                
                for chan = 1:size(data,1)
                    if std(data(chan,:)) == 0
                        data(chan,:) = rand(1,size(data,2));
                    end
                end
               % save clip
                str = sprintf('Hills/seizure-detection-master/seizure-data/Patient_1/Patient_1_interictal_segment_%d.mat',count);
                save(str,'channels','data','freq')
            end
        end
    end
%% clip_test_segments
% Function to clip test segments into 1-second clips
    function  szSize = clip_test_segments(testSeg,freq,channels)
        count = 0;
        
        % For each segment
        for seg = 1:length(testSeg);
            testSeg{seg} = single(testSeg{seg});
            
            % number of clips in segment
            nclips = floor(length(testSeg{seg})/freq);
            szSize(seg) = nclips;
            
            % clipping
            for clip = 1:nclips
                count = count + 1;
                data = testSeg{seg}(freq*(clip-1)+1:freq*clip,:)';
                
                for chan = 1:size(data,1)
                    if std(data(chan,:)) == 0
                        data(chan,:) = rand(1,size(data,2));
                    end
                end
                % save clip
                str = sprintf('Hills/seizure-detection-master/seizure-data/Patient_1/Patient_1_test_segment_%d.mat',count);
                save(str,'channels','data','freq')
            end
        end
    end

%% Evaluation

seni = @(x) find(x == 1);
speci = @(x) find(x == 0);

sensitivity = @(x,y) length(find(x(seni(y)) == 1))/length(seni(y));
specificity = @(x,y) length(find(x(speci(y)) == 0))/length(speci(y));

accuracy = @(x,y) length(find(x==y))/length(x);

%% Start Main

%% Pull ictal/interictal segments from the portal

fprintf('Pulling training data from IEEG.org...\n')

% Sampling rate
fs = round(dataset.sampleRate);

% Annotation times of specified layer (sz_layer)
layerNames = {dataset.annLayer.name};
ind=find(ismember(layerNames,sz_layer));
seizures = dataset.annLayer(ind).getEvents(0);

% Seizure start and stop times (ictal period)
start_times = {seizures.start};
stop_times = {seizures.stop};

% number of seizures
nsz = length(start_times);

% initialize yes and no responses
y = 'y';
n = 'n';

% % ask how many seizures should be trained on
% prompt = sprintf('There are %d annotated seizures in the dataset. How many should be used for training? (>1): ',nsz);
ntrain = tnsz;

% ***************************** Pull ictal segments from IEEG.org
ictal_segments = cell(1,ntrain);

% Seizure duration to train on when no end time is annotated
secSz = tsz;

% Set offset from EEC if desired
offset = offset*1e6;

% Seconds before annotated seizure to classify as ictal for training
preSz = 5;

check = 0;

for i = 1:ntrain
    
    if stop_times{i} == start_times{i}
        ictal_segments{i} = dataset.getvalues(start_times{i}-preSz*1e6+offset,(secSz+preSz)*1e6,1:length(dataset.channels));
        check =1;
    else
    
    % Pull data between EEC and seizure end
    ictal_segments{i} = dataset.getvalues(start_times{i},stop_times{i}-start_times{i},1:length(dataset.channels));    
    end
    
    % Change NaNs and inf to zeros and then add noise to dropped channels
   ictal_segments{i}(isnan(ictal_segments{i}))=0;
   ictal_segments{i}(isinf(ictal_segments{i}))=0;
   
   for k = 1:size(ictal_segments{i},2)
        if std(ictal_segments{i}(:,k)) == 0
            ictal_segments{i}(:,k) = rand(size(ictal_segments{i}(:,k)));
        end
   end

end

% Pull interictal data the same length as ictal data
inter_segments = cell(1,ntrain);

if check == 1
    avgIctalTime = 1e6*(secSz+preSz);
else
% Average ictal time
    avgIctalTime = mean(cell2mat(stop_times(1:ntrain))-cell2mat(start_times(1:ntrain)));
end

% midpoints between seizure end and seizure beginning (interictal
% period)
inter_points = diff([0 start_times{1:ntrain}])/2 + [0 stop_times{1:ntrain-1}];

% Pull interictal segments from IEEG.org
for i = 1:ntrain
        
    inter_segments{i} = dataset.getvalues(start_times{i+1}-5*avgIctalTime,5*avgIctalTime,1:length(dataset.channels));
    
    inter_segments{i}(isnan(inter_segments{i}))=0;
    inter_segments{i}(isinf(inter_segments{i}))=0;
    
    if sum(std(inter_segments{i}))==0
        inter_segments{i} = rand(size(inter_segments{i}));
    end
    
end


% channels
dummy = {dataset.channels.label};
dummy = strrep(dummy,' ','_');
dummy = strrep(dummy,'-','_');

for i = 1:length(dummy)
    if i == 1
        channels = struct(dummy{1},dummy{1});
    else
        channels=setfield(channels,{},dummy{i},dummy{i});
    end
end



fprintf('Done!\n')

fprintf('Clipping training segments...')

% Delete old data files and add path to new files
system('rm -rf ~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1');
system('mkdir ~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1');
addpath(genpath('~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1'));

% Clip training segments into 1-second .mat files
clip_ictal_segments(ictal_segments,fs,channels)
clip_inter_segments(inter_segments,fs,channels)

fprintf('Done!\n')

%% Pull testing data from the portal

custom = 0;
a = 0;
fprintf('Pulling testing data from IEEG.org...\n')

custom = n;

% Ask user what they would like to test on
while (custom == 'y' | custom == 'n') == 0
    if a < 1
        custom = input('Would you like to detect over a custom interval? (y/n): ');
    else
        fprintf('Please answer (y/n)\n')
        custom = input('Would you like to detect over a custom interval? (y/n): ');
    end
    a = a+1;
end

% If user would like to test on a custom data segment, collect input and
% pull
if custom == 'y'
    
    start_time = input('Please enter the start time in \mus: ');
    stop_time = input('Please enter the end time in \mus: ');
    test_segment{1} = dataset.getvalues(start_time,stop_time-start_time,1:length(dataset.channels));
    
    test_segment{1}(isnan(test_segment{1})) = 0;
    test_segment{1}(isinf(test_segment{1})) = 0;
    
    fprintf('Done!\n')
    fprintf('Clipping custom test segment...')
    clip_test_segments(test_segment,fs,channels);
    fprintf('Done!\n')
    
else

    % Pull testing ictal segments from IEEG.org
    ictal_segments = cell(1,length(start_times)-ntrain);
    
    for i = ntrain+1:length(start_times)
        j = i-ntrain;
        if stop_times{i} == start_times{i}
            ictal_segments{j} = dataset.getvalues(start_times{i},secSz*1e6,1:length(dataset.channels));
            check =1;
        else
            % Pull data between EEC and seizure end
            ictal_segments{j} = dataset.getvalues(start_times{i},stop_times{i}-start_times{i},1:length(dataset.channels));   
        end
        
        ictal_segments{j}(isnan(ictal_segments{j})) = 0;
        ictal_segments{j}(isinf(ictal_segments{j})) = 0;
       
        
    end

    % Pull interictal data the same length as ictal data
    inter_segments = cell(1,length(start_times)-ntrain);

    for i = ntrain+1:length(start_times)
        j = i-ntrain;
        inter_segments{j} = dataset.getvalues(start_times{i}-5*avgIctalTime,5*avgIctalTime,1:length(dataset.channels));
        
        inter_segments{j}(isnan(inter_segments{j})) = 0;
        inter_segments{j}(isinf(inter_segments{j})) = 0;
        
    end
    

    fprintf('Done!\n')

    fprintf('Clipping testing segment(s)...')

    
    % Clip into 1-second .mat files
    segments = cell(1,length(inter_segments));
    nictal = zeros(1,length(inter_segments));
    ntotal = zeros(1,length(inter_segments));
    
    for i = 1:length(inter_segments)
        nictal(i) = floor(length(ictal_segments{i})/fs);
        segments{i} = [inter_segments{i}; ictal_segments{i}];
        ntotal(i) = floor(length(segments{i})/fs);
    end
    
    szSize = clip_test_segments(segments,fs,channels);
    
    labels = [];
    for i = 1:length(szSize)
        
        tlabels = ones(1, szSize(i));
        tlabels(1:ntotal(i)-nictal(i)) = 0;
        labels = [labels tlabels];
        
    end
        
      

end

fprintf('Done!\n')

% Prepare to run python algorithm
PATH_PYTHON = '/Users/Tyler/anaconda/lib/python2.7/site-packages';
setenv('PYTHONPATH', PATH_PYTHON); % set env path (PYTHONPATH) for this session

fprintf('Running predict.py...')

% Run python algorithm
! sh /Users/Tyler/Documents/MATLAB/Thesis/Hills/hillsalgo.sh

fprintf('Done!\n')

% load patient data
Hillsprob = csvread('Hills/seizure-detection-master/submissions/submission.csv',1,1);

% Seizure Detections
Hsz = round(Hillsprob(:,1));
times = [];
if custom == 'y'
    times = linspace(start_time,stop_time,length(Hsz));
else
    for i = 1:length(start_times)-ntrain
        if check == 1
            times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain}+secSz*1e6,szSize(i))];
        else
            times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain},szSize(i))];
        end
    end
end

times = times(Hsz == 1);

sens = sensitivity(Hsz,labels');
spec = specificity(Hsz,labels');
acc = accuracy(Hsz,labels');

eval = sprintf('The detection sensitivity is %f, the specificity is %f, the accuracy is %f\n',sens,spec,acc);

fprintf(eval);

% %% Upload detections
% 
% fprintf('Uploading detections to IEEG.org...')
% layerName = 'Hills_Detector_py';
% eventTimesUSec = times';
% eventChannels = cell(1,length(times));
% for i = 1:length(eventChannels)
%     eventChannels{i} = 1:length(dataset.channels);
% end
% label = 'Seizure';
% uploadAnnotations_v3(dataset,layerName,eventTimesUSec,eventChannels,label)
% 
% fprintf('Done!\n')
% addit = 'y';
% 
% while addit == 'y'
% addit = input('Would you like to classify another segment (y/n)?: ');
% if addit == 'n'
%     break
% end
% 
% system('rm -rf ~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1');
% system('mkdir ~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1');
% addpath(genpath('~/Documents/MATLAB/Thesis/Hills/seizure-detection-master/seizure-data/Patient_1'));
% 
% start_time = input('Please enter the start time in \mus: ');
% stop_time = input('Please enter the end time in \mus: ');
% 
% fprintf('Pulling data from IEEG.org...')
% test_segment = {dataset.getvalues(start_time,stop_time-start_time,1:length(dataset.channels))};
%     test_segment(isnan(test_segment)) = 0;
%     test_segment(isinf(test_segment)) = 0;
% fprintf('Done!\n')
% 
% fprintf('Clipping custom test segment...')
% szSize = clip_test_segments(test_segment,fs,channels);
% fprintf('Done!\n')
% 
% fprintf('Running predict.py...')
% % Run python algorithm
% ! sh /Users/Tyler/Documents/MATLAB/Thesis/Hills/hills_predict.sh
% 
% fprintf('Done!\n')
% 
% % load patient data
% Hillsprob = csvread('Hills/seizure-detection-master/submissions/submission.csv',1,1);
% 
% % Seizure Detections
% Hsz = round(Hillsprob(:,1));
% 
% if custom == 'y'
%     times = linspace(start_time,stop_time,length(Hsz));
% else
%     for i = 1:length(start_times)-ntrain
%         if check == 1
%             times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain}+secSz*1e6,szSize(i))];
%         else
%             times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain},szSize(i))];
%         end
%     end
% end
% 
% times = times(Hsz == 1);
% 
% %% Upload detections
% 
% fprintf('Uploading detections to IEEG.org...')
% layerName = 'Hills_Detector_py';
% eventTimesUSec = times';
% eventChannels = cell(1,length(times));
% for i = 1:length(eventChannels)
%     eventChannels{i} = 1:length(dataset.channels);
% end
% label = 'Seizure';
% uploadAnnotations_v3(dataset,layerName,eventTimesUSec,eventChannels,label)
% end
end
