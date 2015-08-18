%% Function to run the Michael Hills Detector over a dataset on the IEEG.org portal

function classified = Hills_Detector2(dataset,sz_layer)
%% Nested Functions

%% extractFeats
    % This function extracts the Time Domain and Frequency Domain Features
    % used in the Michael Hills Algorithm for a clip of data. 
    function feature_vector = extractFeats(data)
        
        % Time Domain
        % Prepare data for decimation
%         if size(data,2)>500
%             data = data(:,1:500);
%         end
        % Save data for use in frequency domain
        fdata = data;
        
        % Resample, normalize, and scale data by channel
        data3 = zeros(size(data,1),400);
        for r = 1:size(data,1)
            
            if std(data(r,:)) == 0
                data(r,:) = rand(1,size(data,2));
                fdata(r,:) = rand(1,size(fdata,2));
            end
            
            ddata = data(r,:);
%             ddata = decimate(data(r,:),5);
%             ddata = spline(1:100,ddata,.25:.25:100);
            ddata = (ddata-mean(ddata))/std(ddata-mean(ddata));
            data3(r,:) = ddata;
        end
        data = data3';
        
        % Calculate the upper triangle of the correlation coefficient
        % matrix and sorted eigen values.
        corrcoefs = corrcoef(data);
        eigen_values = sort(eig(corrcoefs));      
        corrcoefs = nonzeros(triu(corrcoefs))';

        % Frequency Domain
        % FFT magnitudes
        fdata = abs(fft(fdata))';
        
        % log10 of 1-47Hz
        fdata = fdata + 2e-13;
        fdata = log10(fdata(1:47,:));
        logdata = fdata;
        
        % Normalize and scale data
        for r = 1:size(fdata,1)
            ddata = fdata(r,:);
            ddata = (ddata-mean(ddata))/std(ddata-mean(ddata));
            fdata(r,:) = ddata;
        end
        fdata = fdata';
        
        % Compute correlation coefficients and sorted eigen values
        fcorrcoefs = corrcoef(fdata);
        feigen_values = sort(eig(fcorrcoefs));
        fcorrcoefs = nonzeros(triu(fcorrcoefs))';
        
        % Concatenate all features into a row vector
        feature_vector = [corrcoefs eigen_values' reshape(logdata,1,(size(logdata,1)*size(logdata,2))) fcorrcoefs feigen_values'];
    end
    
    %% clip_ictal_segments
    % Function to clip ictal segments into 1-second clips
    function feature_matrix = clip_ictal_segments(ictalSeg,fs)
        feature_matrix = [];
        % For each segment
        for i = 1:length(ictalSeg);
            
            % number of clips in segment
            nclips = floor(length(ictalSeg{i})/fs);
            
            % clipping
            for j = 1:nclips
                data = ictalSeg{i}(fs*(j-1)+1:fs*j,:)';
                feature_matrix_tmp(j,:) = extractFeats(data);
            end
            feature_matrix = [feature_matrix; feature_matrix_tmp];
            feature_matrix_tmp = [];

        end
    end

    %% clip_intertrain_segments
    % Function to clip interictal segments into 1-second clips
    function feature_matrix = clip_intertrain_segments(interSeg,fs)
        feature_matrix = [];
        % for each segment
        for i = 1:length(interSeg);
            
            % number of clips in segment
            nclips(i) = floor(length(interSeg{i})/fs);
            
            % randomly clip for half the number of clips
            ran = randperm(nclips(i));
            for j = 1:round(nclips(i)/2)
                data = interSeg{i}(fs*(ran(j)-1)+1:fs*ran(j),:)';
                feature_matrix_tmp(j,:) = extractFeats(data);
            end
            feature_matrix = [feature_matrix; feature_matrix_tmp];
            feature_matrix_tmp = [];

        end
    end

%% Start Main

%% Training

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
y= 'y';
n= 'n';

prompt = sprintf('There are %d annotated seizures in the dataset. How many should be used for training? (>1): ',nsz);
ntrain = input(prompt);

% Pull ictal segments from IEEG.org
ictal_segments = cell(1,ntrain);

% Seizure duration to train on when no end time is annotated
secSz = 15;

% Set offset from EEC if desired
offset = 0*1e6;

% Seconds before annotated seizure to classify as ictal for training
preSz = 0;

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
        
    inter_segments{i} = dataset.getvalues(inter_points(i)-2.5*avgIctalTime,5*avgIctalTime,1:length(dataset.channels));
    
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


fprintf('Clipping ictal segments and extracting features...')

% clip ictal segments
data1 = clip_ictal_segments(ictal_segments,fs);

% create labels
labels = ones(1,size(data1,1));

fprintf('Done!\n')

fprintf('Clipping interictal segments and extracting features...')

% clip interictal segments
data2 = clip_intertrain_segments(inter_segments,fs);

% concatenate ictal and interictal feature matrices
training_matrix = [data1; data2];

% concatenate labels
labels = [labels zeros(1,size(data2,1))];

clear data data1 data2 inter_segments ictal_segments dummy bchans1 bchans2 bchans

fprintf('Done!\n')

%% Train Model
fprintf('Training model...')

try
    delete(gcp)
catch
end

mypool = parpool(2);

paroptions = statset('UseParallel',true);

% start clock
tic

% train random forest of 3000 trees
model = TreeBagger(3000,training_matrix,labels','NVarToSample',1,'Options',paroptions);

fprintf('Done!  ')

% stop clock
toc

fprintf('\n')

% clear data to free up space
clear training_matrix

%% Pull testing data, detect, and upload detections to the portal
custom = 0;
a = 0;
fprintf('Pulling testing data from IEEG.org...\n')

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

    fprintf('Clipping test set and extracting features...')

    % clip segments
    for k = 1:length(ictal_segments)

        test_matrix_sz{k} = clip_ictal_segments({[inter_segments{k}; ictal_segments{k}]},fs);
        szSize(k) = size(test_matrix_sz{k},1);

    end

    % convert test_matrix_sz cell array to matrix
    test_matrix = cell2mat(test_matrix_sz');

    fprintf('Done!\n')

end

fprintf('Classifying testing set...')

%% Classify clips

% make predictions
classified = predict(model,test_matrix);
classified = cellfun(@str2num,classified);

% create times to upload detections to IEEG.org
times = [];

if custom == 'y'
    times = linspace(start_time,stop_time,size(test_matrix,2));
else
    for i = 1:length(start_times)-ntrain
        
        if check == 1
            times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain}+secSz*1e6,szSize(i))];
        else
            times = [times linspace(start_times{i+ntrain}-5*avgIctalTime,stop_times{i+ntrain},szSize(i))];
        end
    end
end
times = times(classified == 1);
fprintf('Done!\n')

%% Upload detections

fprintf('Uploading detections to IEEG.org...')
layerName = 'Hills_Function';
eventTimesUSec = times';
eventChannels = cell(1,length(times));
for i = 1:length(eventChannels)
    eventChannels{i} = 1:length(dataset.channels);
end
label = 'Seizure';
uploadAnnotations_v2(dataset,layerName,eventTimesUSec,eventChannels,label)



end
