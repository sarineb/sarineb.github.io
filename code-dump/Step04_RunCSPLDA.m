function PublishOutput = Step04_RunCSPLDA(PublishInput)
% Publishing Methods are .m functions that look at result files and produce new result files or publication figures
% PublishInput
%     .ResultFiles - cell array of file names passed to the method
%     .HelperFunctions - cell array of helper function names to be used with EnableHelperFunction
%     .PathSep - path separator '/' in mac/linux and '\' in windows
%     .ProgramStatusIndicator - a handle to the Program Status Indicator on the front panel so that status can be updated from within the reader
% PublishOutput - not currently used
%
% This is the main CSP analysis function. It takes as input the task repetitions, and runs a leave-one-out cross validation to classify force repetitions using the CSP + LDA method. The main method was inspired by EEGLAB?s CSP for BCI, see http://sccn.ucsd.edu/wiki/Minimalist_BCI
% Input: either the output from the last analysis step, or the output from Step03_RunICAeeglab_RemoveComponents.m for ICA-pruned data.
% Output: the output is saved as a mat file and the link to the result is also stored in a datatable. Results for LME purposes also stored in DataTableLME.mat


for i = 1:20
    try
        close(i)
    end
end

BaseDirectory = cd;

try
    
    set(PublishInput.ProgramStatusIndicator,'string','Step04_RunCSPLDA: Running...','foregroundcolor','g');drawnow;
    
    %----------- Your code here !!! ---------------
    file2str = EnableHelperFunction([],'file2str.m');
    
    AnalysisName = 'CSP_UCM';
    
    NumTrials = length(PublishInput.ResultFiles);
    
    path1 = '/Users/amplmember/Google Drive/2015_AMPL_Sarine_UCM/UCM_AnalysisFunctions/';
    path2 = '/Volumes/amplmember/Documents/Experiments/2015_AMPL_Sarine_UCM/2015_AMPL_Sarine_UCM_Analysis/';
    
    % decide on method:
    CSPMethod = 'geneig'; % 'geneig' for generalized eigenvalue method or 'whitening' for whitening transform method
    
    % define/initiate common variables:
    Fs = 100; % to downsample to 100 hz
    EpochLocation = 'last'; % 'first', 'middle', or 'last'
    
    % define which period to analyze:
    DataType = 'prep'; % 'prep' for preparatory period or 'hold' for movement period, or 'rest' for rest period
    
    % decide whether to keep all data or remove middle portion:
    ExcludeMidThird = 'yes'; % 'yes' or 'no', to only keep first and last thirds of data. Only use this for UCM data; use 'no' for two-force data.
    div = 3; % 3 for thirds, use 4 to keep less data...
    
    % define a filter whose frequency response is specified by a function
    % (here: a smooth sine-based function with peaks at 12 and 25 Hz)
    Fmin = 7;
    %     Fmax = 14;
    %     Fmin = 15;
    Fmax = 30;
    flt = @(f)(f>Fmin&f<Fmax).*(1-cos((f-(Fmin+Fmax)/2)/(Fmin-Fmax)*pi*4));
    %     flt = @(f)(f>Fmin&f<Fmax);
    f = flt;
    
    nof = 2; % # filters to use from each end (total = 2*nof)
    n = 200; %  for buffer in 
    
    %---- define channels to use:
    %     channels = [5:7 9:12 15:17 20:31 38 39 41:57 62:64]; % exclude mastoids, T's and some of the frontal electrodes and EOG (44 electrodes)
    % %     channels = [9 10 11 12 15 16 17 20 21 22 23 41 42 43 44 45 46 47 48 49]; % sensorimotor area (20 electrodes)
    channels = [5 6 7 9 10 11 12 15 16 17 20 21 22 23 37 38 39 40 41 42 43 44 45 46 47 48 49]; % sensorimotor area with some frontal (27 electrodes)
    %     channels = [1:12 14:18 20:31 33:64]; % exclude only mastoids and EOG (61 electrodes)
    
    DataTableLME = [];
    
    %---- loop to batch process all selected participants:
    for TrialNum = 1:NumTrials
        clear RepsDS RepsEpoch Reps IndEpoch RepsFLT yy y CrVal S V b w
        
        %---- load file for participant, get subject code, get repetitions
        %in EEG and force, etc.
        load([PublishInput.ResultFiles{TrialNum}]); % load input file
        
        FileNameInd = findstr(PublishInput.ResultFiles{TrialNum},'/');
        FileNameInd = FileNameInd(end);
        FileName = PublishInput.ResultFiles{TrialNum}((FileNameInd+1):end-4); % keep file name
        
        FsUS = round(Result.FT_CMC_IN.CMC{1}.Fs); % original (high) sampling rate
        SubjectCode = FileName(1:4); % subject name
        
        %---- get all repetitions and UCM force values:
        Reps = cellfun(@(x) x(channels,:), Result.FT_CMC_IN.CMC{1}.FTEEGData, 'uni', 0);
        if strcmp(DataType, 'prep')
            Reps = cellfun(@(x) x(channels,:), Result.FT_CMC_IN.CMC{1}.Prep.FTEEGData, 'uni', 0);
        elseif strcmp(DataType, 'rest')
            Reps = cellfun(@(x) x(channels,:), Result.FT_CMC_IN.CMC{1}.Rest.FTEEGData, 'uni', 0);
        end
        
        FUCM = Result.ForceData.ForcesUCM;
        
        %---- set duration of epoch to be analyzed (in seconds)
        EpochDuration = 0.5; % in seconds
%         %-- for prep data:
%         if strcmp(DataType, 'prep')
%             EpochDuration = 0.5;
%         end
        
        %--- get the window of the epoch to be analyzed, in samples:
         wnd = [0 EpochDuration]; wnd = round(Fs*wnd(1)):round(Fs*wnd(2)); wnd = wnd(2:end); % this is the actual window used based on epochduration
        
        
        %---- frequency filtering of EEG and temporal filter estimation:
        EEG = cell2mat(Reps)'; % change from cell to matrix form
        [t,c] = size(EEG); idx = reshape(1:t*c-mod(t*c,n),n,[]);
        FLT = real(ifft(fft(EEG).*repmat(f(Fs*(0:t-1)/t)',1,c)));
        T = FLT(idx)/EEG(idx);
        
        %---- change the filtered data back into cells:
        StartInd = 1;
        for r = 1:length(Reps)
            EndInd = StartInd + size(Reps{r}',1) - 1;
            RepsFLT{r} = FLT(StartInd:EndInd,:)';
            StartInd = EndInd + 1;
        end
        
        
        %---- assign labels based on force output (true labels):
        AllLabels = ones(1,length(Reps)); % initialize labels to 1 for high ucm
        AllLabels(find(FUCM<0)) = -1; % assign -1 to low ucm
        
        %---- downsample the filtered reps -- we don't need 2048 hz for this:
        RepsDS = cellfun(@(x) resample(x',Fs,FsUS,0), RepsFLT, 'uni', 0);
        
        %---- find indices of and extract first, middle or last epoch of specified length:
        if strcmp(EpochLocation,'first')
            IndEpoch = cellfun(@(x) 1:round(EpochDuration*Fs),RepsDS,'uni',0);
        elseif strcmp(EpochLocation,'middle')
            IndEpoch = cellfun(@(x) round((length(x)-EpochDuration*Fs)/2)+1:EpochDuration*Fs+round((length(x)-EpochDuration*Fs)/2),RepsDS,'uni',0);
        elseif strcmp(EpochLocation,'last')
            IndEpoch = cellfun(@(x) (length(x)-round(EpochDuration*Fs)+1):length(x),RepsDS,'uni',0);
        end
        
        RepsEpoch = cellfun(@(x,y) x(y,:), RepsDS, IndEpoch, 'uni', 0);
        
        %---- group reps acc. to class (for training purposes):
        LowUCM = RepsEpoch(find(AllLabels == -1));
        HighUCM = RepsEpoch(find(AllLabels == 1));
        
        %---- remove middle portion of data (refer to 'ExcludeMidThird' and 'div' above for info on this):
        [UCMSort indSort] = sort(FUCM); % sort by ucm value
        LowSortInd = 1:length(UCMSort)/div; % get indices of sorted data belonging to first group
        HighSortInd = ceil(length(UCMSort)-length(UCMSort)/div+1):length(UCMSort); % get indices of sorted data belonging to second group
        %         SortLow = UCMSort(LowSortInd);
        %         SortMid = UCMSort(length(UCMSort)/div+1:length(UCMSort)/div*2);
        %         SortHigh = UCMSort(ceil(length(UCMSort)/div*2+1):end);
        %         SortHigh = UCMSort(HighSortInd);
        
        IndLow = indSort(LowSortInd);
        IndHigh = indSort(HighSortInd);
        if strcmp(ExcludeMidThird, 'yes')
            RepsEpoch = RepsEpoch([IndLow; IndHigh]); % remove the mid portion of reps
        end
        
        %---- cross-validation set up----
        TestCross = 1; % number of reps to use for testing - (not relevant here because it is leave-one-out cv, defaulted to 1)
        NumCross = length(RepsEpoch); % number of cross-validations
        AllInd = 1:NumCross; % for cv loop
        if NumCross < length(Reps) % if mid portion has been removed
            FUCM = [FUCM(IndLow); FUCM(IndHigh)]; % also remove mid portion of FUCM force data
            AllLabels = [AllLabels(IndLow) AllLabels(IndHigh)]; % same with labels
            %             AllLabels = [-1*ones(size(IndLow))' 1*ones(size(IndHigh))'];
        end
        
        %--cross-validation loop:
        cc = 1; % initialize counter
        for ss = AllInd
            EPO = []; X = []; y = []; yy = []; % initialize per cross-validation 
            
            %--- divide into training data and test data:
            TestInd = ss;
            TrainInd = setdiff(AllInd,TestInd);
            TrainData = RepsEpoch(TrainInd);
            TrainFUCM = FUCM(TrainInd);
            TestData = RepsEpoch(TestInd);
            TestFUCM = FUCM(TestInd);
            
            %---- extract the corresponding labels:
            TrainLabels = AllLabels(TrainInd);
            TestLabels = AllLabels(TestInd);
            
            %---- group reps acc. to class:
            LowUCMTrain = TrainData(find(TrainLabels == -1));
            HighUCMTrain = TrainData(find(TrainLabels == 1));
            
            %---- prepare data for CSP:
            EPO{1} = cell2mat(LowUCMTrain');
            EPO{2} = cell2mat(HighUCMTrain');
            
            
            if strcmp(CSPMethod, 'geneig')
                %---- CSP training and filters ----
                % set up the generalized eigenvalue problem:
                [V{cc},D] = eig(cov(EPO{2}),cov(EPO{1})+cov(EPO{2})); % V is the entire Spatial filters matrix (columns); the inverse of its trsfm is the CSP's (columns)
                % extract the few most important filters from each side
                S{cc} = V{cc}(:,[1:nof end-nof+1:end]); % these are the 2*nof Spatial filters [channels x 2*nof]
                
            elseif strcmp(CSPMethod, 'whitening')
                %---- Whitening method:
                TrainReps{1} = LowUCMTrain;
                TrainReps{2} = HighUCMTrain;
                %---- find normalized and averaged covariance matrices:
                for kk = 1:length(TrainReps)
                    X{kk} = cellfun(@(x) x', TrainReps{kk}, 'uni', 0); % transpose the reps to have [channels x samples]
                    C{kk} = cellfun(@(x) x*x'/trace(x*x'), X{kk}, 'uni', 0); % find covariance matrix for each rep
                    Ccat{kk} = cat(3, C{kk}{:}); % concatenate into 3D array to find mean per-element over all trials
                    Cbar{kk} = mean(Ccat{kk},3); % take the mean over all trials (average cov mtx)
                end
                Cc = Cbar{1} + Cbar{2}; % composite spatial covariance
                
                %---- find whitening transform:
                [Uc Lambdac] = eig(Cc);
                [Lambdac,ind] = sort(diag(Lambdac),'descend'); % sort eig. values in desc. order
                Uc=Uc(:,ind); % also sort eig. vectors according to the above
                P = sqrt(inv(diag(Lambdac)))*Uc'; % whitening trsfm
                for kk = 1:length(X)
                    Ss{kk} = P*Cbar{kk}*P'; % whiten average covariance matrices
                end
                
                [B Lambdab] = eig(Ss{1}, Ss{2}); % should be equivalent to eig(Ss{1}) and eig(Ss{2})
                [Lambdab,ind] = sort(diag(Lambdab)); % sort ascending
                B=B(:,ind);
                
                % ---- find spatial filters:
                V{cc}=(B'*P);
                for i=1:length(ind), V{cc}(i,:)=V{cc}(i,:)./norm(V{cc}(i,:)); end % not sure if/why this is needed
                S{cc} = [V{cc}(:,1:nof) ,V{cc}(:,end-nof+1:end)];
            end
            
            %---- get Common Spatial Patterns ----
            A = inv(V{cc}'); % the columns are the patterns
            ALow = A(:,1);
            AHigh = A(:,end);
            
            %---- log-variance feature extraction and LDA classifier weights:
            for k = 1:2
                X{k} = squeeze(log(var(reshape(EPO{k}*S{cc}, length(wnd),[],2*nof))));
            end
            w{cc} = ((mean(X{2})-mean(X{1}))/(cov(X{1})+cov(X{2})))'; % these are beta_1...beta_2nof
            b{cc} = (mean(X{1})+mean(X{2}))*w{cc}/2; % this is -beta_0
            
            
            %---- classification of test data ----
            for tt = 1:length(TestData) % this room only runs once if LOO-CV
                TestRep = TestData{tt};
                yy(tt) = test_bci(TestRep,S{cc},T,w{cc},b{cc}); % classifier output
            end
            
            %---- re-classify training data (for normalization in ROC curve)
            for ntr =  1:length(TrainData)
                TestRep = TrainData{ntr};
                ytr{ss}(ntr) = test_bci(TestRep,S{cc},T,w{cc},b{cc});
            end
            trainlabels{ss} = TrainLabels;
            testlabels{ss} = TestLabels;
            %------
            
            %---- store classification and labels:
            CrVal.y{cc} = yy;
            CrVal.TrueLabel{cc} = TestLabels;
            CrVal.TestInd{cc} = TestInd;
            CrVal.TrainInd{cc} = TrainInd;
            CrVal.TestFUCM{cc} = TestFUCM;
            % CrVal.TestFUCM{cc} = TrainFUCM;
            
            cc = cc+1;
        end
        
        
        %---- plot ROC curves of training data at each leave one out cross
        %validation + find std and mean of training classification outpu to normalize test output:
        figure; hold on; for s = 1:length(ytr)
            outputs=ytr{s}; labels = trainlabels{s};
            outputs(labels==-1) = ytr{s}(labels==-1)-mean(ytr{s}(labels==-1));
            outputs(labels==1) = ytr{s}(labels==1)-mean(ytr{s}(labels==1));
            stdev{s}=std(outputs);
            mean1{s}=mean(ytr{s}(labels==-1));
            outputsn{s} = (ytr{s}-mean1{s})/stdev{s};
            [X,Y,T,AUC(s)] = perfcurve(trainlabels{s},outputsn{s},1);
            plot(X,Y)
        end
        
        % find ROC curve and compute AUC (built-in matlab function)
        [Xt,Yt,Tt,AUCt{TrialNum}] = perfcurve(cell2mat(trainlabels),cell2mat(outputsn),1);
        
        hold on; plot(Xt,Yt,'linewidth',3,'color','k') % plot overall roc curve
        xlabel('false positive rate'); ylabel('true positive rate'); title(['ROC curves of training data of ' SubjectCode '. Overall AUC = ' num2str(AUCt{TrialNum})])
        %----
        yyn{TrialNum} = (cell2mat(CrVal.y)-cell2mat(mean1))./cell2mat(stdev); % normalized output from LOO-cross validation test
        clear outputs stdev mean1 outputsn trainlabels testlabels ytr
        
        %---- plotting results ----
        yval = cell2mat(CrVal.y);
        yplot = yval;
        yplot = yplot/sqrt(mean(yplot.*yplot));
        
        figure('color','w'); plot(cell2mat(CrVal.TestFUCM), yplot, '.', 'markersize', 16);
        xlabel('F_{UCM}', 'fontsize', 14); ylabel('Classification', 'fontsize', 14); title(SubjectCode, 'fontsize', 14)
        
        %----
        % find percentage of correct classifications - this isn't accurate
        % because of bias issues (beta_0 value). refer to ROC curve accuracy for a
        % more normalized comparison between participants.
        Correct = [];
        for cc = 1:length(CrVal.y)
            yy =yplot(cc);
            TestFUCM = CrVal.TestFUCM{cc};
            if (TestFUCM<0 && yy<0) || (TestFUCM>0 && yy>0)
                Correct(cc) = 1;
            else
                Correct(cc) = 0;
            end
        end
        
        Subject.PercentCorrect{TrialNum} = sum(Correct)/length(Correct)*100
        %-----
        
        [sfit,good] = fit(cell2mat(CrVal.TestFUCM)',yplot','a*x+b'); % fit a line to the predictions from low to high ucm
        slopeucm(TrialNum) = sfit.a;
        hold on; plot(cell2mat(CrVal.TestFUCM)',sfit(cell2mat(CrVal.TestFUCM)'),'r')
           
        ForceData = Result.ForceData;
        
        %---- store things into structure:
        clear Result;
        Result.CodeUsed = file2str([mfilename('fullpath') '.m']);
        Result.CSP.W = V;
        Result.CSP.Wj = S;
        Result.CSP.NFilters = nof;
        Result.LDA = [b; w];
        Result.BPFreq = [Fmin Fmax];
        Result.CSP.EpochLocation = EpochLocation;
        Result.CSP.EpochDuration = EpochDuration;
        Result.CSP.ChannelsUsed = channels;
        Result.Fs = Fs;
        Result.CrVal = CrVal;
        Result.Labels = {'-1 = Low UCM'; '1 = High UCM'};
        Result.ForceData = ForceData;
        
        %--- save things:
        path1 = '/Users/amplmember/Google Drive/2015_AMPL_Sarine_UCM/UCM_AnalysisFunctions/';
        path2 = '/Volumes/amplmember/Documents/Experiments/2015_AMPL_Sarine_UCM/2015_AMPL_Sarine_UCM_Analysis/';
        if strcmp(FileName(end-2:end),'ICA')
            ForceNum = FileName(end-4);
            EndVal = [ForceNum '_' FileName(13:end-16) '_' EpochLocation '_' DataType '_ICA']
        else
            EndVal = FileName(end);
            ForceNum = EndVal;
            EndVal = [ForceNum '_' EpochLocation '_' DataType];
        end
        
        SaveFileName = [path2 SubjectCode '_CSP_' num2str(length(channels)) 'ch_' num2str(NumCross) 'reps_' EndVal '.mat'];
        save(SaveFileName,'Result');
        
        %-----------------------
        
        %-----also add results to data table----
        DataTablePath = [path2 'DataTable' num2str(length(channels)) 'Electrodes_' EndVal(3:end)]; % define datatable path
        load(DataTablePath)
        
        % determine the row of current participant:
        RowNum = find(strcmp(DataTable{:,1},[SubjectCode 'S001']));
        if strcmp(ForceNum,'1')
            ForceLevel = 'ThreeN';
        elseif strcmp(ForceNum,'2')
            ForceLevel = 'SixN';
        elseif strcmp(ForceNum,'3')
            ForceLevel = 'SixNVol';
        end
        
        if strcmp(ExcludeMidThird,'yes')
            SlopeColName = ['Slope' ForceLevel 'thirds'];
            ResultColName = [ForceLevel 'thirds'];
        elseif (strcmp(ExcludeMidThird,'no') && ~strcmp(ForceNum,'3'))
            SlopeColName = ['Slope' ForceLevel 'full'];
            ResultColName = [ForceLevel 'full'];
        elseif (strcmp(ExcludeMidThird,'no') && strcmp(ForceNum,'3'))
            SlopeColName = ['Slope' ForceLevel];
            ResultColName = [ForceLevel];
        end
        
        SlopeColNum = find(strcmp(DataTable.Properties.VariableNames,SlopeColName));
        ResultColNum = find(strcmp(DataTable.Properties.VariableNames,ResultColName));
        
        SlopeCol = DataTable{:,SlopeColNum};
        SlopeCol(RowNum,1) = sfit.a;
        DataTable.(SlopeColName) = SlopeCol;
        
        ResultCol = DataTable{:,ResultColNum};
        ResultCol{RowNum,1} = SaveFileName; % this only saves the path to the result in the table entry
        DataTable.(ResultColName) = ResultCol;
        
        
        %---save variance ratio into datatable for stats:
        % determine variance across ucm and ort directions:
        Vucm = var(Result.ForceData.ForcesUCM);
        Vort = var(Result.ForceData.ForcesORT);
        VarRatio = Vucm/Vort;
        RatioColNum = find(strcmp(DataTable.Properties.VariableNames,['RatioVucmVort' ForceLevel]));
        
        VarRatioCol = DataTable{:,RatioColNum};
        VarRatioCol(RowNum,1) = VarRatio;
        DataTable.RatioVucmVort = VarRatioCol;
        
        save(DataTablePath, 'DataTable');
        
        
        %---- save another datatable for LME stats:
        DataTableLMEPath = [path2 'DataTableLME'];
        
        SubjectCol = [];
        for ii = 1:length(yplot)
            SubjectCol = [SubjectCol; SubjectCode];
        end
        
        SubjectTable = table(SubjectCol,cell2mat(CrVal.TestFUCM'), cell2mat(CrVal.TrueLabel'), yplot', yyn{TrialNum}');
        SubjectTable.Properties.VariableNames = {'Subject' 'UCMcontinuous' 'UCMlabel' 'ClassifierOutput' 'NormalizedOutput'};
        DataTableLME = [DataTableLME;SubjectTable];
        
        save(DataTableLMEPath, 'DataTableLME');

        
        %databaser stuff
        set(PublishInput.ProgramStatusIndicator,'string',['Step04_RunCSPLDA: Running... finished ' num2str(TrialNum) ' of ' num2str(NumTrials)],'foregroundcolor','g');drawnow;
    end
        
    %----outside of loop, after batch processing, plot overall ROC curves
    figure('color','w'); hold on;
    Sub = nominal(DataTableLME.Subject);
    AllSubjects= unique(Sub);
    for i = 1:length(AllSubjects)
        ind = find(Sub==AllSubjects(i));
        Input = nominal(DataTableLME.UCMlabel(ind));
        Output = DataTableLME.NormalizedOutput(ind);
        [Xs{i},Ys{i},Ts{i},AUCs{i}] = perfcurve(Input,Output,1); % find ROC curve - also computes area under curve (AUC)
        plot(Xs{i},Ys{i})
    end
    [Xall,Yall,Tall,AUCall] = perfcurve(nominal(DataTableLME.UCMlabel),DataTableLME.NormalizedOutput,1);
    hold on; plot(Xall,Yall,'linewidth',3,'color','k');
    line % plots line from (0,0) to (1,1)
    xlabel('False positive rate'); ylabel('True positive rate'); title(['ROC curves of normalized classifier output. Overall AUC = ' num2str(AUCall)])
    %----
    
    assignin('base', 'AUCs', cell2mat(AUCs));
    
    slopeucm
    assignin('base','yyn',yyn)
    assignin('base','DataTableLME',DataTableLME)
    
    %-----------------------------
catch err
    set(PublishInput.ProgramStatusIndicator,'string','Step04_RunCSPLDA had an error. See Matlab window for details.','foregroundcolor','r');drawnow;
    cd(BaseDirectory);
    rethrow(err);
end

PublishOutput = 1;

function y = test_bci(X,S,T,w,b) % from EEGLAB
% Prediction = test_bci(Raw-Block, Spatial-Flt, Temporal-Flt, Weights, Bias)
% A linear online detector for oscillatory processes.
%
% Input:
%   X : incoming raw sample data [Samples x Channels]
%   S,T: spatio-temporal filter [Channels x Filters], [Samples x Samples]
%   w,b: linear classifier [Filters x 1], [1 x 1]
%
% Output:
%   y : the prediction
%
% Notes:
%   y can be post-processed by sign(y) for pure classification, or 1./(1+exp(-y)) for logistic regression
%
% more info: http://sccn.ucsd.edu/wiki/Minimalist_BCI

% global B; % B is the buffer
B = [];
if any(size(B) ~= [length(T),length(S)])
    B = zeros(length(T),length(S));
end
B = [B;X]; B = B(end-length(T)+1:end,:);
y = log(var(T*(B*S)))*w - b;
