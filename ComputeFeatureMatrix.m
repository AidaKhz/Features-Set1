function featureSelectedMatrix = ComputeFeatureMatrix( directory, filename)

% Initialize
startup;
format('long');
% warning('off'); % removed by Aida
% waveread changed to audio read in matlab 2015 and then raedwav in voicebox
% toolbox
% As file name is stored in cell type of data I used filename{1} istead of
% filename 

%% change by Aida
filename = strcat(directory, filename{1});
[A,B] = readwav(filename);% [y,Fs: samplerate,bits in sample] 

sampleNum = 882; % 20 ml second 
% frameNum = floor(size(A)/sampleNum);

%% changed by Aida


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % AE
    AudioPower = AP(A,size(A),B, sampleNum); % changed by Aida
    %meanAE = mean(AudioPower);% chanded by Aida
    covAE = cov(AudioPower);% chanded by Aida
    %meanDiffAE = mean(diff(AudioPower));% chanded by Aida
    diffAE = diff(AudioPower);
    covDiffAE = cov(diffAE); % chanded by Aida
    maxDiffAE = max(diffAE);
    minDiffAE = min(diffAE);
    clear AudioPower;
    clear diffAE;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Harmonic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % AFF
    standvar = h_mpeg7init(B,[],[],[],[]);
    f_num = size(A,1)/standvar.hopsize;       %hamid changes
    f0 = AFF(A,standvar,ceil(f_num)-1);%485);       %hamid changes
    %[f0,fc]=AFF(A,B);                      %hamid changes
    meanf0 = mean(f0);
    covf0 = cov(f0);
    %meanDifff0 = mean(diff(f0));
    diffF0 = diff(f0);
    covDifff0 = cov(diffF0);
    maxf0 = max(f0);
    minf0= min(f0);
    maxDiffF0 = max(diffF0);
    minDiffF0 = min(diffF0);
    clear f0;
    clear diffF0;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perceptual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TL_SONE
    p = struct('fs',B,'do_visu',0,'do_sone',1,'do_spread',1,'outerear','terhardt','fft_size',882,'hopsize',441,'bark_type','table','dB_max',96);
    [sone, Ntot, p] = TL_SONE(A(:,1),p);
    %figure; plot(sone');  %%Hamid
    %meanTL = mean(Ntot);
    covTL = cov(Ntot);
%     meanDiffTL = mean(diff(Ntot));
    diffTL = diff(Ntot);
    covDiffTL = cov(diff(Ntot));
    maxDiffTL = max(diffTL);
    minDiffTL = min(diffTL);
    
    meanSONE = mean(sone,2); meanSONE = meanSONE(1:8);
    for i=1:8 covSONE(i) = cov(sone(i,:)); end
%     for i=1:8 meanDiffSONE(i) = mean(diff(sone(i,:))); end
    for i=1:8 covDiffSONE(i) = cov(diff(sone(i,:))); end
    for i=1:8 maxDiffSONE(i) = max(diff(sone(i,:))); end
    for i=1:8 minDiffSONE(i) = min(diff(sone(i,:))); end

    clear Ntot;
    clear diffTL;
    clear sone;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ASC
    AudioSpectrumCentroid = ASC(filename,'PT10N1000F',0,[]);
    meanASC = mean(AudioSpectrumCentroid);
    covASC = cov(AudioSpectrumCentroid);
%     meanDiffASC = mean(diff(AudioSpectrumCentroid));
    diffASC = diff(AudioSpectrumCentroid);
    covDiffASC = cov(diff(AudioSpectrumCentroid));
    maxASC = max(AudioSpectrumCentroid);
    minASC = min(AudioSpectrumCentroid);
    maxDiffASC = max(diffASC);
    minDiffASC = min(diffASC);

    clear AudioSpectrumCentroid;
    clear diffASC;


    % ASR
    AudioSpectrumRolloff = ASR(filename,0.020);
    meanASR = mean(AudioSpectrumRolloff);
    covASR = cov(AudioSpectrumRolloff);
%     meanDiffASR = mean(diff(AudioSpectrumRolloff));
    diffASR = diff(AudioSpectrumRolloff);
    covDiffASR = cov(diff(AudioSpectrumRolloff));
    maxASR = max(AudioSpectrumRolloff);
    minASR = min(AudioSpectrumRolloff);
    maxDiffASR = max(diffASR);
    minDiffASR = min(diffASR);
    clear AudioSpectrumRolloff;
    clear diffASR;


    % ASS
    AudioSpectrumSpread = ASS(filename,'PT10N1000F',0,[]);
    meanASS = mean(AudioSpectrumSpread);
    covASS = cov(AudioSpectrumSpread);
%     meanDiffASS = mean(diff(AudioSpectrumSpread));
    diffASS = diff(AudioSpectrumSpread);
    covDiffASS = cov(diff(AudioSpectrumSpread));
    minASS = min(AudioSpectrumSpread);
    maxASS = max(AudioSpectrumSpread);
    maxDiffASS = max(diffASS);
    minDiffASS = min(diffASS);

    clear AudioSpectrumSpread;
    clear diffASS;


    % MFCC
    [ceps,freqresp,fb,fbrecon,freqrecon] = MFCC(A, B);
    [C1 C2] = size(ceps);
    for i=1:C1 for j=1:C2 if(isinf(ceps(i,j)) || isnan(ceps(i,j))) ceps(i,j) = 0; end; end; end;
    meanMFCC = mean(ceps(2:end,:),2); %removing first coefficient of MFCC
    for i=2:C1 covMFCC(i-1) = cov(ceps(i,:)); end
%     for i=1:C1 meanDiffMFCC(i) = mean(diff(ceps(i,:))); end
    for i=2:C1 covDiffMFCC(i-1) = cov(diff(ceps(i,:))); end
    for i=2:C1 maxMFCC(i-1) = max(ceps(i,:)); end
    for i=2:C1 minMFCC(i-1) = min(ceps(i,:)); end
    for i=2:C1 maxDiffMFCC(i-1) = max(diff(ceps(i,:))); end
    for i=2:C1 minDiffMFCC(i-1) = min(diff(ceps(i,:))); end

    clear ceps;
    clear freqresp;
    clear fb;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temporal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % AC
    ACcoefs = AC(A);

    % LAT
    energy_bp = h_energy(A, B, 2000, 2);
    LogAttackTime = LAT(energy_bp, 20); % changed from 50 to 20 by Aida;paper said

    % TC
    TemporalCentroid = TC(energy_bp);

    clear energy_bp;

    % ZCR
    [ZCR1 avZCR] = ZCR(A(:,1),0.020); %changed by Aida from 16ms to 20
    meanZCR = mean(ZCR1);
    covZCR = cov(ZCR1);
    %meanDiffZCR = mean(diff(ZCR1));
    %covDiffZCR = cov(diff(ZCR1));

    clear ZCR1;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Various %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ASF
    [AudioSpectrumFlatness ,lo_edge, hi_edge, XMLFile ] = ASF(A,B,'PT10N1000F',250,500,0,[]); % for 20 ms
    [ASF1 ASF2] = size(AudioSpectrumFlatness);
    meanASF = mean(AudioSpectrumFlatness,1);
    for i=1:ASF2 
        covASF(i) = cov(AudioSpectrumFlatness(i,:)); 
    end
    for i=1:ASF2 
        maxASF(i) = max(AudioSpectrumFlatness(i,:)); 
        minASF(i) = min(AudioSpectrumFlatness(i,:));
    end

    clear AudioSpectrumFlatness;
    
    
    % AudioPower
    AudioAmplitude = AP2(A,size(A),B,sampleNum); % chanded by Aida
    %meanAP = mean(AudioAmplitude);% chanded by Aida
    covAP = cov(AudioAmplitude);% chanded by Aida
    diffAP = diff(AudioAmplitude);
    covDiffAP = cov(diffAP); % chanded by Aida
    maxDiffAP = max(diffAP);
    minDiffAP = min(diffAP);
    
    clear AudioAmplitude;
    clear diffAP;


    %iLPC  By Aida
    lpcord = 2+floor(B/1000);
    %sp = filter([1 -exp(-2*pi*50/B)],1,A);                    % preemphasis zero is at 50 Hz
    sp = A;
    [lpar,lpe,lpk] = lpcauto(sp,lpcord,floor([0.01 0.02]*B));  % 20ms frame with 10ms frame increment(overlap is 10)
    % residual energy
    overlap = lpk(1,2)-lpk(2,1)+1;                               % overlap between adjacent frames
    iLPC = lpcifilt(sp,lpar,lpk(:,1)+floor(overlap/2),0,overlap/4);
    iLPC = iLPC .* 10^3; 
    meanILPC = mean(iLPC);
    covILPC = cov(iLPC);
    diffILPC = diff(iLPC);
    covDiffILPC = cov(diffILPC);
    maxDiffILPC = max(diffILPC);
    minDiffILPC = min(diffILPC);
    meanAbsILPC = mean(abs(iLPC));
    covAbsILPC = cov(abs(iLPC));
    diffAbsILPC = diff(abs(iLPC));
    covDiffAbsILPC = cov(diffAbsILPC);
    maxDiffAbsILPC = max(diffAbsILPC);
    minDiffAbsILPC = min(diffAbsILPC);
    
    kurILPC = kurtosis(iLPC);
    skewILPC = skewness(iLPC);
    
%     fftILPC = signalFFT(A,size(A),sampleNum,20);
%     [FFT1 FFT2] = size(fftILPC);
%     meanFFTILPC = mean(fftILPC,2);
%     for i=1:FFT2 
%         covFFTILPC(i) = cov(fftILPC(i,:)); 
%     end
    fftILPC = abs(fft(iLPC,20));
    
    clear iLPC;
    clear diffILPC;
    clear diffAbsILPC;
    
    fftSignal = abs(fft(A,20));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Put all features together %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[4, 6, 9, 11, 48, 49, 51, 56, 80, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                            %  162, 163, 164]
    featureMatrix = [];
    featureMatrix = [featureMatrix covAE covDiffAE maxDiffAE minDiffAE];%1-4
    featureMatrix = [featureMatrix meanf0 covf0 covDifff0 maxf0 minf0 maxDiffF0 minDiffF0];%5-11
    featureMatrix = [featureMatrix covTL covDiffTL maxDiffTL minDiffTL];%12-15
    featureMatrix = [featureMatrix covSONE covDiffSONE maxDiffSONE minDiffSONE];%16-47
    featureMatrix = [featureMatrix meanASC covASC covDiffASC maxASC minASC maxDiffASC minDiffASC];%48-54
    featureMatrix = [featureMatrix meanASR covASR covDiffASR maxASR minASR maxDiffASR minDiffASR];%55-61
    featureMatrix = [featureMatrix meanASS covASS covDiffASS maxASS minASS maxDiffASS minDiffASS];%62-68
    featureMatrix = [featureMatrix meanMFCC' covMFCC covDiffMFCC maxMFCC minMFCC maxDiffMFCC minDiffMFCC];%69-92-115-138-161-184-207-229
    featureMatrix = [featureMatrix ACcoefs];%'];%230-242
    featureMatrix = [featureMatrix LogAttackTime];%243
    featureMatrix = [featureMatrix TemporalCentroid];%244
    featureMatrix = [featureMatrix meanZCR covZCR];%245-246
    featureMatrix = [featureMatrix meanASF covASF maxASF minASF];%247-262
    featureMatrix = [featureMatrix covAP covDiffAP maxDiffAP minDiffAP]; %263-266 Added By Aida
    featureMatrix = [featureMatrix meanILPC covILPC covDiffILPC maxDiffILPC minDiffILPC];%267-271
    featureMatrix = [featureMatrix meanAbsILPC covAbsILPC covDiffAbsILPC maxDiffAbsILPC minDiffAbsILPC];%272-276
    featureMatrix = [featureMatrix kurILPC skewILPC]; %277-278
    featureMatrix = [featureMatrix fftILPC']; %279-298
    featureMatrix = [featureMatrix fftSignal']; %299-318    
    %     display(meanMFCC') 
 
    

    %featureMatrixes = [featureMatrixes ; featureMatrix];


    %{
    featureSelectedMatrix = [];
    featureSelectedMatrix = [featureSelectedMatrix meanAP covAP meanDiffAP covDifff0 meanTL covTL meanSONE(4) meanSONE(8) covDiffASR meanASS covASS];
    %}


    clear meanAP covAP meanDiffAP covDiffAP;
    clear meanf0 covf0 meanDifff0 covDifff0;
    clear meanTL covTL meanDiffTL covDiffTL;
    clear meanSONE covSONE meanDiffSONE covDiffSONE;
    clear meanASC covASC meanDiffASC covDiffASC;
    clear meanASR covASR meanDiffASR covDiffASR;
    clear meanASS covASS meanDiffASS covDiffASS;
    clear meanMFCC covMFCC meanDiffMFCC covDiffMFCC;
    clear ACcoefs LogAttackTime TemporalCentroid;
    clear meanZCR covZCR meanDiffZCR covDiffZCR;
    clear meanASF covASF meanDiffASF covDiffASF;
    clear meanAE covAE meanDiffAE covDiffAE;
    
     %featureSelectedMatrix = featureMatrixes;
    featureSelectedMatrix = featureMatrix;

end


