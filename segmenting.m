function number_of_segmented_files = segmenting(filename) 

%% File: segmenting.m
%%
%% ------------- segmentation of audio files to 20ms parts--------------------
%%
%% The function of this subroutine is to segment an audio file to smaller parts of 20 ms length 
%% and save them in a subfolder with the name of segmented
%% 
%% Author: Aida


%% Example: number_of_segmented_files = segmenting(filename);
%%          plot(AudioPower_SeriesOfScalar);
                

sampleNum = 882;
frameNum = floor(totalSampleNum/sampleNum);

j = 1;
for i = j:frameNum
   signal = auData(1+(i-1)*sampleNum:i*sampleNum);
   writewav(signal,fs,strcat(filename,'//segmented//',num2str(i)));
end

AudioPower_SeriesOfScalar = meanValues;