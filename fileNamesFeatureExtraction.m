function fileNamesFeatureExtraction( mainDirectory,fileName)
% This function is written instead of GUI_feature_matrices function
%  mainDirectory like this: 'E:\Thesis\Ebrahimi\Thesis Data\Data\
% Sound Description Toolbox\voice_for_featuring';
% mainDirectory is the address of voice files of one subject
% fileName is the name of excel file we want to put features in it (the
% name of person or any name for the features of one subject

mainFolderInfo = dir(mainDirectory);
sizeofMain = size(mainFolderInfo);

% Getting the names of folders in main directory
if ( sizeofMain(1) >= 3)
   folderNames = cell(sizeofMain(1)-2,1); %preallocating empty cell array 
    for i = 3 : sizeofMain(1)
       folderNames{i-2}= mainFolderInfo(i).name;
    end
end
sizeofFileNames = size(folderNames);
a = zeros(1,318);
featureMatrix = [];
featureMatrix = [featureMatrix; a];
for i = 1 : sizeofFileNames(1)
    %[fileName,diffAVG ] = ilpcFilterFeature( folderNames(i) );
    display(folderNames(i));
    Z = ComputeFeatureMatrix(mainDirectory, folderNames(i));
    %Z = ComputeFeatureMatrixSegmentedVoices(mainDirectory, folderNames(i)); % The version for segmenting voice before calculating
    featureMatrix = [featureMatrix; Z];
end

save_feature_matrix(mainDirectory,fileName,featureMatrix);


end

