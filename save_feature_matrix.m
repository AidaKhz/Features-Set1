% Sinartisi save_feature_matrix: Kalei tin sinartisi feature_matrix gia
% na swsei ton pinaka features enos arxeiou ixou ston disko se ena arxeio
% *.sfm.

function [] = save_feature_matrix(mainDirectory,filename,featureMatrix)

sheet = 1;
n=1;
fileID = fopen(strcat(mainDirectory,filename,'.txt'),'w');
for i = 1 : size(featureMatrix,1) % for the number of name of files
    for j = 1 : size(featureMatrix,2)
%         xlRange = strcat('A', num2str(n)); % Indicates which row in excel file
%         xlswrite(strcat(mainDirectory,filename,'.xlsx'),featureMatrix(i,:),sheet,xlRange);
        if j == size(featureMatrix,2)
            fprintf(fileID,'%8.16f ',featureMatrix(i,j));
        else
            fprintf(fileID,'%8.16f , ',featureMatrix(i,j));
        end
    end
    fprintf(fileID,' \n');
end
 fclose(fileID);
end
% fid = fopen([filename '.csv'],'a+t');
% %fid = fopen('Autistics.txt','a+t');
% %Z = ComputeFeatureMatrix(filename);
% %fprintf(fid, '%4.8f\n',Z);
% 
% for i = 1 : size(featureMatrix,1)
%     fprintf(fid, '%4.8f,',featureMatrix(i,:));
%     fprintf(fid, 'A');
% end
%fprintf(fid, 'A');
%fclose(fid);
%clear Z;