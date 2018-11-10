for Z = [10, 15, 20, 25, 30]
    filechar = Z+"Bins.txt";
    fileID = fopen(filechar, 'w');
    outV = 0:0.1:Z;
    fprintf(fileID,'-INFINITY, ');
    for j = 1:length(outV)
        fprintf(fileID, '%f, ', outV(j));
    end
    fprintf(fileID, 'INFINITY');
    fclose(fileID);
    filechar2 = Z+"WEParams.txt";
    fileID = fopen(filechar2, 'w');
    fprintf(fileID,"10 32 10000 %u %u", 10*Z+2, 10*Z+2);
end