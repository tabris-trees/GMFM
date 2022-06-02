function B_data=Read_data(datafolder)
    data = importdata([datafolder,'log.txt']);
    datalist = data(4:end);
    datalength = length(datalist);
    subdatafile_list = erase(datalist,'.gz');
    subdataname = cell2mat(erase(datalist,'20200305dmin.min.gz'));
    latitudes = zeros(datalength,1);
    longitudes = zeros(datalength,1);
    B_xs = zeros(datalength,1);
    B_ys = zeros(datalength,1);
    B_zs = zeros(datalength,1);
    
    for i = 1:datalength
        subdatafile = subdatafile_list{i};
        subdata = fopen([datafolder subdatafile]);
        latitude_local = textscan(subdata,'%s %s %.3f %*[^\n]',1,'Headerlines',4);
        frewind(subdata);
        latitudes(i) = latitude_local{3};
        longitude_local = textscan(subdata,'%s %s %.3f %*[^\n]',1,'Headerlines',5);
        frewind(subdata);
        longitudes(i) = longitude_local{3};
        sub_B_data_all = textscan(subdata,'%s %s %s %.2f %.2f %.2f %*[^\n]','Headerlines',25);
        fclose(subdata);
        B_xs_local = sub_B_data_all{4};
        B_xs_local(B_xs_local==99999) = NaN;
        B_xs(i) = mean(B_xs_local,"omitnan");
        B_ys_local = sub_B_data_all{5};
        B_ys_local(B_ys_local==99999) = NaN;
        B_ys(i) = mean(B_ys_local,"omitnan");
        B_zs_local = sub_B_data_all{6};
        B_zs_local(B_zs_local==99999) = NaN;
        B_zs(i) = mean(B_zs_local,"omitnan");
    end
    
    B_data = [latitudes,longitudes,B_xs,B_ys,B_zs];
    [NAN_row_index, NAN_col_index] = find(isnan(B_data));
    NAN_row = unique(NAN_row_index);
    B_data(NAN_row,:) = [];
    subdataname(NAN_row,:) = [];
    writematrix(B_data,"MMD_20200305.txt",'Delimiter','tab');
    writematrix(subdataname,"STATIONS_NAME.txt",'Delimiter','tab');
end