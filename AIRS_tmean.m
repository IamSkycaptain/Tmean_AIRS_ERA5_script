clc;clear;
first = 0;
workdir = 'E:\EarthData\hdf_standard\2016';
%workdir = 'D:\EarthData\hdf_data\results';
fileID = fopen('E:\EarthData\hdf_standard\result_tmean_pwv\QC_rata\hdf_2016.csv','w');

fprintf(fileID,'%s\n','ymd,doy,lat,lon,h_msl,ts,press,pwv,tm,temp_qc,pwv_qc,press_qc,h_msl_qc,layers,subdoy');
if ~isfolder(workdir)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', workdir);
    uiwait(warndlg(errorMessage));
    return;
end
x = [ 98.29840079,  98.57873878,  98.5635491 ,  98.42412829,98.17172572,  98.14499225,  98.52392419,  99.06538105,99.06294513,  99.39090958,  99.56159711,  98.76389732,98.39998909,  98.25308211,  98.23784721,  98.37052277,98.62441157,  98.86421513,  99.02571391,  99.19482172,99.30847693,  99.47268,  99.59452565,  99.71058974,99.6978771 ,  99.66583124,  99.69455757,  99.77428139,99.91478921, 100.0830995 , 100.20729924, 100.26796931,100.47259482, 100.59326747, 100.76009198, 100.98319707, 101.01230167, 100.96309425, 101.00182246, 101.17336224,101.31446721, 101.54976712, 101.85710781, 102.11453313,101.6849696 , 101.31780865, 100.92781209, 100.50606142,100.29940148, 100.22004819,  99.95726154,  99.86745716,99.71775043,  99.45029485,  99.303078  ,  99.21908727, 99.15997116,  99.22143178,  99.54157347,  99.97702896,99.96932246,  99.98206359, 100.08891373,  99.98820546,100.08411769, 100.26042824, 100.35931946, 100.47157177,100.97462764, 100.90148967, 101.73680297, 102.20030637,102.43955908, 102.69471503, 102.64997608, 102.38442819,102.79857206, 103.23906632, 104.00705273, 104.76927808,105.13087869, 105.43767389, 105.53696459, 105.63760105,105.2484791 , 104.81086913, 104.83406333, 104.0804889 ,103.3129423 , 102.63904534, 102.12125634, 101.622087  ,101.03897126, 101.31041367, 101.16255652, 100.52345339,100.50538974,  99.98210636,  99.68284447,  99.51912092,99.44482043,  99.50521803,  99.19328215,  99.03889293,99.0033505 ,  98.98985549,  98.9019462 ,  98.54298812,98.31621951,  98.02175029,  97.90005873,  97.62770008,97.51884469,  97.60104099,  97.71569472,  98.29840079];
y = [17.0479819 , 16.47975015, 15.85432878, 15.4783534 , 15.14368435,       15.02610422, 14.31016641, 13.69641793, 13.12674405, 12.33616832,       11.7163866 , 10.67790099,  9.30285363,  8.59961004,  8.36024902,        8.21053065,  8.36594527,  8.04775306,  7.81802411,  7.67214633,        7.48444554,  7.31148257,  7.18417899,  7.15332376,  7.04751791,        7.01498262,  6.92804267,  6.83251573,  6.68253089,  6.53731032,        6.6387736 ,  6.58886349,  6.47151007,  6.41300909,  6.35033095,        6.17346263,  6.0131534 ,  5.85635031,  5.66798219,  5.56409881,        5.75396399,  5.842421  ,  5.70581662,  6.15122819,  6.65698211,        6.8870937 ,  6.88355   ,  7.25871958,  8.12291868,  8.4311311 ,        8.62773338,  9.21294207,  9.33073424,  9.19258225,  9.26440829,        9.5808528 ,  9.82355776, 10.43170182, 11.29393305, 12.20246021,       12.67626123, 12.79889208, 13.07993967, 13.27898697, 13.44737365,       13.46889118, 13.48651256, 13.50062512, 13.46708291, 12.64643829,       12.71611997, 12.44843246, 12.2108935 , 12.15973996, 12.62107864,       13.35396482, 13.89688238, 14.34625665, 14.36153547, 14.36721265,       14.22096376, 14.455419  , 14.99737379, 15.66959627, 16.06829816,       16.50681069, 17.46520073, 18.25282297, 18.43297247, 17.92174367,       18.25063622, 17.92231906, 17.54439034, 19.57010206, 19.73047998,       19.62915545, 20.31919451, 20.52304334, 20.40684871, 20.44639241,      20.33363368, 20.154256  , 20.14878533, 20.07740803, 19.95269041,19.85614311, 19.81340609, 19.72482935, 19.72891734, 19.84248363,19.68213979, 18.97673738, 18.59399836, 18.11871785, 17.59760336,17.0479819 ];
thai_boundary = polyshape({x},{y});
%plot(thai_boundary)
%hold on

filelist = dir(fullfile(workdir, '**\*.hdf'));
filelist = filelist(~[filelist.isdir]);
%filePattern = fullfile(workdir, '*.hdf');
%matFiles = dir(filePattern);

%numfile = length(matFiles);
%conn = database('hdf','',''); % ODBC
vgroup = '/L2_Standard_atmospheric&surface_product/Data Fields/';
colnames = {'ymd','doy','lat','lon','TSurfAir','PSurfStd',...
            'totH2OStd','GP_Surface','Tm'};     
for k=1:length(filelist)
    %fn = matFiles(k).name;
    %fullname = fullfile(workdir, fn);
    fullname = fullfile(filelist(k).folder, filelist(k).name);
    fprintf(1, 'Now reading %s\n', fullname);
    %info = hdfinfo(fullname);
    %TAI Time: floating-point elapsed seconds since Jan 1, 1993
    Time(:,:) = hdfread(fullname, '/L2_Standard_atmospheric&surface_product/Geolocation Fields/Time', 'Index', {[1  1],[1  1],[45  30]});
    pieces = regexp(fullname, '\.', 'split');
    yyyy=eval(pieces{2});  %mm=eval(pieces{3}); dd=eval(pieces{4}); granule = eval(pieces{5}); 
    % Start times of granules 00:05:25Z see Manual Appendix A. Level 2 Product Interface Specifications
    %hr = 5.0/60.0+25./(60.*60.)+0.1*(granule-1.);%เวลาเปิดถ่ายภาพ ใชัเวลา 6 นาทีต่อภาพ
 
    TAirStd(:,:,:) = hdfread(fullname, strcat(vgroup,'TAirStd'), 'Index', {[1  1  1],[1  1  1],[45  30  28]});
    TAirStd_QC(:,:,:)=hdfread(fullname, strcat(vgroup,'TAirStd_QC'), 'Index', {[1  1  1],[1  1  1],[45  30  28]});

    H2OMMRLevStd(:,:,:) = hdfread(fullname, strcat(vgroup,'H2OMMRLevStd'), 'Index', {[1  1  1],[1  1  1],[45  30  15]});
    H2OMMRLevStd_QC(:,:,:)=hdfread(fullname, strcat(vgroup,'H2OMMRLevStd_QC'), 'Index', {[1  1  1],[1  1  1],[45  30 15]});
   
    RelHum(:,:,:) = hdfread(fullname, strcat(vgroup,'RelHum'), 'Index', {[1  1  1],[1  1  1],[45  30  15]});
    RelHum_QC(:,:,:)=hdfread(fullname, strcat(vgroup,'RelHum_QC'), 'Index', {[1  1  1],[1  1  1],[45  30 15]});
   
    TSurfAir(:,:) = hdfread(fullname, strcat(vgroup,'TSurfAir'), 'Index', {[1  1],[1  1],[45  30]});
    TSurfAir_QC(:,:) = hdfread(fullname, strcat(vgroup,'TSurfAir_QC'), 'Index', {[1  1],[1  1],[45  30]});

    PSurfStd(:,:) = hdfread(fullname, strcat(vgroup,'PSurfStd'), 'Index', {[1  1],[1  1],[45  30]});
    PSurfStd_QC(:,:) = hdfread(fullname, strcat(vgroup,'PSurfStd_QC'), 'Index', {[1  1],[1  1],[45  30]});

    totH2OStd(:,:) = hdfread(fullname, strcat(vgroup,'totH2OStd'), 'Index', {[1  1],[1  1],[45  30]});
    totH2OStd_QC(:,:)=hdfread(fullname,  strcat(vgroup,'totH2OStd_QC'), 'Index', {[1  1],[1  1],[45  30]});
    %Geopotential Height of surface (m above mean sea level)
    GP_Surface(:,:) = hdfread(fullname, strcat(vgroup,'GP_Surface'), 'Index', {[1  1],[1  1],[45  30]}); %
    GP_Surface_QC(:,:)=hdfread(fullname,  strcat(vgroup,'GP_Surface_QC'), 'Index', {[1  1],[1  1],[45  30]});
    
    pressH2O_dummy = hdfread(fullname, strcat(vgroup,'pressH2O'), 'Fields', 'pressH2O', 'FirstRecord',1 ,'NumRecords',15);  
    pressH2O=pressH2O_dummy{1,1};

    Latitude(:,:) = hdfread(fullname, '/L2_Standard_atmospheric&surface_product/Geolocation Fields/Latitude', 'Index', {[1  1],[1  1],[45  30]});
    Longitude(:,:) = hdfread(fullname, '/L2_Standard_atmospheric&surface_product/Geolocation Fields/Longitude', 'Index', {[1  1],[1  1],[45  30]});      

    %Estimate Tm
    %Reference : M.Wallace, J., & Hobbs, P. V. (2006). Atmospheric Science
    %2nd Edition (J. Helé Ed.). pp-80
    epsilon=0.62198; % epsilon=Mv/Md -> Mv = 18.01 gmol-1
    qc = 2;
    for geotrack=1:45
        for xtrack=1:30     
            lat = Latitude(geotrack,xtrack);
            lon = Longitude(geotrack,xtrack);
            if isinterior(thai_boundary,lon,lat)
                first = 1;
                ts_qc = TSurfAir_QC(geotrack,xtrack);
                pwv_qc = totH2OStd_QC(geotrack,xtrack);
                press_qc= PSurfStd_QC(geotrack,xtrack);
                if ts_qc < qc && pwv_qc < qc
                    %plot(lon,lat,'r.')
                    dt = datetime(1993,1,1)+seconds(Time(geotrack,xtrack)-10.0); %10 leab second after 1993 : UTC = TAI-leapsecond
                    ymd = datestr(dt,'yyyy-mm-dd_HH:MM:SS');
                    doy = decyear(year(dt),month(dt),day(dt),hour(dt),minute(dt),second(dt)); %UTC time
                    subdoy = doy-year(dt);
                    ts = TSurfAir(geotrack,xtrack);
                    press= PSurfStd(geotrack,xtrack);
                    pwv = totH2OStd(geotrack,xtrack);
                    h_msl_qc = GP_Surface_QC(geotrack,xtrack);
                    h_msl = GP_Surface(geotrack,xtrack);       
                    Tm_numerator=0;
                    Tm_denumerator=0;
                    chk = 0;
                    % New Method
                    % Ref Tm model : Improving the weighted mean temperature model: A case study using nine year (2007–2015) radiosonde and COSMIC data in Hong Kong
                    %for layer=1:length(pressH2O)-1
                    %    w_qc = H2OMMRLevStd_QC(geotrack,xtrack,layer);
                    %    t_qc = TAirStd_QC(geotrack,xtrack,layer);
                        
                    %    w_qc1 = H2OMMRLevStd_QC(geotrack,xtrack,layer+1);
                    %    t_qc1 = TAirStd_QC(geotrack,xtrack,layer+1);      
                    %    %fprintf('%.0f, %.0f, %.0f, %.0f\n',w_qc,t_qc,w_qc1,t_qc1);
                    %    if w_qc < qc && t_qc < qc && w_qc1 < qc && t_qc1 < qc
                    %        chk = chk + 1;
                    %        w0=H2OMMRLevStd(geotrack,xtrack,layer)*1e-3;
                    %        p0=pressH2O(layer);
                    %        e0=(w0/(w0+epsilon))*p0;
                    %        t0=TAirStd(geotrack,xtrack,layer);
                            
                    %        w1=H2OMMRLevStd(geotrack,xtrack,layer+1)*1e-3;
                    %        p1=pressH2O(layer+1);
                    %        e1=(w1/(w1+epsilon))*p1;
                    %        t1=TAirStd(geotrack,xtrack,layer+1);                            
                            
                    %        et = e0/t0 + e1/t1;
                    %        et2= e0/t0^2 + e1/t1^2;
                            
                    %        Tm_numerator=Tm_numerator+et;
                    %        Tm_denumerator=Tm_denumerator+et2;    
                    %        %disp([p,w,Tair])
                    %    end                   
                    %end                    
                    % Normal method
                    for layer=1:length(pressH2O)
                        w_qc = H2OMMRLevStd_QC(geotrack,xtrack,layer);
                        t_qc = TAirStd_QC(geotrack,xtrack,layer);
                        if w_qc < qc && t_qc < qc
                            chk = chk + 1;
                            %method 1
                            w=H2OMMRLevStd(geotrack,xtrack,layer)*1e-3;
                            p_partwv=(w/(w+epsilon))*pressH2O(layer);
                            Tair=TAirStd(geotrack,xtrack,layer);
                            Tm_numerator=Tm_numerator+p_partwv/Tair;
                            Tm_denumerator=Tm_denumerator+p_partwv/(Tair^2);    
                            %disp([p,w,Tair])
                       end
                       
                            %method2 is not working
                            %relhum  = RelHum(geotrack,xtrack,layer);
                            %p_partwv=(pressH2O(layer+1)*relhum)/(relhum-(relhum-1)*epsilon);  
                            %fprintf('h20=%.3f, hum=%.3f\n',p_partwv,p_partwv);                     
                    end
                    
                    
                    
                    % Rata
                    %for layer=1:length(pressH2O)-1
                    %       %W=H2OMMRStd(k,geotrack,xtrack,layer)/100;
                    %       %p_partwv=W*((pressH2O(layer)+pressH2O(layer+1))/2)/(W+0.62198);
                    %    w_qc = H2OMMRLevStd_QC(geotrack,xtrack,layer);
                    %    t_qc = TAirStd_QC(geotrack,xtrack,layer+1);
                    %    if w_qc < qc && t_qc < qc
                    %        w=H2OMMRLevStd(geotrack,xtrack,layer)*1e-3;
                    %        p_partwv=(w/(w+epsilon))*pressH2O(layer+1);
                    %        Tair=TAirStd(geotrack,xtrack,layer+1);
                    %        Tm_numerator=Tm_numerator+p_partwv/Tair;
                    %        Tm_denumerator=Tm_denumerator+p_partwv/(Tair^2);  
                    %        %fprintf('%.3f %.9f,%.0f\n',pressH2O(layer),w,w_qc);
                    %    end    
                    %end
                    %disp(chk)
                    tm =Tm_numerator/Tm_denumerator;
                    %disp([ymd,' L: ', num2str(chk),'   tm:  ',num2str(tm),' ts: ', num2str(ts)])
                    %if tm>273
                    fprintf(fileID,'%s,%.10f,%.8f,%.8f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d,%d,%d\n',ymd,doy,lat,lon,h_msl,ts,press,pwv,tm,ts_qc,pwv_qc,press_qc,h_msl_qc,chk,subdoy);
                    %end
                end
            end
        end
    end
    %if k==10
    %    break
    %end
    %break
end
%axis equal
fclose(fileID);
%close(curs);
%close(conn);
  