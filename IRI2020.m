function Out = IRI2020(lat,lon,time,alt,varargin)
% MATLAB wrapper for the IRI model and MEX function
%
% INPUTS
%   latitude(deg), longitude(deg), time(datenum or datetime), altitude(km)
%
% OUTPUTS
%   structure with fields
%   latitude(deg), longitude(deg), time(datetime), altitude(km),
%   electron density (Ne/m^3) and various profile parameters
%
% OPTIONAL INPUTS
%
%   'JF', JF array: use the keyword 'JF' followed by a length 50 numerical
%   vector of option switches as defined in irisub.for. The wrapper
%   interprets JF<=0 as false, JF>0 as true. By default, the wrapper uses a
%   resonable set of options that should work for most purposes. All
%   options related to non-Ne outputs or user defined input parameters will
%   be ignored.
%
%   By default, IRI2016 works in profile mode. Lat and Lon must be the
%   same size, and the output will have Ne at every altitude for each
%   Lat/Lon pair. Time can either have size 1, and all profiles will be run
%   for the same time, or else it can be the same size as lat/lon in which
%   case each corresponding [lat,lon,time] will be run
%
%       There is no limit on how many altitudes can be given for a profile,
%       unlike the base IRI. This function will break calls to IRI_sub up
%       if more than 1000 altitudes are requested
%
%   'sat' mode: in this mode, all inputs must be the same size, including
%   altitude
%
%   'map' mode: in this mode, all unique values of inputs are expanded in a
%   4d grid. The output Ne is of a size [n time,n lat,n lon,n alt]
%
%   Ben Reid 2020

ErrorStr = 'MATLAB:IRI2020';

% PATH TO MEX FILE
[RootDir, ~, ~] = fileparts(mfilename('fullpath'));

% PATH TO NRLMSISE2.1 FORTRAN CODE
ModelDir = fullfile(RootDir,'IRI-2020');

ModelDir = fullfile(ModelDir,filesep);

if isempty(ModelDir)
    error(ErrorStr,[...
        'IRI2016 model directory not configured! Please edit',...
        ' IRI2016.m to specify the appropriate path'])
elseif ~exist(ModelDir,'dir')
    error(ErrorStr,[...
        'Invalid model directory. Directory does not exist. Please edit',...
        ' IRI2016.m to specify the appropriate path'])
end
JF = ones(50,1);

JF([4,5,6,22,23,30,33,34,35,39,40,47])=-1;

%    JF switches to turn off/on (.true./.false.) several options
%
%    i       .true.                  .false.          standard version
%    -----------------------------------------------------------------
%    1    Ne computed            Ne not computed                     t
%    2    Te, Ti computed        Te, Ti not computed                 t
%    3    Ne & Ni computed       Ni not computed                     t
%    4    B0,B1 - Bil-2000       B0,B1 - other models jf(31)     false
%    5    foF2 - CCIR            foF2 - URSI                     false
%    6    Ni - DS-1995 & DY-1985 Ni - RBV-2010 & TBT-2015        false
%    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t
%    8    foF2 from model        foF2 or NmF2 - user input           t
%    9    hmF2 from model        hmF2 or M3000F2 - user input        t
%   10    Te - Standard          Te - Using Te/Ne correlation        t
%   11    Ne - Standard Profile  Ne - Lay-function formalism         t
%   12    Messages to unit 6     to messages.txt on unit 11          t
%   13    foF1 from model        foF1 or NmF1 - user input           t
%   14    hmF1 from model        hmF1 - user input (only Lay version)t
%   15    foE  from model        foE or NmE - user input             t
%   16    hmE  from model        hmE - user input                    t
%   17    Rz12 from file         Rz12 - user input                   t
%   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
%   19    F1 probability model   only if foF1>0 and not NIGHT        t
%   20    standard F1            standard F1 plus L condition        t
% (19,20) = (t,t) f1-prob, (t,f) f1-prob-L, (f,t) old F1, (f,f) no F1
%   21    ion drift computed     ion drift not computed              t
%   22    ion densities in %     ion densities in m-3                t
%   23    Te_tops (Bil-1985)     Te_topside (TBT-2012)           false
%   24    D-region: IRI-1990     FT-2001 and DRS-1995                t
%   25    F107D from APF107.DAT  F107D user input (oarr(41))         t
%   26    foF2 storm model       no storm updating                   t
%   27    IG12 from file         IG12 - user                         t
%   28    spread-F probability 	 not computed                        t
%   29    IRI01-topside          new options as def. by JF(30)       t
%   30    IRI01-topside corr.    NeQuick topside model   	     false
% (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick, (t,f) IRIcor2
%   31    B0,B1 ABT-2009	     B0 Gulyaeva-1987 h0.5               t
% (4,31) = (t,t) Bil-00, (f,t) ABT-09, (f,f) Gul-87, (t,f) not used
%   32    F10.7_81 from file     F10.7_81 - user input (oarr(46))    t
%   33    Auroral boundary model on/off  true/false	             false
%   34    Messages on            Messages off                        t
%   35    foE storm model        no foE storm updating           false
%   36    hmF2 w/out foF2_storm  with foF2-storm                     t
%   37    topside w/out foF2-storm  with foF2-storm                  t
%   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
%   39    hmF2 (M3000F2)         new models                      false
%   40    hmF2 AMTB-model        Shubin-COSMIC model                 t
% (39,40) = (t,t) hmF2-old, (f,t) AMTB, (f,f) Shubin, (t,f) not used
%   41    Use COV=F10.7_365      COV=f(IG12) (IRI before Oct 2015)   t
%   42    Te with PF10.7 dep.	 w/o PF10.7 dependance               t
%   43    B0 from model          B0 user input in OARR(10)           t
%   44    B1 from model          B1 user input in OARR(35)           t
%   45    HNEA=65/80km dya/night HNEA user input in OARR(89)         t
%   46    HNEE=2000km 	         HNEE user input in OARR(90)         t
%   47    CGM computation on 	 CGM computation off             false
%   48    Ti  Tru-2021           Bil-1981                            t
%      ....
%   50
%   ------------------------------------------------------------------

mode='default';
Update=false;
Storm=false;
F107=[];
F107_81=[];
i=1;
while i <= numel(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        if strcmpi(varargin{i},'default')
            mode = default;
        elseif strcmpi(varargin{i},'sat')
            mode = 'sat';
        elseif strcmpi(varargin{i},'map')
            mode='map';
        elseif strcmpi(varargin{i},'update')
            Update=true;
        elseif strcmpi(varargin{i},'jf')
            JF = varargin{i+1};
            i=i+1;
        elseif strcmpi(varargin{i},'F107')
            F107 = varargin{i+1};
            i=i+1;
        elseif strcmpi(varargin{i},'F107_81')
            F107_81 = varargin{i+1};
            i=i+1;
        else
            error(ErrorStr,[...
                'Invalid option: ',varargin{i}]);
        end
    else
        error(ErrorStr,[...
            'Invalid option.'])
    end
    i=i+1;
end

if Update % update indices.dat file in iri directory
    read_apf(fullfile(ModelDir,'apf107.dat'),true);
    read_ig_rz(fullfile(ModelDir,'ig_rz.dat'),true);
end

if any((lat(:)<-90) | (lat(:)>90))
    error(ErrorStr,...
        'Input latitude out of range [-90,90]')
end

if any((lon(:)<-180) | (lon(:)>360))
    error(ErrorStr,...
        'Input longitude out of range [-180,360]')
end

try

    OldDir = pwd;
    cd(ModelDir);

    if ~exist('mex_iri2020','file')
        warning(ErrorStr,...
            'Required MEX file mex_iri2020 does not exist. Attempting to compile.')

        mex('-compatibleArrayDims','-DO3','-Dstatic','-DfPIC', ...
            '-Dffast-math','-Dmarch=native', ...
            '-output',fullfile(RootDir,'mex_iri2020'), ...
            fullfile(ModelDir,'irisub.F'), ...
            fullfile(ModelDir,'irifun.F'), ...
            fullfile(ModelDir,'iritec.F'),...
            fullfile(ModelDir,'iridreg.F'), ...
            fullfile(ModelDir,'iriflip.F'), ...
            fullfile(ModelDir,'cira.F'), ...
            fullfile(ModelDir,'igrf.F'),...
            fullfile(ModelDir,'rocdrift.F'),...
            fullfile(RootDir,'mex_iri.F'))
    end


    % check to make sure the IRI will have the needed files
    if ~all(arrayfun(@(n) exist(sprintf('ccir%02i.asc',n),'file'),11:22))
        error(ErrorStr,[...
            'CCIR files required by IRI are missing from the model folder',...
            ' Check to make sure IRI2020.m is configured correctly.'])
    end
    if ~all(arrayfun(@(n) exist(sprintf('mcsat%02i.dat',n),'file'),11:22))
        error(ErrorStr,[...
            'mcsat files required by IRI are missing from the model folder',...
            ' Check to make sure IRI2020.m is configured correctly.'])
    end
    if ~all(arrayfun(@(n) exist(sprintf('ursi%02i.asc',n),'file'),11:22))
        error(ErrorStr,[...
            'URSI files required by IRI are missing from the model folder',...
            ' Check to make sure IRI2020.m is configured correctly.'])
    end
    if ~exist('ig_rz.dat','file')
        error(ErrorStr,[...
            'File ig_rz.dat required by IRI is missing from the model folder',...
            ' Check to make sure IRI2020.m is configured correctly, or', ...
            ' use the ''update'' keyword to download.'])
    end
    if ~exist('apf107.dat','file')
        error(ErrorStr,[...
            'File apf107.dat required by IRI is missing from the model folder',...
            ' Check to make sure IRI2020.m is configured correctly, or', ...
            ' use the ''update'' keyword to download.'])
    end


    time = datetime(datenum(time),'ConvertFrom','datenum');

    if strcmpi(mode,'default')
        tvec=[numel(lat),numel(lon)];
        tvec=tvec==tvec';
        if ~all(tvec(:))
            error(ErrorStr,[...
                'In "default" mode lat and lon must have the same size'])
        end

        if numel(time) == 1
            time = repmat(time, size(lon));
        elseif numel(lon) == 1
            lon=repmat(lon,size(time));
            lat=repmat(lat,size(time));
        else
            error(ErrorStr,[...
                'In "default" mode, the number of input times must be either',...
                ' equal to one, or the number of lat/lon pairs'])
        end

        Out.dates = time(:)';
        Out.lon = lon(:);
        Out.lat = lat(:);
        Out.alt = alt(:);

        vbeg=[];
        vend=[];
        vstp=[];

        split_alt;

        lon = repmat(lon(:)',numel(vstp),1);
        lat = repmat(lat(:)',numel(vstp),1);
        time = repmat(time(:)',numel(vstp),1);
        vbeg = repmat(vbeg(:),1,size(lon,2));
        vend = repmat(vend(:),1,size(lon,2));
        vstp = repmat(vstp(:),1,size(lon,2));
        nout = round((vend-vbeg)./vstp)+1;

        nalt = sum(nout(:,1),1);
        nummax = max(nout(:));
        nper = size(nout,1);
        nprof = size(nout,2);

        outfIndex = false(nummax,nprof.*nper);
        for i = 1:numel(nout)
            outfIndex(1:nout(i),i) = true;
        end
        outfIndex = find(outfIndex);
        %         outfIndex = [(1:nalt)']+[(0:(nprof-1)).*nummax.*nper];
        OutFFmt = @(X) X(outfIndex);
        oarIndex = 1:nper:(nper*nprof);

        oarFmt = @(X,io) reshape(X(io,oarIndex),nprof,1);

    elseif strcmpi(mode,'sat')
        tvec=[numel(lat),numel(lon),numel(alt),numel(time)];
        tvec=tvec==tvec';
        if ~all(tvec(:))
            error(ErrorStr,[...
                'In "sat" mode lat, lon, alt, time must all have the same size'])
        end

        Out.dates = time(:)';
        Out.lon = lon(:);
        Out.lat = lat(:);
        Out.alt = alt(:);
        vbeg=alt(:);
        vend=alt(:);
        vstp=alt(:);

        OutFFmt = @(X) X(:);
        oarFmt = @(X,io) X(io,:)';

    elseif strcmpi(mode,'map')

        lat = unique(lat);
        lon = unique(lon);
        alt = unique(alt);
        time = unique(time);

        vbeg=[];
        vend=[];
        vstp=[];

        split_alt;

        nlon = numel(lon);
        nlat=numel(lat);
        nalt=numel(alt);
        ntime=numel(time);


        Out.lat = lat(:);
        Out.lon = lon(:);
        Out.alt = alt(:);
        Out.dates = time(:)';

        [vbeg,~,~,~] = ndgrid(vbeg(:),lat(:),lon(:),time(:));
        [vend,~,~,~] = ndgrid(vend(:),lat(:),lon(:),time(:));
        [vstp,lat,lon,time] = ndgrid(vstp(:),lat(:),lon(:),time(:));


        nout = ceil((vend-vbeg)./vstp)+1;

        nalt = sum(nout(:,1),1);
        nummax = max(nout(:));
        nper = size(nout,1);
        nprof = nlon*nlat*ntime;

        outfIndex = [(1:nalt)']+[(0:(nprof-1)).*nummax.*nper];


        OutFFmt = @(X) permute(reshape(X(outfIndex),nalt,nlat,nlon,ntime),[4,2,3,1]);
        oarIndex = 1:nper:(nper*nprof);

        oarFmt = @(X,io) permute(reshape(X(io,oarIndex),nlat,nlon,ntime),[3,1,2]);


    end

    noutf = 0;
    if JF(1) > 0; noutf = noutf+1; end
    if JF(2) > 0; noutf = noutf+3; end
    if JF(3) > 0; noutf = noutf+5; end
    if JF(3) > 0 && JF(6) < 0; noutf = noutf+2; end
    if JF(24) < 0; noutf = noutf+1; end

    % no user inputs
    JF([8,9,10,13,14,15,16,17,25,27,32,43,44,45,46]) = 1;



    % no outputs but our outputs
    JF([34,38]) = -1;

    UT = 24.*mod(datenum(time(:)),1);
    mmdd = 100.*time(:).Month + time(:).Day;
    yr = time(:).Year;

    ind = [zeros(size(lat(:))),lat(:),lon(:),yr(:),mmdd(:),UT(:)+25,vbeg(:),vend(:),vstp(:)];
    UInput = false;
    if ~isempty(F107)
        JF(25) = -1;
        UInput=true;
    end

    if ~isempty(F107_81)
        JF(32) = -1;
        UInput=true;
    end

    if UInput
        UInput = zeros(size(ind,1),16);
        if ~isempty(F107)
            UInput(:,10) = F107;
        end
        if ~isempty(F107_81)
            UInput(:,12) = F107_81;
        end

        [outf,oar] = mex_iri2020(double(JF),double(ind'), double(UInput'));

    else
        [outf,oar] = mex_iri2020(double(JF),double(ind'));

    end
    k = 1;

    if JF(1) > 0 % output electron density
        Out.dens = OutFFmt(outf(:,k:noutf:end));
        k = k+1;
    end

    if JF(2) > 0% output temperatures
        fn = {'Tn','Ti','Te'};
        for i = 1:numel(fn)
            Out.(fn{i}) = OutFFmt(outf(:,k:noutf:end));
            k = k+1;
        end
    end

    if JF(3) > 0% output ion densities
        fn = {'O_p','H_p','He_p','O2_p','NO_p'};
        for i = 1:numel(fn)
            Out.(fn{i}) = OutFFmt(outf(:,k:noutf:end));
            k = k+1;
        end

        if JF(6) < 0% output other ions
            fn = {'Cluster_p','N_p'};
            for i = 1:numel(fn)
                Out.(fn{i}) = OutFFmt(outf(:,k:noutf:end));
                k = k+1;
            end
        end
    end


    iid=[     1,    2,      3,     4,    5,    6,  10, 35, ...
        23, 25, 36, 33, 39, 41, 46, 83];
    iin={'NmF2','hmF2','NmF1','hmF1','NmE','hmE','B0','B1', ...
        'solzen','dip','M3000F2','Rz12','IG12','F107','F107_81','Kp'};
    for ji = 1:numel(iid)
        Out.(iin{ji}) = oarFmt(oar,iid(ji));
        Out.(iin{ji})(Out.(iin{ji})==-1)=nan;
    end

    Out.mlat = oarFmt(oar,49);
    Out.mlon = oarFmt(oar,50);
    Out.mlt  = oarFmt(oar,54);

    %            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS
    %
    %      #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
    %      #OARR(3) = NMF1/M-3           #OARR(4) = HMF1/KM
    %      #OARR(5) = NME/M-3            #OARR(6) = HME/KM
    %       OARR(7) = NMD/M-3             OARR(8) = HMD/KM
    %       OARR(9) = HHALF/KM           #OARR(10) = B0/KM
    %       OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
    %       OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
    %      #OARR(15) = TE-MOD(300KM)     #OARR(16) = TE-MOD(400KM)/K
    %       OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
    %       OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
    %       OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
    %       OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
    %       OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
    %       OARR(27) = MODIFIED DIP LAT.  OARR(28) = Geographic latitude
    %       OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
    %       OARR(31) = ISEASON (1=spring) OARR(32) = Geographic longitude
    %      #OARR(33) = Rz12               OARR(34) = Covington Index
    %      #OARR(35) = B1                 OARR(36) = M(3000)F2
    %      $OARR(37) = TEC/m-2           $OARR(38) = TEC_top/TEC*100.
    %      #OARR(39) = gind (IG12)        OARR(40) = F1 probability
    %      #OARR(41) = F10.7 daily        OARR(42) = c1 (F1 shape)
    %       OARR(43) = daynr              OARR(44) = equatorial vertical
    %       OARR(45) = foF2_storm/foF2_quiet         ion drift in m/s
    %      #OARR(46) = F10.7_81           OARR(47) = foE_storm/foE_quiet
    %       OARR(48) = spread-F probability
    %       OARR(49) = Geomag. latitude   OARR(50) = Geomag. longitude
    %       OARR(51) = ap at current time OARR(52) = daily ap
    %       OARR(53) = invdip/degree      OARR(54) = MLT-Te
    %       OARR(55) = CGM-latitude       OARR(56) = CGM-longitude
    %       OARR(57) = CGM-MLT            OARR(58) = CGM lat eq. aurl bodry
    %       OARR(59) = CGM-lati(MLT=0)    OARR(60) = CGM-lati for MLT=1
    %       OARR(61) = CGM-lati(MLT=2)    OARR(62) = CGM-lati for MLT=3
    %       OARR(63) = CGM-lati(MLT=4)    OARR(64) = CGM-lati for MLT=5
    %       OARR(65) = CGM-lati(MLT=6)    OARR(66) = CGM-lati for MLT=7
    %       OARR(67) = CGM-lati(MLT=8)    OARR(68) = CGM-lati for MLT=9
    %       OARR(69) = CGM-lati(MLT=10)   OARR(70) = CGM-lati for MLT=11
    %       OARR(71) = CGM-lati(MLT=12)   OARR(72) = CGM-lati for MLT=13
    %       OARR(73) = CGM-lati(MLT=14)   OARR(74) = CGM-lati for MLT=15
    %       OARR(75) = CGM-lati(MLT=16)   OARR(76) = CGM-lati for MLT=17
    %       OARR(77) = CGM-lati(MLT=18)   OARR(78) = CGM-lati for MLT=19
    %       OARR(79) = CGM-lati(MLT=20)   OARR(80) = CGM-lati for MLT=21
    %       OARR(81) = CGM-lati(MLT=22)   OARR(82) = CGM-lati for MLT=23
    %       OARR(83) = Kp at current time OARR(84) = magnetic declination
    %       OARR(85) = L-value            OARR(86) = dipole moment
    %       OARR(87) = SAX300             OARR(88) = SUX300
    %      #OARR(89) = HNEA              #OARR(90) = HNEE

    Out.JF = JF;

    cd(OldDir);

catch err

    cd(OldDir)
    rethrow(err)

end


    function split_alt
        if numel(alt)>1
            dalt = diff(alt(:)');
            dalt(abs(dalt)<1e-3) = 0;
            hind=find(abs([1,diff(dalt),1])>1e-3);
            hind(1)=0;
            hind([0,diff(hind)]==1)=[];

            while any(diff(hind)>1000)
                ib = find(diff(hind)>1000,1);
                hind=[hind(1:ib),hind(ib)+999,hind((ib+1):end)];
            end

            vbeg = alt(hind(1:end-1)+1);
            vend = alt(hind(2:end));
            vstp = dalt(hind(1:end-1)+1);
        else
            vbeg=alt;
            vend=alt;
            vstp=1;
        end
    end

end