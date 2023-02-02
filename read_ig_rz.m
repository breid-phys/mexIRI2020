function ig_rz = read_ig_rz(file, update, remote_url)
%READ_IG_RZ read and update an ig_rz.dat file as needed for the IRI
%
%   Syntax
%   ------
%   IG_RZ = READ_IG_RZ(FILE)
%   IG_RZ = READ_IG_RZ(FILE,UPDATE)
%   IG_RZ = READ_IG_RZ(FILE,UPDATE,REMOTE_URL)
%
%   Description
%   ------------
%   IG_RZ = READ_IG_RZ(FILE) will read the solar and magnetic indices given in
%   an ig_rz.dat file as needed for geophysical models like the
%   International Reference Ionosphere (IRI). If no path is given, it will
%   look in the current folder by default.
%
%   IG_RZ = READ_IG_RZ(FILE,UPDATE) will attempt to update the existing
%   ig_rz.dat file from a remote source. If the file has already been
%   updated today, it will not re-download it. If the download does not
%   complete successfully, the existing data will not be overwritten.
%
%   IG_RZ = READ_IG_RZ(FILE,UPDATE,REMOTE_URL) will use a remote source other 
%   than the default
%
%   Input Arguments
%   ---------------
%
%   Name  Description                                     Data Type
%   ----  --------------------                            ---------
%   FILE  Path to the ig_rz.dat file                     string
%
%   UPDATE
%         Option to update local file                     boolean
%
%   REMOTE_URL
%         URL of remote source data file                  string
%
%   See also http://irimodel.org/indices/IRI-Format-indices-files.pdf

if ~exist('file', 'var') || isempty(file)
    file = fullfile(pwd,'ig_rz.dat');
end

if ~exist('update', 'var') || isempty(update)
    update = false;
end

if ~exist('remote_url', 'var') || isempty(remote_url)
    remote_url = 'https://chain-new.chain-project.net/echaim_downloads/ig_rz.dat';
end


if update % don't update if file already updated today
    if exist(file,'file')
        fdat=dir(file);
        if fdat.datenum<(now-1)
            update = true;
        else
            update = false;
        end
    end
end

ig_rz = [];

while update
    update = false;

    try
        data = webread(remote_url);
    catch Err
        warning('ig_rz:DownloadError', ...
            'While downloading the indices file the following error occured:\n %s', Err.message)
        break
    end

    try
        ig_rz = t_read_ig_rz(data');
    catch Err
        warning('ig_rz:CorruptedUpdate',...
            'While parsing the updated indices file the following error occured:\n %s', Err.message)
        break
    end

    fid = fopen(file,'w');
    fprintf(fid,'%s',data);
    fclose(fid);

end

if isempty(ig_rz) % did not update file
    ig_rz = t_read_ig_rz();
end


    function C = t_read_ig_rz(fp)

        if ~exist('fp','var')
            % reads ig_rz.dat
            fp = fopen(file,'r');
            D  = fread(fp,inf,'uint8=>uint8');
            fclose(fp);
        else
            D = fp;

        end

        D  = char(strread(char(D),'%s','delimiter','\n','whitespace',''));
        Dates = strsplit(D(3,:),',');
        DateStart = datetime(str2num(Dates{2}),str2num(Dates{1}),1)-days(14);
        DateEnd = datetime(str2num(Dates{4}),str2num(Dates{3}),1) + days(28+14);
        Weeks = DateStart:days(7):DateEnd;
        Months = unique([Weeks.Year(:),Weeks.Month(:)],'rows');
        Q = str2num(reshape(D(5:end,:)',1,[]));
        Q = Q(:);

        C.IG12 = Q(1:size(Months,1));
        C.Rz12 = Q(size(Months,1) + (1:size(Months,1)));
        C.date = datetime(Months(:,1),Months(:,2),ones(size(Months(:,1))));

    end

end