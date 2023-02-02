function apf = read_apf(file, update, remote_url)
%READ_APF read and update an apf107.dat file as needed for the IRI
%
%   Syntax
%   ------
%   APF = READ_APF(FILE)
%   APF = READ_APF(FILE,UPDATE)
%   APF = READ_APF(FILE,UPDATE,REMOTE_URL)
%
%   Description
%   ------------
%   APF = READ_APF(FILE) will read the solar and magnetic indices given in
%   an apf107.dat file as needed for geophysical models like the
%   International Reference Ionosphere (IRI). If no path is given, it will
%   look in the current folder by default.
%
%   APF = READ_APF(FILE,UPDATE) will attempt to update the existing
%   apf107.dat file from a remote source. If the file has already been
%   updated today, it will not re-download it. If the download does not
%   complete successfully, the existing data will not be overwritten.
%
%   APF = READ_APF(FILE,UPDATE,REMOTE_URL) will use a remote source other 
%   than the default
%
%   Input Arguments
%   ---------------
%
%   Name  Description                                     Data Type
%   ----  --------------------                            ---------
%   FILE  Path to the apf107.dat file                     string
%
%   UPDATE
%         Option to update local file                     boolean
%
%   REMOTE_URL
%         URL of remote source data file                  string
%
%   See also http://irimodel.org/indices/IRI-Format-indices-files.pdf

if ~exist('file', 'var') || isempty(file)
    file = fullfile(pwd,'apf107.dat');
end

if ~exist('update', 'var') || isempty(update)
    update = false;
end

if ~exist('remote_url', 'var') || isempty(remote_url)
    remote_url = 'https://chain-new.chain-project.net/echaim_downloads/apf107.dat';
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

apf = [];

while update
    update = false;

    try
        data = webread(remote_url);
    catch Err
        warning('apf107:DownloadError', ...
            'While downloading the indices file the following error occured:\n %s', Err.message)
        break
    end

    try
        apf = t_read_apf(data');
    catch Err
        warning('apf107:CorruptedUpdate',...
            'While parsing the updated indices file the following error occured:\n %s', Err.message)
        break
    end

    fid = fopen(file,'w');
    fprintf(fid,'%s',data);
    fclose(fid);

end

if isempty(apf) % did not update file
    apf = t_read_apf();
end


    function C = t_read_apf(fp)

        formatSpec = '%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%*3d%5.1f%5.1f%5.1f';
        formatSize = cellfun(@str2double,regexp(formatSpec,'(?<=\%\**)[\d]*','match'),'UniformOutput',true);

        dw=arrayfun(@(n,m) (1:n) + m,formatSize,[0,cumsum(formatSize(1:end-1))],'UniformOutput',false);

        if ~exist('fp','var')
            % reads apf107.dat
            fp = fopen(file,'r');
            D  = fread(fp,inf,'uint8=>uint8');
            fclose(fp);
        else
            D = fp;
        end

        D  = char(strread(char(D),'%s','delimiter','\n','whitespace',''));

        del = repmat(',',size(D,1),1);

        D = cell2mat(cellfun(@(n) [D(:,n),del],dw,'UniformOutput',false));
        C = textscan(D',formatSpec,'Delimiter',',');

        C = cell2struct(C,{'year','month','day','Ap03','Ap06','Ap09','Ap12','Ap15','Ap18','Ap21','Ap24','Ap','F107','F107_81','F107_365'},2);
        C.year=double(C.year)+1900+((double(C.year)<50).*100);
        C.date =  datetime(C.year,double(C.month),double(C.day));
        C = rmfield(C,{'year','month','day'});
    end

end