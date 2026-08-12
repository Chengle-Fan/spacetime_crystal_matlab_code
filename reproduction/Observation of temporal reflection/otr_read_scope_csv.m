function scope = otr_read_scope_csv(filename)
%OTR_READ_SCOPE_CSV Read Tektronix CSV exports with embedded metadata.
%   SCOPE = OTR_READ_SCOPE_CSV(FILENAME) locates the TIME/CH* header rather
%   than assuming a fixed metadata length.  It also supports exports such as
%   Fig. 2c source data, where a second TIME/MATH (frequency spectrum) block
%   is placed to the right of an empty separator column.
%
%   Fields of SCOPE:
%       time, channels, channelNames       primary time-domain block
%       frequency, spectra, spectrumNames optional right-hand block
%       sampleInterval, metadata, rawData, rawHeaders, filename
%
%   No toolbox functions are used; compatible with MATLAB R2019b+.

    if nargin < 1 || isempty(filename)
        error('otr:MissingFilename', 'A Tektronix CSV filename is required.');
    end
    if exist(filename, 'file') ~= 2
        error('otr:FileNotFound', 'Cannot find scope CSV: %s', filename);
    end

    fid = fopen(filename, 'rt');
    if fid < 0
        error('otr:FileOpenFailed', 'Cannot open scope CSV: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    metadata = cell(0, 1);
    headers = {};
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if ~isempty(line) && double(line(1)) == 65279
            line(1) = [];
        elseif numel(line) >= 3 && isequal(double(line(1:3)), [239 187 191])
            line(1:3) = [];
        end
        fields = regexp(line, ',', 'split');
        first = '';
        if ~isempty(fields)
            first = upper(strtrim(strrep(fields{1}, '"', '')));
        end
        if strcmp(first, 'TIME')
            headers = fields;
            break;
        end
        metadata{end + 1, 1} = line; %#ok<AGROW>
    end
    if isempty(headers)
        error('otr:ScopeHeaderMissing', 'No TIME header found in %s.', filename);
    end

    headers = cellfun(@(value) strtrim(strrep(value,'"','')),headers, ...
        'UniformOutput',false);
    % Tektronix exports in the paper contain 100,000 rows.  TEXTSCAN keeps
    % this reader toolbox-free while avoiding a slow cell-array append and
    % per-token SSCAN loop for every record.
    dataPosition = ftell(fid);
    firstDataLine = fgetl(fid);
    if ~ischar(firstDataLine)
        error('otr:NoScopeData', 'No numeric scope records found in %s.', filename);
    end
    firstFields = regexp(firstDataLine, ',', 'split');
    maxColumns = max(numel(headers),numel(firstFields));
    fseek(fid,dataPosition,'bof');
    format = repmat('%f',1,maxColumns);
    try
        parsed = textscan(fid,format,'Delimiter',',', ...
            'CollectOutput',true,'TreatAsEmpty',{'','""'}, ...
            'EmptyValue',NaN,'ReturnOnError',false);
        raw = parsed{1};
    catch
        % Older releases do not accept an empty token in TreatAsEmpty;
        % EmptyValue already handles consecutive delimiters.
        fseek(fid,dataPosition,'bof');
        parsed = textscan(fid,format,'Delimiter',',', ...
            'CollectOutput',true,'EmptyValue',NaN,'ReturnOnError',false);
        raw = parsed{1};
    end
    if isempty(raw)
        error('otr:NoScopeData', 'No numeric scope records found in %s.', filename);
    end
    headers(end + 1:maxColumns) = {''};

    primaryTime = find(strcmpi(headers, 'TIME'), 1, 'first');
    if isempty(primaryTime)
        primaryTime = 1;
    end
    rightTime = find(strcmpi(headers, 'TIME'));
    rightTime = rightTime(rightTime > primaryTime);
    if isempty(rightTime)
        primaryEnd = maxColumns;
    else
        primaryEnd = rightTime(1) - 1;
    end

    channelColumns = [];
    for jj = primaryTime + 1:primaryEnd
        name = upper(headers{jj});
        if strncmp(name, 'CH', 2) || (~isempty(name) && any(isfinite(raw(:, jj))))
            channelColumns(end + 1) = jj; %#ok<AGROW>
        end
    end
    validTime = isfinite(raw(:, primaryTime));
    time = raw(validTime, primaryTime);
    channels = raw(validTime, channelColumns);
    channelNames = headers(channelColumns);

    frequency = zeros(0, 1);
    spectra = zeros(0, 0);
    spectrumNames = cell(1, 0);
    if ~isempty(rightTime)
        freqColumn = rightTime(1);
        spectrumColumns = [];
        for jj = freqColumn + 1:maxColumns
            if any(isfinite(raw(:, jj)))
                spectrumColumns(end + 1) = jj; %#ok<AGROW>
            end
        end
        validFrequency = isfinite(raw(:, freqColumn));
        frequency = raw(validFrequency, freqColumn);
        spectra = raw(validFrequency, spectrumColumns);
        spectrumNames = headers(spectrumColumns);
    end

    scope = struct();
    scope.time = time;
    scope.channels = channels;
    scope.channelNames = channelNames;
    scope.frequency = frequency;
    scope.spectra = spectra;
    scope.spectrumNames = spectrumNames;
    scope.sampleInterval = local_sample_interval(time);
    scope.metadata = metadata;
    scope.rawData = raw;
    scope.rawHeaders = headers;
    scope.filename = filename;
end

function dt = local_sample_interval(time)
    delta = diff(time(:));
    delta = delta(isfinite(delta) & delta > 0);
    if isempty(delta)
        dt = NaN;
    else
        dt = median(delta);
    end
end
