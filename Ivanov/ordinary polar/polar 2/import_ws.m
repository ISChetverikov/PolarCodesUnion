function import_ws(varargin)
    for i=1:length(varargin)
        evalin('caller', [varargin{i} ' = evalin(''base'', ''' varargin{i} ''');']);
    end
end
