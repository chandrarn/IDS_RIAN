function [ HRES ] = PhGetCineAuxData( CH, ImageNumber, SelectionCode, DataSize, pData)
% pData is a libpointer to the desired information
[HRES, dummyAll] = calllib('phfile','PhGetCineAuxData', CH, ImageNumber, SelectionCode, DataSize, pData);
OutputError(HRES);
end
