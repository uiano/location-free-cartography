[fileName,fileAddress] = uiputfile('*.mat');

MATLABStructMode = 1;
    if exist('measurementLocationXyz','var')
        save([fileAddress,fileName],'wall','ceillFloor','ceilingLevel','groundLevel','Tx','boundary','mesh_','Rx','measurementLocationXyz','MATLABStructMode')
    else
%         measurementLocationXyz = [];
        errordlg('Could not find ''measurementLocationXyz'' in workspace.',...
            'Variable not found! This variable saved empty.');
        save([fileAddress,fileName],'wall','ceillFloor','ceilingLevel','groundLevel','Tx','boundary','mesh_','Rx','MATLABStructMode')
    end
    
    disp('File saved sucessfully.')
% elseif strcmpi(saveType,'n')
%     save([fileAddress,fileName],'wall','ceillFloor','ceilingLevel','groundLevel','Tx','boundary','mesh_','Rx')
% end