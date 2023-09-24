function dataStruct = b4_collateRs(dataStruct,timepoint,param,dataType,paramVals,rawAmplitude)

    if ~isfield(dataStruct,'timepoint');dataStruct.timepoint=[];end
    dataStruct.timepoint = [dataStruct.timepoint,timepoint];
    if ~isfield(dataStruct,'param');dataStruct.param=[];end
    dataStruct.param = [dataStruct.param,param];
    if ~isfield(dataStruct,'dataType');dataStruct.dataType=[];end
    dataStruct.dataType = [dataStruct.dataType,dataType];
    if ~isfield(dataStruct,'r');dataStruct.r=[];end
    dataStruct.r = [dataStruct.r,corr(paramVals', rawAmplitude')];
    
    if ~isfield(dataStruct,'paramVals');dataStruct.paramVals={};end
    dataStruct.paramVals = [dataStruct.paramVals,{paramVals}];
    if ~isfield(dataStruct,'rawAmplitude');dataStruct.rawAmplitude={};end
    dataStruct.rawAmplitude = [dataStruct.rawAmplitude,{rawAmplitude}];

end
