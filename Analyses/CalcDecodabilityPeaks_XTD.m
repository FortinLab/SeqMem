%% Calc Decoding Peaks
function [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodabilityPeaks_XTD(data)
tm_PeakNdx = cell(size(data,1), size(data,1), size(data,1));
tm_PeakVal = cell(size(data,1), size(data,1), size(data,1));
tm_PeakWid = cell(size(data,1), size(data,1), size(data,1));
for odrPos = 1:size(data,1)
    for histOSpos = 1:size(data,1)
        for trlPos = 1:size(data,1)
            temp_Data = squeeze(data(odrPos,histOSpos,trlPos,:));
            if ~isempty(temp_Data)
                temp_PeakNdx = nan(size(data,1),size(temp_Data{1},2),size(temp_Data{1},3));
                temp_PeakVal = nan(size(data,1),size(temp_Data{1},2),size(temp_Data{1},3));
                temp_PeakWid = nan(size(data,1),size(temp_Data{1},2),size(temp_Data{1},3));
                for pos = 1:size(temp_Data,1)
                    cur_temp_Data = temp_Data{pos};
                    for t = 1:size(cur_temp_Data,2)
                        for trl = 1:size(cur_temp_Data,3)
                            [pks,loc,wid,prom] = findpeaks(cur_temp_Data(:,t));
                            if ~isempty(pks)
                                featWeightPeaks = pks.*wid.*prom;
                                temp_PeakNdx(pos,t,trl) = loc(featWeightPeaks==max(featWeightPeaks));
                                temp_PeakVal(pos,t,trl) = pks(featWeightPeaks==max(featWeightPeaks));
                                temp_PeakWid(pos,t,trl) = wid(featWeightPeaks==max(featWeightPeaks));
                            end
                        end
                    end
                end
                tm_PeakNdx{odrPos,histOSpos,trlPos} = temp_PeakNdx;
                tm_PeakVal{odrPos,histOSpos,trlPos} = temp_PeakVal;
                tm_PeakWid{odrPos,histOSpos,trlPos} = temp_PeakWid;
            end
        end
    end
end
end