function [ out ] = data_in(Conn, node,sihi)
% data_in: generate output object of time and data from some node of a
% connection object initiated tree. `Sihi' initiatites sihi smoothing,
% where sihi is the injector frequency if not turned off (sihi=0)
try
out.t = Conn.get(['dim_of(' node ')']).getDoubleArray;

out.y = Conn.get(node).getDoubleArray;
if sihi ~= 0 
    out.y = sihi_smooth(out.y,out.t,sihi);
end


catch
    display(['NO DATA FOUND, NODE: ' node]);
    out.t=NaN;
    out.y=NaN;
end

end

