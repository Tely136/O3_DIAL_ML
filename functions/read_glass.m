function ds = read_glass(file)
    dt = hdfread(file,"DATETIME");
    dt_start = hdfread(file,"DATETIME.START");
    dt_stop = hdfread(file,"DATETIME.STOP");

    ds.dt = mj2k(dt);
    ds.dt_start = mj2k(dt_start);
    ds.dt_stop = mj2k(dt_stop);

    ds.alt = double(hdfread(file,"ALTITUDE"));
    ds.alt_inst = hdfread(file,"ALTITUDE.INSTRUMENT");
    ds.N = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL");
    ds.p = hdfread(file,"PRESSURE_INDEPENDENT"); % hPa
    ds.t = hdfread(file,"TEMPERATURE_INDEPENDENT"); % k
end


function out = mj2k(in)

out = datetime(round(in .* 24 .* 60 .* 60),'ConvertFrom','epochtime','Epoch','2000-01-01');

end