function data = lidarReadFcn(filename)
    tmp = load(filename);

    data.N_o3 = tmp.N_o3_fr_mg;
    data.on = tmp.mg_on_fr;
    data.off = tmp.mg_off_fr;
end