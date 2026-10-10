%[text] checking 'create\_new\_database' reads GLASS HDF files in specified folder
%[text]     creates table with time of data, file location, and id's of good profiles
%[text]     checking rm\_cld will filter out profiles with clouds from good profiles
%[text]     the table is save to specified location to be used later
clc; clearvars;

addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\functions")); %[control:filebrowser:470d]{"position":[17,84]}
addpath(genpath("C:\Users\tely1\OneDrive\Matlab Onedrive\O3_Lidar\proc_functions")); %[control:filebrowser:41bc]{"position":[17,82]}
addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\functions")); %[control:filebrowser:5b94]{"position":[17,79]}

folder = "C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\GLASS"; %[control:filebrowser:0877]{"position":[10,68]}

files = dir(fullfile(folder,"*.hdf"));

info = hdfinfo(fullfile(folder,files(1).name));
n = length(info.SDS);
for i = 1:n %[output:group:856bc5d8]
    SDS_temp = info.SDS(i);

    disp(SDS_temp.Name) %[output:0e9341a7]
end %[output:group:856bc5d8]
%%
n_files = length(files);

save_name = "all_glass_profiles.mat"; %[control:editfield:1521]{"position":[13,37]}
save_folder = "C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_DIAL_ML\dir"; %[control:filebrowser:2ba8]{"position":[15,71]}

create_new_database = true; %[control:checkbox:782a]{"position":[23,27]}
rm_cld = false; %[control:checkbox:1db7]{"position":[10,15]}

if create_new_database==1
    glass_dataset = table;

    count = 0;
    for k = 1:n_files
        file = fullfile(folder, files(k).name);
    
        dt = hdfread(file,"DATETIME");
        dt_start = hdfread(file,"DATETIME.START");
        dt_stop = hdfread(file,"DATETIME.STOP");
        N = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL");
        p = hdfread(file,"PRESSURE_INDEPENDENT"); % hPa
        t = hdfread(file,"TEMPERATURE_INDEPENDENT"); % k
    
        p = p';
        t = t';
    
        N(N<0) = NaN;
        N = N';
        
        max_row = size(N,1);
        
        N_good = NaN(size(N));
        n_prof = size(N,2);
        id_good = [];
        
        ct_good = 0;
        for i = 1:n_prof
            temp_prof = N(:,i);
    
            if all(isnan(temp_prof))
                continue
            end
        
            good = true;
            if rm_cld==1
                id_nan = find(isnan(temp_prof));
                for j = 1:length(id_nan)
                    if id_nan(j) == max_row
                        break
                    end
            
                    if ~isnan(temp_prof(id_nan(j)+1))
                        good = false;
                        break;
                    end
                end
            end
            
            if good == true
                N_good(:,ct_good+1) = temp_prof;
                id_good(end+1) = i;
                ct_good = ct_good + 1;
            end
        end
    
        if ct_good==0
            continue
        end
        count = count+1;

        id_good = {id_good};
        time = dt(1);

        dataset_temp = table(time, file, id_good);
        glass_dataset = [glass_dataset; dataset_temp];
    end
    save(fullfile(save_folder,save_name),"glass_dataset");

else
    load(fullfile(save_folder,save_name))
end
%%
id =18; %[control:slider:4d1b]{"position":[5,7]}

test_glass = read_glass(glass_dataset.file(id));
idx_good = cell2mat(glass_dataset.id_good(id));

N_glass = test_glass.N';
alt = test_glass.alt;
alt_inst = test_glass.alt_inst;
time = test_glass.dt;

N_glass(N_glass<0) = NaN;
%%
n_prof = length(idx_good);
idx = 6; %[control:slider:842b]{"position":[7,8]}
prof = idx_good(idx);
N_glass_prof = N_glass(:,prof);

figure; %[output:1cd59973]
plot(N_glass_prof,alt) %[output:1cd59973]
grid on %[output:1cd59973]
%%
% Extend altitude grid and ozone vector down to lidar
alt_lower = alt_inst:30:min(alt);
alt_ext = [alt_lower alt];

N_ext = [zeros(length(alt_lower), n_prof); N_glass(:,idx_good)];
%%
% Fill missing values in GLASS data
N_ext = fillmissing(N_ext,"linear",1,"EndValues",0);
%%
% Ozone Absorption Cross Section
BDM = load("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\cross_sections\ozone\bdm\bdm_287_299.mat"); %[control:filebrowser:0fea]{"position":[12,110]}

o3_287 = BDM.BDM_c287(1) .* 1e-4; % use first temp dependent cross section for now
o3_299 = BDM.BDM_c299(1) .* 1e-4;

do3 = BDM.BDM_d_287_299(1) .* 1e-4;
%%
% Air Number Density
Nair = airNumberDensity(alt_ext);

sigma_m_on = rayleigh_cross_section(287.2e-9);   
sigma_m_off = rayleigh_cross_section(299.1e-9);   
dsigma_m = sigma_m_on - sigma_m_off;

D_No3 = dsigma_m .* Nair ./ do3;
%%
% Prepare parameters for signal calculation
E0_287 = 10.6e-3;
E0_299 = 6.3e-3;

eta_287 = 1/53.1572075; % Efficiency factor
eta_299 = 1/47.0176125; % 

etaQ = 0.7;

t_r = .2;   % Radius of the receiver in meters
o_r = 0.1;  % Radius of the obstruction
tau = 6e-9; % Laser pulse duration

% Calculate the signal using the defined parameters
[N_287,K_287,G_287,B_287,T_287] = signal(alt_ext,alt_inst,287.2e-9,E0_287,t_r,o_r,tau,eta_287,etaQ,Nair,o3_287,N_ext);  %[output:693537c4]
[N_299,K_299,G_299,B_299,T_299] = signal(alt_ext,alt_inst,299.1e-9,E0_299,t_r,o_r,tau,eta_299,etaQ,Nair,o3_299,N_ext); 

N_287(1:length(alt_lower),:) = [];
N_299(1:length(alt_lower),:) = [];

save("./sim_prof.mat","N_287","N_299","alt");
%%
figure;
plot(alt_ext,T_287(:,idx))
xlabel("altitude (m)")
ylabel("Transmittance")
grid on
figure;
semilogy(alt,N_287(:,idx),'DisplayName','287')
hold on
semilogy(alt,N_299(:,idx),'DisplayName','299')
grid on
legend
xlim([0 2e3])
%%
dr = NaN(length(alt)-1,1);
for i = 1:length(alt)-1
    dr(i) = alt(i+1) - alt(i);
end

s = log(N_299./N_287);
S = diff(s)./dr;
N_sim = S./(2*do3)' - D_No3(length(alt_lower)+1:end-1)';

figure;
plot(N_sim(:,idx),alt(1:end-1),'DisplayName','Simulated','LineWidth',2)
hold on
plot(N_glass_prof,alt,'DisplayName','GLASS','LineStyle','--','LineWidth',2)
legend
%%
Clim = [.6e18 2e18];

figure
imagesc(time(:,idx_good),alt,N_glass(:,idx_good))
colorbar
set(gca,'YDir','normal')
xlabel('Time (UTC)')
title('GLASS')
clim(Clim)
figure
imagesc(time(:,idx),alt,N_sim)
colorbar
set(gca,'YDir','normal')
xlabel('Time (UTC)')
title('Simulated')
clim(Clim)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[control:filebrowser:470d]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:41bc]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:5b94]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:filebrowser:0877]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:editfield:1521]
%   data: {"defaultValue":"\"\"","label":"save_name","run":"Nothing","valueType":"String"}
%---
%[control:filebrowser:2ba8]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[control:checkbox:782a]
%   data: {"defaultValue":false,"label":"crete_new_dataset","run":"Nothing"}
%---
%[control:checkbox:1db7]
%   data: {"defaultValue":false,"label":"rm_cld","run":"Nothing"}
%---
%[control:slider:4d1b]
%   data: {"defaultValue":1,"label":"Slider","max":162,"maxLinkedVariable":"n_files","min":1,"run":"SectionToEnd","runOn":"ValueChanging","step":1}
%---
%[control:slider:842b]
%   data: {"defaultValue":1,"label":"Slider","max":10,"maxLinkedVariable":"n_prof","min":1,"run":"SectionToEnd","runOn":"ValueChanged","step":1}
%---
%[control:filebrowser:0fea]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"","run":"Section"}
%---
%[output:0e9341a7]
%   data: {"dataType":"text","outputData":{"text":"LATITUDE.INSTRUMENT\nLONGITUDE.INSTRUMENT\nALTITUDE.INSTRUMENT\nDATETIME\nDATETIME.START\nDATETIME.STOP\nINTEGRATION.TIME\nALTITUDE\nO3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL\nO3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD\nO3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD\nO3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD\nO3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM\nO3.MIXING.RATIO.VOLUME_DERIVED\nO3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD\nO3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD\nO3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD\nPRESSURE_INDEPENDENT\nTEMPERATURE_INDEPENDENT\nPRESSURE_INDEPENDENT_SOURCE\nTEMPERATURE_INDEPENDENT_SOURCE\n","truncated":false}}
%---
%[output:1cd59973]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAQAElEQVR4Aey9DZBV13mm+wUuHIKAlrGGvyCJZohkEYKv+Mc3EtwpMddVbjEZ0IUqIocolhxqMGOoSYlcD+XIvpQqqFwFHkwKS3JkYg1VUGKSi1ozzkiVgMiEf1SWCbJkLj8W5q9kSQ3ouAUDnvMs9Xe8erMP3Q3dffbZ56V4z1p77b3P2etdu3s\/\/X1r79Nn6NChv5Lkgc4BnQM6B3QO6BzQOVBL50Af0z85IAfkgByQA3Kgiw5o82o7IICp9gjo8+WAHJADckAOyIEuOyCA6bJl2kEOyAE5UH0HdARyoN4dEMDU+xmg\/ssBOSAH5IAcqEEHBDA1OGg6ZDlQfQd0BHJADsiB6joggKmu\/\/p0OSAH5IAckANy4CYcEMDchGnapfoO6AjkgByQA3Kgvh0QwNT3+Kv3ckAOyAE5IAdq0gEBzE0Nm3aSA3JADsgBOSAHqumAAKaa7uuz5YAckANyQA7UkwPd2FcBTDeaqbeSA3JADsgBOSAHescBAUzv+KxPkQNyQA7Igeo7oCPIkQMCmBwNproiB+SAHJADcqBeHBDA1MtIq59yQA5U3wEdgRyQA93mgACm26zUG8kBOSAH5IAckAO95YAAprec1ufIgeo7oCOQA3JADuTGAQFMboZSHZEDckAOyAE5UD8OCGDqZ6yr31MdgRyQA3JADsiBbnJAANNNRupt5IAckANyQA7Igd5zINMA09jYaNu2bbPdu3cHrVmzpp0zCxYssJ07d4Z1bMf28QabNm0K69h\/2bJl8SrjvWhH1Nut1IIckANyQA7IATmQaQcyDTBPPfVUMG\/RokWGxo0bZ8uWLQttwMrChQttz549YR2NS5YsoQhiu1GjRtny5ctt8+bNNnfuXJs+fXpYB\/hMmjTJ1q5dG0SdtrBSL3JADsgBOSAH5EDkQDarmQUYYAMA2bFjhx0\/fjyI+pQpU4KTU6dOtSFDhtjBgwfL6wAcwIYN2O706dO2d+9ea25utgsXLti0adNYZZMnTw7L+\/fvt61btxrb0RZW6kUOyAE5IAfkgBzIvAOZBZhKzg0aNMiAlOHDh4dNTp48GUpeAJphw4aF9Wx39uxZmssaPXp0qI8YMcIuXboUwCc0lF5oKxX6LwfkgByQAxlzQIcjB9IcyCzAnD9\/PkRJZs+eHYAEaKHukEJniKqwHfVz587Z5cuXqZZ16tSpUCeCA7CEhbaXGG7ietvqdsXIkSONNJM0ST6UUo86D3Qe6BzQOVBL5wDXMMvhv8wCDNDxzDPPhDQRc1i+\/\/3v289\/\/vNeHwIGftWqVbZhwwZJHugcqMtzQD\/7+v1X2+cA1zCuZb1+Ae3hD8wswNBv5q\/MmTPHZs6cabNmzbJf\/vKXISrjUZc4GkNKqX\/\/\/uxWlqeMiN6QUiqvKFXilFFcL61q959Bh7RXr15tS5curQs9\/\/zzwYN66jNjq37Xzzmu8TbTz3f9\/D7nGsa1LPxiz9FLpgEm6TOgQSqI6ExayshTSqxnu+T+nlJKSxmltcX7nzlzxg4dOlQXeuONN0LX66nPjG1W+82x9aTU7\/r52eY80njX53iHX+o5e8k0wPAcF39GC7dFjxkzxl555ZUwBNxB1NraavPmzQtzZJgfc\/To0fLE3AMHDoT5GtzN1NTUFFJR+\/btC\/ty5xJ3OHHrNKJOW1ipFzkgB+SAHJADciDzDmQaYP7yL\/8yQAgPm+M5MC+99FK47RlXibJs2bLFZsyYEZ7zQtvGjRspgtavXx9uj163bl14Tsz27dvDLdWs5NZp\/hJZsWKFIeq0sU4ye\/\/99+2\/\/\/f\/Hkqz+nFE\/X6\/fga71FONt8a7dBrofw07kGmAiefAMA8GKIm9BjqYG8O6+fPnl6Mvvs3ixYvD\/BnWJ\/dduXJleR1130el2QcffBAAhrKe\/KC\/gBul+p1\/BxhnjXf+x9l7WK\/j7f3PY5lpgMFwSQ7IATkgB+SAHJADSQcEMElHtCwH5IAckANyoPYdyH0PBDC5H2J1UA7IATkgB+RA\/hwQwORvTNUjOSAH5ED1HdARyIEedkAA08MG6+3lgByQA3JADsiB7ndAANP9nuod5YAcqL4DOgI5IAdy7oAAJucDrO7JATkgB+SAHMijAwKYPI6q+lR9B3QEckAOyAE50KMOCGB61F69uRyQA3JADsgBOdATDghgesLV6r+njkAOyAE5IAfkQK4dEMDkenjVOTkgB+SAHJAD+XSgZwAmn16pV3JADsgBOSAH5EBGHBDAZGQgdBhyQA7IATkgB+RA5x0QwHTeK20pB+SAHJADckAOZMQBAUxGBkKHIQfkgByovgM6AjlQOw4IYGpnrHSkckAOyAE5IAfkQJsDApg2I1TIATlQfQd0BHJADsiBzjoggOmsU9pODsgBOSAH5IAcyIwDApjMDIUOpPoO6AjkgByQA3KgVhwQwNTKSOk45YAckANyQA7IgbIDmQeYTZs22e7du4Ool4+8VFmwYIHt3LkzrNu2bZs1NjaWWn\/9n+1932XLlv16Ram2Zs2asB\/rqZeaqv5fByAH5IAckANyQA50zoFMAwxgMWjQIFu0aJEtX77cRo0aZbTRNWBl4cKFtmfPnrCetiVLllAEASxsz36bN2+2uXPn2vTp08M6wGfSpEm2du3aIOq0hZV6kQNyQA7IgVw50KdhpF0ZMDSo312TrN9dk60wsSlo4ANPGBrU9OeGGv5go7k+tXS7sW+uzMhRZyKAyV6vRowYYZcuXbLjx4\/b3r177fTp00YbRzp16lQbMmSIHTx4MKzfsWOHjRs3rhyFmTJlStie\/Zqbm+3ChQs2bdo0drXJkyeH5f3799vWrVvDdrSFlXqRA3JADsiBzDgAQLgADxTDB9DhisED+Ljja\/sNDS2ByMAvbbVj079u9vAz1vDoRhtcAhZUmPhwCWQetn53T7a+JdCh41dbzhi6cvKgXSvVaZOy50CmAebs2bNGBIZoC9ETIioHDhwILg4fPjyUJ0+eDCUvAM2wYcMCxLAf+9PuGj16dKgCQQ5GoaH0QlupqPj\/\/vvvNyI1aMyYMdbQ0JBbDR482G677bYAiHnuZ7Jv6veQ3J7TybFmWePdu+M9cHijoRhACqUoSFr0IwkfAAgCPBDggRqmPmKD7\/mcDb3rM0GDrGho4IdHbcDJ163P6+vKGrzj\/7XPHH3BBr38Fbv64hftyoaHgq6V6i5rXmmo76515uJc6VAZux5wjeJahUaOHFnxulbrKzINMCtXrjQiK6SAnnnmGXvuueds\/fr1Zc+Jqpw\/fz4snzt3zi5fvhzq\/nLq1KlQJYIDsISFtpcYbuJ62+rriscff9w2bNgQ9Nhjj9mECRNyq\/Hjx9vYsWPtvvvuy20f08ZP\/dZ4p50XeWvrjvP83ikPGhr7e79vdz70x0EjvvCndsf8p4M8EkLUo9\/S1wxRRzGADH1omcUQckfhqqGhv3zXPn3ihzbi7c1Bd\/7oO4bG7v2moXt3LjdEHbEOjTmxzVxjP9htv\/2rE7\/WoFabMPSa3TOywcbf+elc\/27jGuXXq1WrVl13PctLQ6YBhkm4s2fPDnNc\/uiP\/siY8+JzYHp7AFavXm1Lly4NeuGFF+zw4cO51VtvvWXHjh2zI0eO5LaPaeOnfmu8086LGmjr0s9ppfP8yLu\/MPT2hf72098YY8c+NdNOjJkf9O5nv2KIFMzbs9aFVAx12s7eu8jQ+U991t7\/zTvtvY\/72vs\/+0lQy\/6X7GLzN4JaXlxiqPi9BYaIgFC6Wv7zEntv29fKOvvKt+zd1\/4q6Ng\/\/q2htw+8buhmxqRSv2\/mvbK+D9cov149\/\/zzvX257LXPyyzAeMqICAwRFESdkBjrcMhTRtRJKfXv359qWZ4yamxsDKmo8opSJU4ZxfXSqtT\/Z86csUOHDgWdOHHCWlpacisiWx999JFdvHgxt31MGz\/1W+Oddl7UettFG2io2DDOWu9+wC7d+\/v2\/4\/+N3b5\/1pt1rTG+jz6gxAh6VsqEXNErj24vLTtg1a8fZxdKu1fCUje3zDX3nt6atAHpToCRC6VoAUVdz1nH7\/ZHHTlZwcNFc8dN9TbvtbTzzfXKL9evfHGG6nXtDw0ZhZgOjI3LWXECUpKCdhJpox4P08ppaWM0trYR5IDcqADB7S6qg74BFfmkyAmtHoKxyew+vwR5o0MmNgUwOTj1tYwUfXjN19uFyFxKAFG0I2ARBNcqzr0df\/hmQUY7h7irqPZpRQSERREnTbWcQdRa+kHcN68eWHSLuuOHj0a7khiVJns69GapqamMCF13759rAp3LjEhmFunEXXuZgor9SIH5IAcyJgDQIpPfgVQEJASAwpwMvCBL4e7abiDBsXpGwcTUjbMG2GyahwlITqCBCUZG3wdTkUHMgswHPFTTz1FYUziRURVFi9eHNqIsmzZssVmzJgR1tO4ceNGiiAm+wI769atC3Notm\/fHm7FZiW3ThNeW7FihSHqtLFOqjkHdMByIDcOJEElhpSGtlt\/+7Xd7psEFFI5REwQYILi9I3AJDeniTrS5kCmAQZImT9\/vs2cOTPI4aXt2MMzXGbNmhXWsR3b+zpKtvd9ARraXNzh5Ouoe7tKOSAH5EBvOJCEFW4d9lRPHEmJoygOKZ7WiQGlN45ZnyEHsuRApgEmS0ZVPBatkANyQA50wgGAJZ6jEsNK34aRxkPTgBXu1CGKgpJRlE58jDaRA3XjgACmboZaHZUDcqA3HQBYmLfCg9pIBQEszFMhBeTpH4eVOKLCPJTePE59lhyolgO3+rkCmFt1UPvLATkgB9occGhhki3AwryVwsSHw90+wIqngDyyIlhpM06FHLgJBwQwN2GadpEDckAO4ADAgsqRlke\/G75nhyhLcdez4cFtngoSrOBYlqRjqXUHBDC1PoI6fjkgB7rdAaAEASaFiU2GiKog0kFMuI1vYfZIC89U4XZloKW467nw4LZuPzi9oRyQA8EBAUywQS9yQA7UkwPAid3zUHgEfr85Ky2GkiSYMG\/F7wpisi3zV3zCrU+6jaGls7cr15Pf6qsc6AkHBDA94areUw7Igcw4AKx4JMUjKMxPAUp4VP6VwtAwRwUoKZbSPg4lzFkBTHzeClEVn2zrc1i4jZnUkKAlM8OtA6kjBwQwdTTY6qoc+LUD+a0BLH7nj0dTSPEALERQAA7gJPlEWqCkWEr7OJSwncAkv+eJelb7Dghgan8M1QM5IAdKDpTB5dHvWqHtzh+PphBJ8QgKkAKclHbRfzkgB2rYAQFMDQ9eLR+6jl0OdKcDwAtzVQAXJtICK0RUPJqiSEp3uq33kgPZcEAAk41x0FHIATlwCw7wDcvcugy4EGG5hbfSrnJADtSIA3UKMDUyOjpMOSAHOuVA65vN4VH8zH3p1A7aSA7IgZp3QABT80OoDsgBOUCKiDuEuNuIibs8pwWYQXJHDsiBbnQgQ28lgMnQYOhQ5IAcuDUHgBjuMGIeDDDDnUcCSB8OTAAAEABJREFUmlvzVHvLgaw6IIDJ6sjouOSAHLgpB7jDiHkwwAx3H1UCmqsPLLf37v689btr0k19jnaqigP6UDlQdkAAU7ZCFTkgB\/LmAKmlSkAzcHij\/WLM580efsY8SsMTeUk79btrct6sUH\/kQO4cEMDkbkjVITkgByo5EAPNxy991e7dudyuvvjF8KWLpJ3Yj7QTD75zqAFoUIAaNpDkgBzIhAMCmEwMgw5CDsiBajkQQw1pJ746gNRTcdezBtQALg41TA5GfCUBUFOY2FRKQSlaU62x0+fWtwMCmPoef\/W+thzQ0faSA0BNcddzhmKoAWgQX0kA1PDwvDhaoxRULw2QPkYOlBzILMBMnz7dXn31Vdu9e3c70ca60rHbggULbOfOnWH9tm3brLGxkeayNm3aFNbxHsuWLSu3U1mzZk15HXXaJDkgB+RAJQcqQQ3RGr6yALBhX8Amhpo4YkPUhogOYltJDsiBm3cgswCzd+9emzNnjs2cOTNo0aJFdvr0aTt06JCxDlhZuHCh7dmzx1iHBUuWLKEIAlhGjRply5cvt82bN9vcuXMtBp9JkybZ2rVrg6gDQ2FHvVR2QGvkgBxo5wBQg\/jKgmS0xqHmysmDRsSmMPFhc7iJAcehpt0ba0EOyIEOHcgswCSP3OFk48aNYdXUqVNtyJAhdvDgQTt+\/Ljt2LHDxo0bV47CTJkyJQAPsNPc3GwXLlywadOmhX0nT54clvfv329bt24N29EWVupFDsgBOXCLDsRQc6n5G0Yaiq858Pk13NpdbJtj41DDpGFExAZ5OgrAKUxsMr7v6RYPS7vLgVw5UBMAQ+Rk0qRJAVKAFUZg+PDhFHby5MlQ8gLQDBs2LEDMoEGD7OzZszSXNXr06FAfMWKEXbp0KYBPaCi90FYqKv6\/\/\/77jUgNGjNmjDU0NORWgwcPtttuuy0AYp77meyb+j0kt+d0cqxZrtZ4D7aiDWw5av3e3Bp0ZcND4U4oe\/nJ0vIWG3DydRv44dEwFnf87r8KURvm2gxdut0AG4CGFBS3gdOPrqpa\/e7qcXb39vXUb65RXKvQyJEjK17Xan1FTQDMvHnzQsSESEpsOFGV8+fPh6Zz587Z5cuXQ91fTp06FapAD8ASFtpeYriJ622rrysef\/xx27BhQ9Bjjz1mEyZMyK3Gjx9vY8eOtfvuuy+3fUwbP\/Vb4512XvRG2\/g7P233DrlsYz\/YHTTmxDa780ffCeJWbzR27zdt2Ac\/CkBDCmrgl0oAtPQ1o0REbO6Y\/7SN+MKf2r1THqz4s6vzPP\/nOdcov16tWrXKzCyX\/zIPMI2NjSE1dPTo0XYRk94ejdWrV9vSpUuDXnjhBTt8+HBu9dZbb9mxY8fsyJEjue1j2vip3xrvtPMiK21vH3jdzr7yLSNiU\/zeAiMNxTyblv0vWVBLi7WMmBYezvfTn\/604s+uzvP8n+dco\/x69fzzz\/f25bLXPi\/zAENKaMCAAWGuS9IVTxnRTkqpf\/\/+VMvylBEQREqpvKJUiVNGcb20KvX\/mTNnwgRiJhGfOHHCWvhlkVMR2froo4\/s4sWLue5ncgzVb4138pzI6nLx3HHjCcM+eZiUUr+7J9vVljPG3JoPfvaTij+71TjPs+BjPfWbaxTXKvTGG2+kXtPy0Jh5gGHiLWASz3XB+LSUEScoKaW0lBH7eEopLWWU1sY+khyQA3Igqw4wsRcBLxxj34aRVpj4sCb8YoaUewcyDzBEUbh9mruJ4tHgDqLW1lZjfgwRltmzZ1ucZjpw4ECYdMsE4KampjAhdd++feEtuHOJW6y5dRpRpy2s1IsckAM5c6A2uwOYoMLEphKUNBmTd3kCMHNdEHcsMbEX0UOiL+jah6dZlORA7h3IPMBUSu8QZdmyZYvNmDHDeM4LI+W3WFNfv359uD163bp14Tkx27dvD8+PYR23ThNaW7FihSHqtLFOkgNyQA70pANACSLt43ACmCDAhDuNYjjhDiRutS6UIitEWIAUxPwXxFwYbs\/mNm3ELdvcxt2TfdB7y4EsOJB5gFm8eLGhNLOAjlmzZoUH3c2fP\/+6Sb7s5w\/CA2ji91i5cmXYj\/XU43Wqy4HudEDvVR8OJKHEIyaVwIQ7iRxOSAE5nPBEXwcT4IQn\/QImCDjhuTKI+S+IuTD14bB6KQfaO5B5gGl\/uFqSA3JADlTHAQAFeeTEwYSoiUdMHEoAEyImMZjwRF7ABDmYJCMngElx13PmYAKcKJpSnfHWp2bfAQFM9sfoFo9Qu8sBOXAjB4AS5GBSmNh0w\/kmQAopHeCE9\/V0DlCCiJikgQlwApgggQnOSXLg1hwQwNyaf9pbDsiBDDoAkKAYSgAToiaIqEmfR39gb89aZ31LJRNhARMiJ6hwg\/kmAArpHOQpHYcSgUkGTwYdUm4d6HGAya1z6pgckAM97gAQggARVChFR1ByfglA4pNfPZ2ThBKPmvhck77nDtuItzfblVf\/IjwULhk9cThR5KTHh1kfIAduygEBzE3Zpp3kgBzoyAHAwwV8IODDFUOIR0VuBCLJCAkpHGCE4yCNwxyT4q5njTkmyIGEiElaSufKq2us4ew+s3deCw+FU\/QEJ6UMOaBD6cABAUwHBmm1HKhnB2IAcfCgjOED6EBEQBAREEQExAV8INIziGhIoZSmcQgBRIAQ5HfhACEIEEEOIkkYiSMlxbYJsJ7SYQIsqucxVN\/lQF4dEMDkdWTVLzkQORCDSBwJAUSuPrDczt67yAqPfNscRAAQFAMI4IEcPoAOBHQgIiA3go80APF5JJUgxEFE0ZFoMHurqs+RAxl3QACT8QHS4cmBNAdiICl0Yl5IDCIeCXEQGXzP5+zKgKHG9+sAIjGEEPlAwAe6UfSDuSJxBCQJH4qEpI2k2uSAHLhZBwQwN+uc9pMDPeBAEkyIkDA\/BMXRkRhIiIqgQpSSAURQJRiJQYRvNr7zR9+xvrvWWRJCiHwg4AP1QJcrvaXa5YAckAM3dEAAc0N7tFIO9IwDgIqncmI4SYJJEkqSQOJRkRhIKqVjABEkEOmZMdW7ygE50LsOCGB61299Wq040M3HWSileZKgQiqHNA4TWYmWMIfEJ606mDBHBDmUxCkawUg3D5LeTg7IgZpyQABTU8Olg61FB4i0pIEKkAKcIFI3yOeNKEpSiyOtY5YDcqA3HRDA9Kbbnf8sbZkjB4iUkPrhjh2iLddaTofe9W0YFUq9yAE5IAfkQNcdEMB03TPtIQe67ACpHyIugAzRGCbdkkLiVmXE81OYpMukXUTKichNn4aRXf4s7SAH5IAcqAcH0gGmHnquPsqBXnbgWssZA2SYcIsAGm5RLu561gAbDgdoiQGHSb0x4DCPBsBBDjnsI9DBPUkOyIF6ckAAU0+jrb5mygGAhvRScddzAWyYqIscbhxwmNgL4DDRN6Sh7ppsMeQQyYlBx6M5RHSSwOPQI+DJ1Kmgg8mRA+pK7zkggOk9r\/VJcqDTDgA3CMBhYi+QwyRfAAcBOQjIQURyHHQcdvgwgKcw8eEy8HjqyoGH6E6fR39gx6Z\/3SgdeiqBj6I9uCrJATmQBQcEMFkYBR2DHLhJB4AcFINOEna4ywnYQcAOAngQ0GNvv2pDzu4zbuP2KE8l8ImjPcCPR3tIaSEAB91kd7TbLTugN5AD9eOAAKZ+xlo9lQMG7CCABxHdubb\/r+2Okz8sP4mXCA9KAx+HH8AHebQHaCGtBeAg4CYWoBOLSE8sIj4uQCipwsQmc\/FZsUiHuTTEckAO1I8DmQeYNWvW2O7du4O2bdtmjY2N5dFZsGCB7dy5M3UdG23atCmsY\/9ly5bRVFb8vtTLK1SRA3LgOgeAHpeDD\/ATR3uI8AA4CLhxFdsmKQM7LiI9Lj6MiI+rUEp5uYAiROrLBSDFIh3miqGJegxN1GNouvrAcjt77yLrN2dlKcX2xHUqVIAmwRIjJsmB6juQaYABOmbMmGFr1661mTNn2qVLl+ypp54KrgEyCxcutD179tiiRYtC25IlS0LJC\/uOGjXKli9fbps3b7a5c+fa9OnTWWWAz6RJk8L78t7UaQsr9SIH5MBNO+CQA9y4im2TlOOS+Twuoj2xiPy4gKJYwFFSpMJcDk2UxQQ4eYoMcKKDA4c3hi+xvFIYWoruPFwWwIQcmCjTgAlAcgFHsRyUPKpEGUeVHI6IJAFEHI8kB+RA1xzINMBMmTIlAMrWrVtDrxYvXmyIhalTp9qQIUPs4MGDdvz4cduxY4eNGzeuHKFh39OnT9vevXutubnZLly4YNOmTWNXmzx5cljev3+\/8d5sR1tYqZcadkCHnncHHJDikoiQy6GJspgAp0vN3zAXwPTxS1+1O3\/0HbPmlebARBkDE\/VKwAQkuTyyRBmDkkeVKAttkaUYjgAjokeAkAMQ8OPAY\/onB+RARQcyCzBES4ignDp1KvXghw8fHtpPnjwZSl4AmmHDhgWIGTRokJ09e5bmskaPHh3qI0aMCNEcwCc0lF5oKxUV\/99\/\/\/1GpAaNGTPGGhoacqvBgwfbbbfdFgAxz\/1M9k39HpLbczo51ix3drwHW9FiDWw5amjAyV3m6vfmVnP13bUuzCeiBI5c1178oqHi9xYYImoEALW+2RwmUPdtGGmIpzUPKKWvAB2Os7vV2X539+dW+\/3qqd9co7hWoZEj8\/swzMwCDCRx+fJlu+OOOyrOcyGqcv78eTa1c+fOGduHhbYXhx9AhfRTW3MoYriJ62Flysvjjz9uGzZsCHrsscdswoQJ1ykvbePHj7exY8fafffdl9s+po2V+q3xTjsvaLt3yoPmGvt7v2+uOx\/6Y3ON+MKfGrpj\/tPmIpqCBn6pBDhLX7N+JVFHRF9ITwErQEu\/1vcNNZzdZyPe3mz37lzeIz9\/Os\/zf55zjfLr1apVq1KuaPloyjTA9O\/fP\/wA\/9Ef\/VF5nstTTz1VFedXr15tS5cuDXrhhRfs8OHDudVbb71lx44dsyNHjuS2j2njp37X\/ngfefcXFuvtC\/0N\/fQ3xpjr2KdmGvrnf\/F5+6f+M+ynv\/WwnRgz39797FeCeCaO6+1Z6wz5MqVvR8kkYPTzEf+nnf\/UZ4Pe\/8077b2P+wa9\/7OfGGrZ\/5IRbXERfUFEYq5seMgQdfTetq\/Zu6\/9VY\/97Ok8r\/3zPO33V9zGNcqvV88\/\/3xVrpld\/NCb2jzTAEOPmNtCBAVRJ61Eeol1njKiTkoJ4KHu8pRRY2OjDSqllLydMk4ZxXXWpenMmTN26NChoBMnTlhLS0tuRWTro48+sosXL+a2j2njp35Xb7wv2kBzFRvGGWq9+wFDVyYuMMRdQ8ia1piLh+8hIhuo76M\/sFj28DOGrj243Fytdz9Yet8H7drw37EP\/2fB3r\/czy6VPh\/QQBff+ScDOJADByXA4fJ5McyRQcydicUcG+RzbpiPw7wcl8\/ZKZ473us\/YzrPq3eep\/3e6Yk2rlF+vXrjjTfSLp1zBvoAABAASURBVGm5aMsswJAa4gct6TJtrEtLGfk6YCeZMuJ9PKWUljJKa2MfSQ7IgZtzgLtrEHfaFCZ+8hwXJqci0ipMWmXyKmIiq4vUCiK9ggptk19JszBHhKPhTiLEpFkEYLgcMigdNCgBDeSgQbTDJ\/HGsAF0ABzIgYPSoYPSJxFzLFIOHFAXatKBzAIMEHL06FGbPXt2mJRLFIU6bazjDqLW1labN29eWB+vYyQOHDgQJt0SrWlqagoTUvft28eqcOcSkRxunUbUuZsprNSLHJADXXYgBhUABTiJgQQQYVKqQwjwwd06Dh2UAAcCNhCwgRw4KAENBGQgIAMBGC4Aw+WgQdnlTmkHOSAHMu1AZgEG11auXGkAC89xQURVaGMdELNlyxbjOTGso23jxo0UQevXrzduj163bl2YP7N9+\/ZwSzUruXWa8NqKFSsMUaeNdZIckAOddwBw4fkmwApREyCFKInDiQNJDCIxgDh0UCaho\/NHUfNbqgNyQA7chAOZBhj6A7DwEDvkz4ChHQEds2bNCg+5mz9\/fngeDO0utmc\/BNB4O2X8vtRpk+SAHLjeASAFeSqICAspIE\/9kOIp7nrWiJp4lIToiEOJoh\/Xe6oWOSAHbt2BzAPMrXdR7yAH5AAA4gJEUKFtXord85C9d\/fnwyP1HU4AFNJADikeYSEVRBqIFJCnfYCW4q7nwvcsyWk5IAfkQG85IIDpLaf1OXKgEw44ZHgJaLgcOChJ27iADgR0uIAPBIAgAMRFqgcBI6jfnD+zCyOm2dXhE8KD1IATFE+OTaaC4ghLJ7qlTeSAHJAD3e6AAKbbLdUbdtGBzG\/uMOGlAwUlMBHLoYISqIjlcEEJXLgADJdDhpeAhgvYcBXa7syhJCLCvBOMBDwQc1BiAPFoCSCCSPcg5qbwDJKxe78ZnhIbz08hqkIaCDE\/RakgHJbkgBzIigMCmKyMhI7jphxwqPASqHAV2lIklAAFqgQUQAXPE+FBZZQOFJQOE146UFA6UHhZSAEL4AIBFi4HjCRkVAINhw2AA5G2iQV4ICIjrhhAHEIAEQSMoJsyXTvJATkgBzLggAAmA4OgQzCLAaTQBh4xcAAYKBm1cKjwEqhwOVRQFtrAwqMVAAW+O1BQ9j132Iac3Wf29qvtnpqahIo4gtFZsEjCBZABYCDgIhaA4QIyXByvJAfkgByQA584IID5xAe9VskBICWOcgAfAAcCOhw4AAxExAIBFS6AIhZQgYhUuNKiFUmouPLqGrvj5A\/t2v6\/thgoqDtQeOlQQVkl6\/SxckAOyIGqOlDtDxfAVHsE6vzzgRIU28AyuvbhaUOh3lKql8R2LFNe+\/BMaf0Zu1pqR0m4YBtJDsgBOSAH8umAACaf41ozvSK6QXSEiAkRleKuZ0vRj5eNOSJXW86EfpDuIRqDeFAa0RlEtAZ5+ohIDvI0EyknxLwXIj2uQilF5fNkKD19dWXAUEO+HD5cL3JADmTUAR1WvTsggKn3MyAj\/ScVA8wUdz1niDkiiDQPAnKQp4QAHpenj2IAIs2UBCDgB8XwEwPQwC9tNSbx8mWAQBEwVEkOSWkl0FRJwFSaHK4qlUBXmgCwSnIQi8uMDLcOQw7IATlwyw4IYG7ZQr1BNRwAeFyeOooByCEI+EEdwQ8QZC8\/aXy5X+tLX71uEi9w5Cq2RYmApKQ8cgQ8xXKPiCalqdA2yThZAlwI6EoTAFZJQFhSaUDGtzi\/PWudcfdVV4AsDcRoqwRhtBdK0a+kkgAWAxd19y5rpY5HDsiB6joggKmu\/\/r0Kjjg4EPp8PNJecgGfnjU+p77ZwOGKqnYFiVKK4kapQmIupEArDR5xCmt9AhUWgmQpckhLC6v7vtr+\/SJH4a7r7oCZGkgRlshBcaAMNQZCOsIujoDWYBULOApqT5T\/zA8gbjwu01WSAEr2pJwlVwGsG6kKpze+kg5UDcOCGDqZqjV0Xw5YAaAVdInQHbQkmUalHHXld99lYSyNBijraswlgZgtCXhKwldMWhRL6ZEv5JRL8YZkIpVSEAVMNV32h\/aL8Z83q49uNzSwIq2StEtb0\/CVnI5LeIVt6XBWNxWKRXp7TGkxXU8kORA3h0QwOR9hNU\/OVBDDlxtOd3uaP1OM2+82sJdZ59M7va2a4l9Urdpu6PtWlsJ9KB+re+bb58s\/f17sowhK63OYwRi9bl9lMWK18X1njxmvbccyIoDApisjESNHYcON\/sOVEptJNMgV4f\/jrWMmGZpqZRk2oXl+C\/9uO5RgbiMownU4+gD9WTEgmWPbsQl0ZBYyYgKF+8kADBCSSiJl\/t9\/H5IGQIyydQZy0R7XER\/OlIyepS2nIw4JZeJSt1IaWnGtDb6LsmBvDsggMn7CKt\/3eZAZ4HAAaGQMq8CAEgqhoC4HoOA14GApACBNAEDaYrBgPqAR75tZ+9dlJpKKSRSLywDC6gzwJCEA4AgDQSSF\/uOLuxpF+20tBbprkriwYUj3t5sfXetC3e+FW8wtykt9ZZsS6br0pavlSJIN1K3nax6IzlQBw7UKMDUwcioi9c5kAQIBwUvCwlg6AgUHAq85C4cbqOmBBKSUJAGA7QBAWmKIwZeZ+5FIQEFaTAAHGBAHDGgngQCjxSkQQFtSTDw5RgQrr74RePLHPlSx+Rf\/2mg4G2dBYYkGCQv\/CwnL\/bJizxeSHJADsiB2AEBTOyG6j3uQBpsABppkQcgAjlIAAuxktDgkOBlR6BAZ4ECV9\/ou5AAAwAgTQ4ByTKGAupJGPBlB4C4TIMB2tKiB0kgYBkISFMSDHw5CQjMBcEPSQ7IgZw7kKPuCWByNJhZ70qhFCFx6HDI8IhEHIVwoPBoQwwRN4IGBwQvY0DwOlAQKwYEUgrx3ThpQECbQ0CyTEJB1sdDxycH5IAcqGUHBDC1PHo1duzASqVD5u4QB5drLacNAQi0XWv7ziPKq6V1LtZfi+YUVHpvtcsBOZAZB3QgcqDbHMg0wEyfPt1effVV2717d1lr1qwpd37BggW2c+fOsG7btm3W2NhYXkdl06ZNYR37L1u2jKayeB\/aEfXyClV6zIGWF\/\/EPIJCVKUYPdPjaglEmPeBChMfNmDHozQetaGMU0ieWqIk1RTL57XEqSlSVbEKpYiQi9QWd+MUbx9n\/e6aZMn5Nj1mit5YDsgBOSAHbsqBTAPM3XffbZcvX7bly5fbzJkzg1auXBk6CqwsXLjQ9uzZY4sWLQptS5YsCSUvAMuoUaPCvps3b7a5c+caQMQ6wGfSpEm2du3aIOq0sU7qOQeIlhA1QaRiitFdH6RyPLVDusfTQMwlieUARAkEoWIEQsxdIfUEECGAyFUogZErBiRACTjibpx3P\/sVs4efsRiUqANJsWJYou7ARBlDE\/UYmrxeSMBTv7sml8Dp14oBqudGpA7fWV2WA3IgNw5kGmCGDx9ura2tdv78+esMnzp1qg0ZMsQOHjxox48ftx07dti4cePKUZgpU6bY6dOnbe\/evdbc3GwXLlywadOmhfeZPHlyWN6\/f79t3bo1bEdbWKmXTDkA9MQCflxAECpGIET9UvM3zOVQRAkYuRyQvASS\/G4cSgDJBSTFKt4AmJLQBDwVEuCUBk8AVCygyRWDk9eBpqSAp1jAUyyHp2Rp9zzU7jkwSZhiOQYq6pk6SXQwckAO1KUDmQaY0aNHVxwU4IaVJ0+epAgCaIYNGxYgZtCgQXb27NnQ7i\/+fiNGjLBLly4F8PF1tHk9rbz\/\/vuNSA0aM2aMNTQ05FaDBw+22267LQBinvsZ922wFe22a5fs9r4f26BffWQDW46WNeDkLovV782tFovniMSy5lKUMNK1F79obTJuVU4KYIplLz9psfq8vs5i9XtzS+l4Xm8nvsNpUKkPLvo29K7PmGvwPZ+zhqmPtBMghfrN+bN2z4GJQcrrDlNeOkjFZRKofDmGKq\/HYOX1JFjFy4UoYuV1wCpNA4c3WlL4kVQ9nud4oH4Pye3vbsYXcY3iWoVGjhyZdknLRVumAQaoIA1ECmj37t2WnOdCVMWjM+fOnQvppnhUTp06FRaJ0AAsYaHtJYabuN62+rri8ccftw0bNgQ99thjNmHChNxq\/PjxNnbsWLvvvvty28e08atWv8ff+WmLde+Qyxbrt391wmKN\/WC3JTXmxDaLdeePvmOxeM5LUvfuXG5o0oH\/aI+8\/z2jTG7Dcvw+XucBcLE+feKHNuyDH12nob981+4oXG0nfsE6WMVlErBYBrAQab6kHK6S5cAvbbWk+i19zZJq\/cO\/sZ9M+n\/sf\/7fz7fbHsi6Y\/7TNuILfxp050N\/bGN\/7\/eD0s6bWmur1nlebZ\/qqd9co\/x6tWrVquuuZ3lpyCzANDY22qBSFOWdd94Jc198nstTTz1VFe9Xr15tS5cuDXrhhRfs8OHDtadOHvNbb71lx44dsyNHjuS2j2njp34fsbcPvH6djv3j31pS7772Vxbr7CvfsjS9t+1rlhTpvDQVv7fAkkpGq3w5uZ0ve9ovWZICJPXnan2z2ZgrdeXy5fD7pPXjj+3KgKFBxdvHhZTaL8Z8PnzZI08pZm4UenvWOjvyv6+s6Z8Lnef5\/73GNcqvV88\/\/3w4x\/P4klmAIWoyf\/58W7x4cfCdZea5EJHxybieMmIDUkr9+\/enWpanjByGyitKFaI7pSL8j+uhIeXlzJkzdujQoaATJ05YS0tLbkVk66OPPrKLFy\/mto9p46d+d+94X7SBFqvYMM6Sar37AUvqysQFFuvqA8stqX4PrbRYhUe+bWjQw39urjgyQ\/SGSI5rQCklxbOHHFr6Noxs96WOAW5OHgyQA+y4ACAmiqedP7XSpvO8e8\/zLI471yi\/Xr3xxhspV7R8NPUmwHSLY\/zwkTZKSxn5OmAnmTLiwz2llJYySmtjH0kO5MUBJt+6fO4I80mSk3jjuSc+P8VL0iuxfJ4LZTwfhrrPl\/EyBgqvAxZJFaJJz9QBDQAjFpOlYwEcgEUsoi5JeWQmnrRNVIfJ3D7BmzKOEPmEcMpi24TxvJwT6occqGUHMgswRFl4Bgy3Q2Mwy9wKffTo0TD5ljuIuENp3rx5YdLu7Nmzzdex\/YEDB8KkW\/ZramoKE1L37dvHqnDnEpEcbp1G1LmbKazUixyoQQcAkkIpqgB8JGHD4cJBgjIGiHgSLxGKQgQQSXiIoYF6DA7FXc9aJWCIwQF4iAU8xAIgkoqBgjowkVSxDS685A61pPwONr+zTV+hUIMne10esjqd5kBmAYbbn5l8BLQwgXfdunXhdmd\/DgxRli1bttiMGTOMSb50buPGjRRB69evD9uzH\/Nntm\/fHm6pZiW3ThNeW7FihSHqtLFOkgO16EC\/uycZkQzgo9\/dk41oBYCBiEo4WKSBBBEIJvNSJsGBZYDBlYQGlh0YKCsBQxIcHCAoa9FvHbMckAPVdyCzAIM1QMycOXPCJF4eZOfzYViHgI5Zs2aF9cyXAWrwwU43AAAQAElEQVRod7E9+yGAxtspASHaEXXaJDlQSw7E6SCOG1i59PI3zKEDuEAxWKSBBPtKcuBGDmidHMiiA5kGmCwapmOSA93pgEMIJWkgVCilglykhFCcFkqmhEgHkfrhuAAUSkkOyAE5kHcHBDB5H2H175YcACxiARguhwwvAQ3ksEGZNuGVCa4u5qO4ABFEKshVaJuPEqeFmHcSp4SYT8IcEiIvt9TZzO6sA5MDckAOXO+AAOZ6T9SSMQdigPC6Q4SXhShqQR2QiAVMxIrBwut9Hv2BHZv+9fDAszTAADQADJdDBiURkEICNpiHgpWkdhDgwXyUYmKya3JeCjDiAkpczEMhJYR8rgkRF80jwWVJDsiBenNAAFNvI97J\/vJ8DIcFLx0WvCwkoIHlGBqox9DgdQcGL0mJxHJ48BJwSMohwksgIlahDSa89AgGUIGwAaiIZW+\/akPO7rOr+\/663d00DhheEvFADhmUDhleAhsugMNVbLtTxgGEEghBgAji2CQ5IAfkgBy4sQMCmBv7U9W1Dg5eOjh4WUgABMAQy4HBSwcGLytBA49cJxLRtxSRiMHBYcHLGBi8XugAHICHGBqoe2SC6AQqJiIUpEscHuISiHABEbEcJOLSgcJLhwovr+3\/a7vj5A+NErBwARexgAxU1ZNDHy4H5IAcqHMHBDC9cAKkAQegkQYWDhVEHz6Bh+3mpYODlw4NXnYED8BCrBgcYmjgiwP5npvWl75qMTBQd2DwMoYG6jEweN2BwUsHhrgstkUmvHR48DIGCK8DEa5eGEZ9hByQA3JADmTIAQFMLwxGGnAwZ4LHmSPSG31uHxWOBKhAxbYoBNAQy8HBS6DB5cAQlw4NlDEweL0YgYPDwsc\/braGUiql77l\/NocFLx0YvAwHrRc5IAfkgByQA73sQIcA08vHk8uPA0CKJSBBwAlRkGRHSa0AMgANAnCIrDj8UPIdLw2Pftfa6Q82mqeEPKLjJVGepAqJtBPLHiHy8urw3zG+0K7fXZPM01dxmTx2LcsBOSAH5IAc6G0HBDC94DjRi2JbpINICBESoiYeRaEEcmIx7yNWsQRAzA+J5TAEECEgKFYhmo8CECGgKKmGR0sQFGnAI982vnnXHn6mnL7yNBYl6a1K8hRYWumglSwduJJlEr6Sy4UUGIvb+t012TpSDGZMXEZxG\/VeOEX0EXJADuTPAfWohx0QwPSwwTd6e0\/DUAI5scrpnDebjXqxDYDi0tNAXgJHsQAlF8CUJuApqasvftHG7v2mpc2BAbJisIrrxRTIcuBKwhbA5d7E0OV1olGFCMC8DoS5kiCWXE6CWdoyQOYa+KWt4Tbq5OTlSrAWt6cBW7ItCW5py0mISy4nIc6XCwmYqwRuAFksHwOVckAOyIFac0AAU2sj1s3HCzyliS+565syBwbIAqjSVEyBLG9zyIrLGLaSdQevZJkGYcm2JJBVWgbGYtnLT9qdP\/pOGdxiOKtUL5agDTmoVSrTAA6Ic\/mwOsCllWlQVwnm0mCNNoc1L4Ew7jp7e9a6ds+\/ScKXLyehC4DyY1eZMQd0OHIg5w4IYHI+wPXYvTQgS2sDxtrrkA388Kg5uKVBWrKteANoi9fF4JZWTwJc2nIS5lhOwhvLMZR5HQArtsEWZWspsocAqz4lUKXfMUylARRtQFSsPg2fTD6vx\/NMfZYDcqC6Dghgquu\/Pr3GHIjTL3E9LWWTTOv4MlGLWMk0kS8nox0seyQkLomixCLSkhTpNY\/WUDJRHAEj19ombftQOMhQAjixgB4XUASM+X6JUotyQA7IgR51QADTo\/bqzavhQBpYODxQOjw4KFACB9a0Jkxe5isFHBBiMKDuqZdkmQQGloGGNBUSc3uACNS3YaTFwjsgIhYwkUxTARJJeeQlLpOpNKI16MqGh+zencvt2otfNKI6sZKRIIDFRTSKY5TkgByQA9VwQABTDdf1mT3qQMOj37WhSz95AGBD291VMUgQgUAegQAe+tw+yvp9\/L4x9yeGhCQYsBxDQVxPAgLLAEJSMSDE9SQssOywEJfFRNoKkEiqfWrsYHieTzKN1qODoDeXA3JADvSwAwKYHjZYb9\/7DhQTcz0AEkQkI3k07SIe98yxlhHTDLABcBDgQxmL5\/F4FMfLwsQmi7eh\/sn7PFFqby+2jZWWfvK2OJqUrCf7omU5IAfkQD05IIDpmdHWu1bRAaIRxbYoBZELIhmIaIdHQ4iOoDiC4l+hcOXVvzAiLahYgqFkysZhCCBCdDUGIa8XEqkioAYBRbE8SpRWJlNV8TIprUryFFilkpSZy1NnhUe+XfGhiKTZXA5tlcpCCebS5FAWl4Iyzh5JDsiBm3FAAHMzrmmfmnfA0ylxqsW\/QsHeec2AIFRsA6G4BIpiAUdpAphiOTzFJRCVphis0urAVSUVU6ArhrAkgHnqrHjueHlcHcKSJem2QgLMfBk4QzGcxfXOAFoakKVBmMOXlw5XXqbBVblzqsgBOZBhBzp\/aAKYznulLeVAtzvgIJUsY7BKqwNXlVRMga64LYYv6ldeXWN8eWffXessDcTithjIkvUYzJL1zkBaGpAVU2AsCWAMioMWgIUAqaQAoUIpOsT2khyQA7XvQM0AzLJly+zVV1+16dOnl11fsGCB7dy503bv3m3btm2zxsbG8joqmzZtCutYz\/60udasWVNeR93bVcoBOdD9DiQBjeUkmCWBDFDhlu0rJw8ZqToUH5lDCxOwEeDibfF2cZ30VbysuhxIOqDl2nGgJgAGaJk7d247V4GVhQsX2p49e2zRokVh3ZIlS0LJC8AyatQoW758uW3evNnYn\/dhHeAzadIkW7t2bRB12lgnyQE50HsOMAcGqCAyQvqHlBCRElJJPt+H1JOnogql9JWDCkCDAB1SZB7BiVNuRH48GkTEiIhT7\/VOnyQH5EBPOlATADNv3jy7cuVKOx+mTp1qQ4YMsYMHD9rx48dtx44dNm7cuHIUZsqUKXb69Gnbu3evNTc324ULF2zatGnhPSZPnhyW9+\/fb1u3bg3b0RZW6kUOyIFucyAGlEIpfcPcFAcVYAVIAVBI9xA9AUjSYMRBBAhBpLUuNX\/DXMVS2swjOHFkh0hPt3WmV95IHyIH5EBnHcg8wBAZueuuu+z1119v16fhw4eH5ZMnT4aSF4Bm2LBhAWIGDRpkZ8+epbms0aNHh\/qIESPs0qVLAXxCQ+mFtlKh\/3Kg7hy4MmCoIWADERFxAR0u4MMFhLgAERfRE5dHUQAUIihASqEUQXFQAVaImhAxiaEkDUYEInV3WqrDcqBDBzIPMF\/4whfsn\/7pn+yjjz66rjNEVc6fPx\/az507Z5cvXw51fzl16lSoEqEBWMJC20sMN3G9bfV1xf3332+kmtCYMWOsoaEhtxo8eLDddtttIcKV534m+9bb\/R44vNGScnDwslCKWiTlEEHpEBGXDhNx6VDhJXDhir+Fe+jS7QZwuAAPV8PUR8w1+J7P2dC7PhM0yIrm4juVBpx83RC3pdvLTxriab886RdZ80pjwjAacHKXDWw5WpWfpd4e7+T5Vq1l9XtIVc633hxvrlFcq9DIkSOvu57lpSHTAMM8Foxev349RVX1+OOP24YNG4Iee+wxmzBhQm41fvx4Gzt2rN1333257WPa+Hm\/x93\/f9i9Ux60sb\/3+0F3PvTHFmvEF\/7UXHfMf9qSisGBOoDg6rf0NXN5W1w6OHjp8BCXDhGUgARymKDkF+UdhasWa+gv37VhH\/yoLO46co079f\/ZrAv\/1f7lkWftzh99J2js3m8a4isGXCzH8m3HnNhmscZ+sNvQb\/\/qhN075HJQmt\/VbvPxrtfzXP3O7+9wrlF+vVq1alVVr509+eGZBZjGxkb73Oc+Z6+88krF\/nvKiA1IKfXv359qWZ4y4r0GlVJK5RWlSpwyiuulVan\/V69ebUuXLg164YUX7PDhw7nVW2+9ZceOHbMjR470ch9719Of\/sYYQyfGzA\/fgfSTcY\/ZS0O\/ZD+e8B\/s2PSvh7Z3P\/sVO3vvonY6\/6nPmuv937zT0Hsf9zXX+z\/7icVq2f+SuUiZxCJ9klTxewssFtGLWPE6rzMnJNZ7275mSZ195Vvmeve1vzLXz\/9hk11853\/Yzw\/8Nzv2j38b9PaB1w3pPO\/dc7I3\/K6Xn++kl\/XUb65Rfr16\/vnnU69peWjMLMAwSZe7iFasWGHcBs2dRkDIM888Y8yLSUsZeUopLWXEYHlKKS1llNbGPq4zZ87YoUOHgk6cOGEtLS25FT6Ssrt48WJu+8j4fXz3A3btweVGWbx9nCEf72TJnS7ok9t6DxoTTWO1\/rjZXFd+xvrm0jafiNuAke\/L+5RV2pbtY\/FAuVgca0+qXsY76aH6ne+f73oeb65Rfr164403kr\/OcrOcWYDh7qBZs2bZzJkzg7gVmnksTz75ZLhziDuIWltbjTuUiLDMnj3bjh49Wp6Ye+DAgTBnhVunm5qawnyOffv2hYHjziXgCBBC1GkLK\/VSNw4QCeHuFkQUhLkapFWuvPoXVtz1bDsx4RT1bRhpqDDxYUNMTEVxioe6p4DikvklSfk8lLj0eSrJknRUrHjei9eZF5OmQmIujc+xobw6\/HcCvDGBN6m6ORnUUTkgB2rOgcwCTEdOEmXZsmWLzZgxw4Abtt+4cSNFEPNmuI163bp14Tkx27dvD7dUsxI4gk6J7iDqtLFOqh8HuLPF9UkE5JA1nC1B7juvleDluXby23U9TcNdM8ifMRKXAFGagKRYAFRSgFMc2fE6ERsAysUoAVJJFdrAykvgCgFVsWKwGvDIt0O6rO+jP7CbBawYrLzuUEUZQ1WhDaaAJ+TQRJ8kOSAHatKBqhx0zQAMQDJnzpwyhOAW0OFRmvnz55ejL6xDixcvDtEbojjsT5tr5cqV5XXUvV2lHLhVBxyKkuUnkHTQvPz4zU9STHFZ3PVcO3DyZQcoLx2kkiVQFSsGK68n4erqi18ME3ZbX\/qq3SxgOVh5iYcxXAEqhRJcxTDlEOXQ5FGoOPIE\/PBekhyQA3Ig6UDNAEzywLUsB+RA1xyoFOngGTBoQKFQfsM+t19\/6yVAUt6gVOnTMMpcpcXy\/76lNFt5oa3S5\/ZRbbUbFwCQb3Htw9MhXefLKuVARQe0oi4dEMDU5bCr0z3tgMNCXBKFSKrQlk6JyzjdQgQiKU\/RxGUctaDu0Yy49EhHsvQUkj38TPkZMHG6yeuFUgQlFo\/0RwBLLLwFRGJ5KsxTZnGkh4iQR4fi6JFHl3g\/SQ7IATmQdEAAk3REy3Ig4YBDCPBRaAMOIMPBApBIQkMSElj2lElcOhzEZSECBQABdQQIzJFxSKAs7nrWHBaSZQwP1Jm8zPNdSCUBE8iBIi5juIjrDhpx6akuL4ttqTFPl3kajZJUW8LyWlvU8coBOVAFBwQwVTBdH1lbDgwoQYsDiIMGczloRw4YRBwACcQt08UERAALRxJJ6gAAEABJREFUwEGsGBC8HsOB12M4oO5gEJfFNkjw0mEhWQIN7XXI+rW+b4CEq7ZGSEcrB+RAPToggKnHUVefu+QAMEIUo1gCEhdtCFhBwAtREmAGATZAjgMPJZEXQAg1PPrdUrrmu0bkBhHFcXlkh5JIT6xCCaZcRISC7ppsHiXqUse0sRyQA3Kghh0QwNTw4OnQe8cBohJEMYpRhMMjH0RDEJESj6DEERaiLi4gyEWaJxYA5AKEXIUonZQGRA2PbiyB0Mby7c8+5wUocjkYUQJFKIYiu+eh8ByY3nFTnyIH5IAc6B4HBDDd42M9v4v6nnAA4HHFqRogyBXDEHUHIkqAyAUYuRyQKGNIop6EpDQ44jCTYAQU9ZvzZ+E5MHxHk0MPkFMoRXvYR5IDckAOZNEBAUwWR0XHJAc6cMAByUsHpTRAAoqQQxFlDEV8zxKTeHkOjIMPoBPgppSe6uBQtFoOyAE5UBUHah9gqmKbPlQO5MsBJvH2PffP5YfoATktL\/6JDXr4z83TTcy3yVev1Rs5IAdq2QEBTC2Pno5dDvSgA0R3Lr38DQNciMYw3yZOMTnYkGpiG9SDh6O3lgNyoJsdqPW3E8DU+gjq+OVADzpAaopoDPNsUJxiAlgAG7\/DCsBJTiJOmzgM8CD2R9xB1YNd0FvLATmQUwcEMDkdWHVLDnSnA0RjEBOOEVCDmFCMgBvEZGLutHLQ4c4q5tOgwsSHS+moLxvAgwAexG3laeDj8BNHegoTm0JESODTnaNbrffS58qBW3NAAHNr\/mlvOSAHSg4AN4iIDROJgRwUTx5m4jCwg4Adl0NPcdezFoNP6W3DdyEVbgA+pLQAHSBHkRwck+RA\/TgggKmfsVZP5UBmHAB2XEnoKYbn7fz6qxCYh+OQ4w8PJLJDZ4js8OBAT2XR1llpOzkgB2rbAQFMbY+fjl4OZMYBIiAuUjyFUrrHRYQE8XwZF5ETF5EUTyNRklZykWZCpJ2AlT63j7JrH542YAaRskKAT2bM0IHIATnQ4w4IYHrcYn2AHEhzoHbagBJAJAaQGDwcPhw4KB04gA5UKKWBEABC1AQRRXGROiqWUkiAiIuoC\/JUE6kn0lCI+TekpxApK0Qkp3Zc1ZHKATlwqw4IYG7VQe0vB3LqAFEUQAUgIUUTA4iDB98DBXzcCDpi8HD4iAEECCmGtNFzBoi4ABLkqaac2qxuyQE5cJMOCGBu0rha303HLwc64wARE7YjYkIJsAAVQAbggYAPlhHrkKADtyQ5IAd60gEBTE+6q\/eWAzXsACBC9IQUDhEWIi2ADNEYUkTMVSF9RGqpMLGphnuqQ5cDcqAWHagSwHTOqunTp9urr75qu3fvDlqzZk27HRcsWGA7d+4M67Zt22aNjY3t1m\/atCmsY\/9ly5a1W8d70Y6ot1upBTkgB8oOEE0hukKkhdSPQw3zU4AaUk3McwFmmKAL0NBWfgNV5IAckAM94EBmAQYYefLJJ+3QoUM2c+ZMW7t2rc2YMcMcRFi\/cOFC27Nnjy1atChYs2TJklDywnajRo2y5cuX2+bNm23u3LkGELEO8Jk0aVJ4T96XOm2sk+SAHOjYAaCGCI1DDVEaYCYtQsM8Ggcb4IZoDYDD5OCOP0lbyAE50M4BLZQdyCzAHD9+3ObPn28rV64MB7t\/\/3577733bPTo0WF56tSpNmTIEDt48KCx7Y4dO2zcuHHlKMyUKVPs9OnTtnfvXmtubrYLFy7YtGnTwr6TJ08Oy7zn1q1bw3a0hZV6kQNyoMsOADQOMx6hKe765MF0TPgFbAoTHy4\/iZcUFJODSUMhojfIYScNeAQ9XR4W7SAHcu1AZgEm6XoMLKwbPnw4hZ08eTKUvAA0w4YNCxAzaNAgO3v2LM1lOfyMGDHCLl26FMDHV9LmdZVyQA7cmgMONMVdz9ml5m8YqSfuQAJuEBEbUlDMrSm2gQ4RHIedNOCJoQfYQQ48RHZQYeInXzVwa0evvW\/ggFbJgcw4kHmAIVXE\/JYVK1aESAkRE3ePqMr58+fD4rlz5+zy5cuh7i+nTp0KVSI0AEtYaHuJ4Saut62+rrj\/\/vuNVBMaM2aMNTQ05FaDBw+22267LUS48tzPZN\/U7yG9dk4PtqINbDlqA07usn5vbi2r7651Zs2lqGtJ1178ol3Z8FDQ1VId2ctPWp\/X15X2ez1oUOl9ht71GWuY+ki76I5HdYCahqnzU\/ul8e698U7+rFVjuZ7Gm2sU1yo0cuTI665neWnok\/WOAB+kkpjnQlSFibnVOObHH3\/cNmzYEPTYY4\/ZhAkTcqvx48fb2LFj7b777sttH9PGT\/3O7niPv\/PThu4dctl++1cnbOwHu4PGnNhmd\/7oOzZ27zft3p3LQzni7c3WcHafDSgUAtTc8bv\/KvU81nhnd7zTfj5vta2exptrlF+vVq1aVY1LZq98Zp9e+ZRu+BBAhnkuTMz1ybieMuLtSSn179+falmeMmpsbDTgp7yiVIlTRnG9tCr1\/+rVq23p0qVBL7zwgh0+fDi3euutt+zYsWN25MiR3PYxbfzU72yO95F3f2Ho7Qv97ae\/MaasE2PmG3r3s1+xY9O\/bm\/PWhfKs\/cuspYR06z144+N9NSZV76Veh5rvLM53mk\/m93RVk\/jzTXKr1fPP\/986jUtD401AzCYDZB42igtZeTrgJ1kyoj9PaWUljJKa2Mf15kzZ8IdUdwVdeLECWtpacmt8PGjjz6yixcv5raPaePXiX7n0o\/e7vdFG2io2DDOWu9+IOjKxAV29YHlZk1rgvo8+gPrt\/S1oL6lOrKHn7FrDy4PuvLZhVa8fZxdKr3X+z\/7iV1855+M+TTMq0HMs2HOTXHXc1Y8dzx13Hq732nnXDXa1O\/8\/17jGsW1Cr3xxht+GctdmVmAIWrC3Bd\/RgtRF\/J5R48eDZNvuYOotbXV5s2bFybtzp4923wdo3TgwIEwZ4X9mpqawnyOffv2sSrcuUQkh1unEXXuZgor9SIH5EC3O8BdRUy4ZW4Kdx8hJuXy\/BgejFeY+LD1bfgkV89EXp74W9z1bICSGEyY\/OtwAqAwOZhJwohn1XBrN+r2DugN5YAcyJwDmQUYoijPPPNMgBAeNrdu3bowiddvq2b9li1bwrNheM4Lzm7cuJEiaP369WF79mP+zPbt28Mt1axkIjBkysRgRJ021mVaOjg5UMMO+NcS0IUYUrj7CFhp\/XFzKeXzyfchsVwsRU+AEgSUIO5uYn9JDsgBOZBZgGFoeIbLnDlzwoPseJjd4sWLaS4L6Jg1a1ZYz0RfoKa8slRhe\/ZDAE2pqfwfEKIdUS+vUEUOyIFud4AISbEUUWl9s9mIrlz78LT1uX2UEXkhAkMkBhGVQURoiNYgbpd2EcVxEdVB3GnkKkxsKr3nJ7dS89wYxAPzULd3Sm8oB+RAVR3oCsBU9UD14XJADtS2A8VSRAWQIe2DSAEhUkKI9JCLeSykjhARGheRG4QTpJxQoZR+Qmkg5DCUBkTMt2HCb5+pf2gOQAAP4v0lOSAHsu2AACbb46OjkwN14wDpIRfpIlJHCPBxAUAIAHIBQQgIcjkIUQJDCBhCDkMYy0Rgu3dOCWC+HATwICI+rJfkQPc4oHfpCQcEMD3hqt5TDsiBXnOA9BAiclJoSyENKJWIqAxtRFgG\/G5TgJTCxE++0oA5OVcGDL1u8jBprmIp3dVrHdAHyQE5cFMOCGBuyjbtJAfkQG85UCjBCADCfBfEHBgiJMyPQaSHEJET5tGgQglSECkmRNoJMf+GCAyRGJ7oy0PwWl\/6qhG5IYpDVIcID5Gf3upfb3yOPkMO5NEBAUweR1V9kgM5cuCTKMrDRsTEYQQQIUoCiJAeQqSLABEEjCCABAElruKu5wxA+fjHzTbww6PW99w\/58gtdUUO1I8DApj6GWv1VA5UyYFb+1hAxBXDiIMI82UQ82du7ZO0txyQA7XkgACmlkZLxyoH5IAckANyQA4EBwQwwQa95NkB9U0OyAE5IAfy54AAJn9jqh7JATkgB+SAHMi9AwKYHh9ifYAckANyQA7IATnQ3Q4IYLrbUb2fHJADckAOyAE5cOsOdPAOApgODNJqOSAH5IAckANyIHsOCGCyNyY6IjkgB+SAHKi+AzqCjDsggMn4AOnw5IAckANyQA7IgesdEMBc74la5IAckAPVd0BHIAfkwA0dEMDc0B6tlANyQA7IATkgB7LogAAmi6OiY5ID1XdARyAH5IAcyLQDAphMD48OTg7IATkgB+SAHEhzQACT5oraqu+AjkAOyAE5IAfkwA0cyDTALFiwwHbu3Gm7d+8O2rRpU7uuxOu3bdtmjY2N7dazve+7bNmyduvWrFkT3pP11Nut1IIckANyQA7IATmQaQcyCzDTp0+3J554wl566SWbOXOmLV++3EaNGmUOG8DKwoULbc+ePbZo0aJg8pIlS0LJC8DC9uy3efNmmzt3rvGerAN8Jk2aZGvXrg2iThvr2qRCDsgBOSAH5IAcyLADmQWYvXv32pw5c2z9+vXBPpZPnz5tI0aMCMtTp061IUOG2MGDB+348eO2Y8cOGzduXDkKM2XKFGN79mtubrYLFy7YtGnTwr6TJ08Oy\/v377etW7eG7WgLK\/UiB+SAHJADckAO3KQDvbdbZgGmIwuGDx8eNjl58mQoeQFohg0bFiBm0KBBdvbsWZrLGj16dKgDQZcuXQrgExpKL7SVCv2XA3JADsgBOSAHasCBmgEYUjxjxoyxAwcOlG0lqnL+\/PmwfO7cObt8+XKo+8upU6dClQgNwBIW2l5iuInrbauvK+6\/\/34j1YQ4joaGBsurBg8ebLfddluIcOW1j2n9Ur+H5Pac1ng3lMc27+d52ljTVk\/95hrFtQqNHDnyuutZXhpqAmCYu8J8mBMnTpRTSr09AI8\/\/rht2LAh6LHHHrMJEybkVuPHj7exY8fafffdl9s+po2f+q3xTjsv8tam8zz\/5znXKL9erVq1qrcvl732eZkHGOBl9erVYZ7K4sWL2xnjKSMaSSn179+falmeMmpsbDRSSuUVpUqcMorrpVWp\/zmGpUuXGnrhhRfs8OHDudVbb71lx44dsyNHjuS2j2njp35rvNPOi661Zf\/3gs7z\/J\/nXKO4VqHnn38+9ZqWh8ZMA4zDy6FDhywJL2kpI08ppaWMGCxPKaWljNLa2Md15swZ4zgQkaCWlhbLq\/Dxo48+sosXL+a2j2ljp35rvNPOi7y16TzP\/3nONYprFXrjjTf8Mpa7MrMAQ9TkySefDJGXlStXXmc8dxC1trbavHnzwqTd2bNn29GjR8sTc5krQ\/4PCGpqagrzOfbt2xfehzuXuMWaeTWIOm1hpV7kQA07oEOXA3JADtSLA5kFGG6TvuOOO+yee+4pP3COh875A+uIsmzZssVmzJhhPOeFAdu4cSNFELdfcxv1unXrwnNitm\/fbtxSzUpunYZMV6xYYYg6bayT5AhYqdEAABAASURBVIAckANyQA7Igew7kFmAAShmzZoVHmLHg+xc8+fPL0dZ4m3idredtJPvB9B4OyVRHV9HnTbpVh3Q\/nJADsgBOSAHeseBzAJM73RfnyIH5IAckANyQA7UogO5AphaHAAdsxyQA3JADsgBOdB1BwQwXfdMe8gBOSAH5IAcyJMDNdkXAUxNDpsOWg7IATkgB+RAfTsggKnv8Vfv5YAckAPVd0BHIAduwgEBzE2Ypl3kgByQA3JADsiB6joggKmu\/\/p0OSAHqu+AjkAOyIEadEAAU4ODpkOWA3JADsgBOVDvDghg6v0MUP+r74COINcO9GkYmev+qXNyoFoOCGCq5bw+Vw7IgVw6ALAUJjbZoKY\/t4Y\/2GhDl2432nLZWXVKDlTRAQFMFc3PyEfrMOSAHOhGBwaXwAUNKEHM1ZYzdrH5G3atVHbjR+it5IAcKDkggCmZoP9yQA7Ige5yoLjrOQNceD8gZuADX7ZPlaIwRGOIygx84AkjQtPvrsmGFJ3BKUkOdN2B6gNM149Ze8gBOSAHMuvAlZ8dtJYX\/8Te3zA3RF8+fvNlQ0BN34aRBtAQoWl4dKMhUkx3fG2\/IUAHxbAD8Ah0MjvcOrAqOiCAqaL5+mg5IAfy6QApI\/Txm81GRAZdKqWSWv7zEnvv6alBAA5qeXFJAJ3W0rZXTh60ax+etn53TzaP3gA8DSXYcdBxwFHkJp\/nTld6Ve\/bCmDq\/QxQ\/+WAHOg1B4AOoitBj363FIFBG42IDMCCgJe0AyKCg3xd34ZRXlUpB+rSAQFMXQ67Oi0H5EBvOdCnlDZizgvzXxpK0AKEINJKiEm+LqIxLqIzyCM2H5RSUi4iOaSqeqsP6Z+jVjlQXQcEMNX1X58uB+RAlR0AMGIx38QFeMRiPooLIHGFiMofbAy3TVMnzeMi9UMaiMjKpZe\/YaSSEGklRJrJBZS4SEGhKtujj5cDmXVAAJPZodGByQE5cCMHKq1zGAFCHD6AjiRs9Hn0B\/b2rHXWt1QCGa6GR0sg0iZSO7EKEx82F0DCpFzEsRBVcTGXhegKIrrikRPghG0lOSAHbt0BAcyte6h3kANyIGMOACOAiMMHERDmlyDAAzl4+KEDHS4m1BZ3PWsIACGtQzrHQcRLUjkuoiqxirueK+3\/nBFd8c9QKQfkQPc5UBMA09jYaNu2bbNly5a16\/mCBQts586dtnv37rCe7eINNm3aFNaxPrnvmjVryuuox\/upLgc6dkBbVMMBj67cqOS4gI5iG4AAIw4mlB4loWRbF1DjAnSAHgQEAUMNj37XN1UpB+RABhzIPMAAJd\/61rds1Kj2M+5pX7hwoe3Zs8cWLVoUrFyyZEkoeQFY2Gf58uW2efNmmzt3rk2fPp1VBvhMmjTJ1q5dG0SdtrBSL3LgBg7c6MJZaR2pjM6qMLGplKLomkiPdFaeRonLqw8st7P3LrJ+c1aGx9\/H66gzp6Mr8rkfnSl59klXRGSlMwI6gA8EjDiYUBJ5cVUaauAGAT8IIGL+SqXt1S4H5EDvO5BpgAEqvv\/97wdXLl26FEp\/mTp1qg0ZMsQOHjxox48ftx07dti4ceMMsGGbKVOm2OnTp23v3r3W3NxsFy5csGnTprHKJk+eHJb3799vW7duDdvRFlbWyEstHCYX9OSFlQtiV9WViyfbdubC6dvEF89+S18LcyIo4\/a43pmLZ3Kbhkd\/PaeiozoX3q6Ki3QhmptxozoXcJdfxAcOb7QrA4balcJQ87a49HONC3pnRJQDMf+jI3mUhBJI6IpI63RFpIBiXX3xizZ27zeNMm7nrh9PEV1q\/mTSLWkgzV\/xM0GlHMiGA5kGmJaWFnv66aftmWeeuc6t4cOHh7aTJ0+GkheAZtiwYQFiBg0aZGfPnqW5rNGjR4f6iBEjDCACfEJD6YW2UqH\/3ehA34ZR5hdTLrKIv4a7Kr\/gdraML74d1TvT3UoXbS7SXRF\/yXOh7qw6upjHF2+iA11VsW2ORrGUahnx9mbrs2tdeKBapc9lu47kF3zKYvn9P5kL0tEykNAVARRdEXf0JNWv9X1LtnXmnNA2ckAOVN+BmwSY3jnwv\/u7vzNU6dOIqpw\/fz6sPnfunF2+fDnU\/eXUqVOhCqgALGGh7SWGm7jetvq64v777zdSTWjMmDHW0NCQWw0ePNhuu+22EOG6lX4ObDlq10p\/5aIrGx6yNPHXb5rs5SctTX1eX2dJ9Xtzi1VS4eQuS2rgh0ctTQ3Fn9m\/uHImrOPC5hpQKFiaiFx0pBi6ADcgrrPqKBLTUTSns+sHPPJtOzb969fdjZOMJHVmOY5W3Wrdo2TdURKZS8qa1ti7n\/2KFUr9T66Llz\/10LJc\/ax318\/3rfxuqMa+9dRvrlFcq9DIkSOvu57lpSHTAJMlkx9\/\/HHbsGFD0GOPPWYTJkzIrcaPH29jx461++67r8f7OP7OT1ua7h1y2dL02786YUmN\/WC3VdKYE9ssqTt\/9B1L07888l2bdfG\/2meOvhBSC6QXblX37lxu3aVbPZZ4\/7j\/9HfWhf9q\/\/LIs6m+xNum1YnexPr0iR9arIaz+yxWGjwCi8mf974NI1NTWjfTHoNkXC\/ePi58RhJE4wvs0Ls+0+M\/B735+6Q3f757s18dfVZm+t0L1w6uUX69WrVqVfJHKzfLNQ0wnjJiNEgp9e\/fn2pZnjJqbGw0UkrlFaVKnDKK66VVqf9Xr15tS5cuDXrhhRfs8OHDudVbb71lx44dsyNHjuS2j2njl\/V+v33gdesuHfvHvzXXyd0v28V3\/of9\/MB\/K7f5uo7K87d\/NkwAZhKw6xdjPm+xWkZMs1hAQ1LMwUn9wYsa41RenLorllJgScXpNerF7y2wWEQDf\/MH82zSgf9ofbd+qd06tvNboynf2\/a1XP0cZP08T\/vZ7I62euo31yi\/Xj3\/\/PPRT1G+qjULMGkpI08ppaWMGDZPKaWljNLa2Md15swZO3ToUNCJEyeM+Tl5FT5+9NFHdvHixVz3Mzl+6nfXx\/u9H\/+9xfDAPJ+kYtiI6zGQeN1\/3tLKOPISR1HSUnLJ9NnAL221WDzE7tLcDfbjCf\/BPvo3f2ksp4k0U0fiLq7O6MrEBdYZtd79gHWkYsM464wu2kBLyH5x+X8z\/Xy35Pp3G9cov1698cYbaT9OuWirWYDhDqLW1labN29emLQ7e\/ZsO3r0aLgjiZE5cOBAmLPCrdNNTU1hPse+fftYFe5c4hZr7nJC1LmbKazUixyQA5124OM3m0sA81xZTN5NiihGmvxOn7jkDqA0xXcJxXWiK2nqaCIyd0fZ26\/akFJ6i5LlWA5aDlZppZsUg1Vcd8hi7hNKA620to7mPrE+CWiVltPmLQFyzP3x41cpB2rVgZoFGKIsW7ZssRkzZoTnvDAAGzdupAhav359uD163bp14Tkx27dvD7dUs5Jbp6HTFStWGKJOG+skOZAXB\/o0jLQ0JZ9Jc3X47xgpnX53TbLCTTyHhn2St8tXXH7gCYvXdfaW+uSFnos4it8rrg\/43SZLU2Hiw6U+fiK7d45dKKW3KON26g4fXgIgSfm6SqXDTHefT0mYctiiTEa\/4uiYQx2T4Iec++SPue4+Nr2fHOhNB2oCYHiWy5w5cwwoic0BOmbNmmUzZ860+fPnl6Mvvs3ixYvDOtYn9125cmV5HXXfR6Uc6AkHChEYxBdar6ddyOM7YeJ62l05aXf8pP31TVvyr3XuQgp\/kT\/8jAEFNyMHjEIECJ2pxxd\/v+B3pvQxSl7MO1rmIu8i8lIpAhNHY6g7CDgEdLVMixJ11BZHmuKoVByxoh5Ht5LRr2J0KzvRsqAfN4dJ1e6hSjlQqw7UBMDUqrk67sw5UJUDKpTghQu8g0Eh5SLPhTx54eZg0y7IXIC5qMYq7nq24jNc0i6U8cXx6ovXP9AtvmB2pc4F9WYVX4g7W09esLuyfG3\/X9sdJ39olMXoQn+jegCAUtqsq2VXnlfj28bPp+FckOSAHGjvgACmvR9akgM94gBwUixBBorBgzptqPXHzeYqtl1Q4wtlsW1\/Sk8VFNu2o4y3jet+QYzL+OJInduYKV09YoLeVA7IATnQjQ4IYLrRzA7fShvUpQPABFGQZOSFqAzyyExcJtM8LJP+SSotdURbWpqJtjgV5XXuoOEW6OR3IXl6q1JZKEWWKqnfXZOtktLm5XhbXZ4g6rQckAM35YAA5qZs005yoGsOEP1IplYqpWbi9E5cB4KSSpuLUSxFaojsJEXqKZmSohc8xI3nsCS\/C6mQkuqiDehCMXAl6wBXJSUhLF4GvioJALuRHMgoHbqAKOCIfkpyQA7UrgNpRy6ASXNFbXKgig54GidZAkFJEd1JqhilleJ62vwQ5pp8\/NJXwxN4rXmlsVzWi39iLSnqClQVSzAVy1NfXgJVriRcsRwPA2m4G4l5RC4ACwFRABJQBNTE76e6HJADte2AAKa2xy83R89fyTcj\/sLuqgoVUh92z0PhibGF322ytG24AHZFaXcWdaaNCEJXdaPIRKV1XNQR37799qx1Rsmyiwt\/RwIQkoqjMUBErLRbkfvcPspQ2skMxCCHnLTSYehGJRB15eShtI9QmxzohAPaJIsOCGCyOCptx+QX9I4u0GkX22TbjS68yYtq2pyIji6olS6S3u4XxUplRxfKSuuTF8\/OLMcX2Ljeb86fhcfiX3tweertxFyICxXSKmntHg2gvFHkILmubfiNC3dnFV\/Yk6mjSstc1BHPBeG7jK68+hepdzIl01Zpy3FUJlmvlCrz9mRqLW25HBX6z0ssWU+LLCXbiqWoFNEr91alHJADte+AACbDY9jw6HeNC3fDoxvtRoovwpXqXHwrKflXMct8d43dM8eoIy7CiL+U0+Q2xhfSuM5fx8VSOiFtzkayLe0C2VFb8qLZmWW\/gHrJ9+PwxYuU3pYs0y6unWlLXnQ7s5y8CHd2uVi6WHdFH\/+47bkg77xmyXQUy1z4O1Iy3RUv+7mh8uYd0J5yQA5c74AA5npPMtNy6eVvmF+4kxd5loslIKgkgCGpGCjievKv\/BsZkIwWJJeBnUoCoCoBVtw+6OE\/t46UjCjF+3udz7ux2j8Vts\/UP7T37v68UcbvX0hJOREdu5FPWicH5IAckAM964AApmf9vaV3j\/\/q5S\/hpIo3+Es77a\/1Sn\/1JyMIRCDSIhE3imrcCLSALVSMgCsJVyw7VF378LShSuYBTUSDOlIlkPL2JNz0nfaH4VuUKeN1DkRxyXtUOr78tatHckAOyIHsOSCAqfKY8Jd8rHi+S\/yXfxwRSM5ZqbTc0bwVX+\/zVLzs8+gP7Nj0rxultzF3hXRWJXmKK77IJ+sxFAAASSWBBFBJU3LIkhEklh2GvASQYsUwRR3AiueCOJB5mYS3Ygkek8ehZTkgB+SAHOg9BwQwPeQ18AFYOCRQxjAAEKAkEDgIUMYAUIgmj8YX+rQLvLd517ig30h+kff\/dremAAANBUlEQVTJnmnfEVNsi55woe+M\/MLfUZkEg7Tl5DyU5HIygsRyMtqUjEgVSwASK0S3orkgcfSLejyng7p7q1IOyAE5IAeq44AApgd9BySY8OrAwTLq6COBDaCCiEGxDRwcLnyZkjkyxbYLsUNF8sLNcvLiXWnZ34vvhqn0HTHhQt+J74Lhot8ZAQMdqSO\/tF4OyAE5IAfqzwEBTIdjfnMbcKEHHogGeMQgji4QmXDoKLZBCiXQ4vM\/gJ048kIKJo7KEKVxeSSHqI7LIz5Ef26uF9pLDsgBOSAH5EA2HRDA9OK4xJEGohNADiq2RVEoiY4APi4AyOUgRJmEIYAIORRRErWhe0SBKCU5IAfkgByQA73mQA9\/kACmhw3uqbdPwhBAhAAiVwCil79hAFBPHYfeVw7IATkgB+RANRwQwFTD9V78TKCmFz9OHyUH5IAcyIoDOo6cOyCAyfkAq3tyQA7IATkgB\/LogAAmj6OqPskBOVB9B3QEckAO9KgDdQswa9assd27dwdR71GX9eZyQA7IATkgB+RAtzpQlwCzYMECmzRpkq1duzaIOm3d6mwNv9mnPvUp+9f\/+l8bZQ13o8uHTn9z1O9O91\/9\/lSnvcrDhhrv+hrvPJyzlfpQlwAzefJku3Dhgu3fv9+2bt1qp0+fNtoqmVRv7UOHDg0AQ1lPfae\/AAyl+p1\/BxhnjXf+x9l7WK\/j7f3PY1mXADNixAi7dOmSHT9+vDymtJUXVLl1B\/QOckAOyAE5IAd60IG6BBj8PHv2LEVQXA8NKS\/3339\/SDuRbsq7Ro4cGRyopz4zpup3\/ZzjGm8z\/XxPqovf6f57LfxSz9lLXgGm24bpzJkzdujQIXv88cdtw4YNdaFVq1YF\/+qpz4yt+l0\/57jG2+rqd1q9jzfXMK5l4Rd7jl7qFmDilFFcT44tg7569WpbunSpJA90Dugc0Dmgc6DmzgGuYVzLkte39OXaaa1LgElLGaW1+TAy8BCsdChEo+SDfNA5oHNA50DtnANcw\/x6lqeyLgHm4MGDNmrUKOPWaUSdtjwNrPoiB+SAHKhFB3TMcqCzDtQlwHDrNH89rFixwhB12jprmraTA3JADsgBOSAHqutAXQIMlq9cudJmzpwZRJ02SQ7IgXp3QP2XA3KgVhyoW4CplQHSccoBOSAH5IAckAPXOyCAud4TtciBqjmgD5YDckAOyIHOOSCAKfnERN6dO3eGL3bctm2bNTY2llqv\/8+XPvoXQMYl7ddvnf2WzvabnsTb4hXLtNeiOHb6wBjeaLzp27Jly8J5wbavvvqqTZ8+nebciXMeL+hv7jrX1iH6Rh\/pa1tTxaIr21Z8k4ys6ExfNm3aVD7Pd5Z+F\/IzkpHDv+nD6Ey\/+d3NzzbKy893Z\/rtpvKzwM8EPnhbLZV1DzAM4MKFC23Pnj22aNGiMHZLliwJZfKFuTI+b4Zy8+bN4SsJ\/st\/+S\/JTTO\/nN7v9H5z0X7iiSfspZdeCnOG8Ipl2jPf0cQBdqXf\/BJ\/5JFHjHFmvLdv325PPvlkRcBNfFTNLOLJt771rXBnXs0cdBcPlF\/q\/vPd0a5d2baj96r2+s70hYvXoEGDwu8\/zvNa\/vl2vzvTb7bhiczLly8Pv9f4Trxa\/\/mmT509z\/GKax134VKvRdU9wEydOtWGDBli3EbNdyPt2LHDxo0b1+FFiov33LlzjYva3r17a27su9Lvu+++O\/Tv3LlzocQrKt5OvVbUlX7zBZ\/vvfeeNTc3h+5RXr582XiP0JCDFyDt+9\/\/fugJ3w8WKjl7IboAiJ44caLDnnVl2w7frMobdLYv\/GE2f\/788nfD1fLPN5Z3tt\/r16+3OXPmmP\/+PnDgQLgWDBs2jLepOXW2394xrmETJkwIf4R7W62VmQKYapg3fPjw8LEnT54MJS8ATUcn8bx58+zChQvlixv71ZK60m\/3xvfhwk7f9+\/fX0tdDsfqffA+0diZ8WY71L9\/f\/P3YLnW1dLSYk8\/\/bQ988wztd6VisfPhWnWrFn2s5\/9rOI2vqIr2\/o+WS3z1JeueHyz\/Z4yZYoRhXGg6cpnZmHbrvb73\/27f2e7du0K17EsHP\/NHEPdAwymcTE+f\/48VSPKwF\/ZYaHCC+RK6HFHKVpD1KbCZplv7my\/+YH+8pe\/bLNnzw55cr56If6LLfMdTRxgZ\/t96tSpkFZpamoK70BJuHX06NFhOQ8vf\/d3f2coD32p1Af+0q60LtnelW2T+2Zt+Wb78oUvfCFc1GrxDxTGoKv9Ju3CHJgxY8bYK6+8wlvUoqwr\/Sbyyh9j\/\/AP\/1CTffWDFsC4E10op02bFrbet29fKPP+wslOmgFgI0cO6edlwtuNxo5fCK+\/\/nqYG8AvuM997nPhL7Qb7aN1cqCWHWA+DBfyLVu2lFNKtdyfzhw7P+f8XmOOH995x++7zuxXq9sw5+3f\/tt\/a3\/zN39j\/od7rfZFAFMauTiFQHoAMi01V\/zPX+C1HGr0jnW236SMknNBiGKQRvP3qqWys\/2mT8wP4JcbWrduXciRE5lhnSQH8uQA8DJjxozwjftdfjJ5Doxgjhu\/5\/h9l4PuVOwCkWSyDHkY47oHmLSUERfnSmQKvTLJ90Zf\/ljxzMnQiq72O0OHfkuHciv99knLvMctHYR2lgMZcwB4IS3OXTh5uLDdjL3Me+SPm7z\/gcJcn3vuuSdMB+AOS9LiDz74oG3atOlmbKvqPnUPMOR5W1tbjWgCcMI8j6NHj1YMn3KSDxgwINy1VNWRu8UP70q\/uSvhjjvuMMidj6VkmXaWa0ld6TehZHLizHni3OB2eyJv9foLPsPjrEO7BQeYA0Lk5bnnnivfkXMLb1czu9LvOBXONYCDz\/vUgMWLF4fbxokqc8s1v9NIldNO\/2tJdQ8wTMIl38sPMDTK4G3cuJEiiL9MYjLlr\/COUkxhx4y\/dKXfXLDJD3OyMxeEW1I3bNhgtGe8m9cdXlf7ffjwYSN15OfGU089dd17qqF2HUj+fNduT7p25HG\/+Yuc32l8sS0\/3y4u8F171+xvHfebuS98kS8\/3\/SZa0BeIS7ud\/ZHqfNHWPcAg1VciLnNEiJN3l3DHIiYTNmWZwdQsm8tiz50tt\/8sOMPYh\/2rdW+c+z0gb50NN6MP9uh5Lbl\/uegwp1mnNeMcw66k9oFxjI5hrTFP9++I+3JbX1drZVpfaHN+03J+Z1UrZ8L9DE5hrTRXx9Dlr3f\/E7gd4Ovq9WSPnXUb+8bf9CxLft4Wy2VAphaGi0dqxyQA3JADsgBORAcEMAEG2r6RQcvB+SAHJADcqDuHBDA1N2Qq8NyQA7IATkgB2rfgVsHmNr3QD2QA3JADsgBOSAHaswBAUyNDZgOVw7IATkgB\/LhQHf0gsc8cKcsj3rg\/Vjm9nDurKKdtrxKAJPXkVW\/5IAckANyINcO8Kwqvoh10KBB5X7yPBtuD+exF9wezzbllTmrCGByNqDqjhyQA3Kgcw5oq6w5AGz4wzPjYyOSQkQFxc\/nGTx4sD399NN26dKl8ubxk4T5yoCTJ0+W1+WtIoDJ24iqP3JADsgBOVBzDgAvfJkkUZP44AEWHve\/fPly44Gac+fONdJEbPO9733PPvzwQ6pl8fyeESNGhG354l2e8VRembOKACZnA6ruyIFacUDHKQfkwCcO8KRcnoR84sSJTxqiV56UzOP+ARG+cJLv6ps2bVq0Rfsq0RrAhQf0sS8A1H6L\/CwJYPIzluqJHJADckAOZNABoivxV9SwHKeK3nnnHWPOCuARHz4Tc5nfkvzy4NGjR8ebtav\/8pe\/NP\/C2eR+7TbMwYIAJgeDqC7cjAPaRw7IATnQOw7wFQW\/+Zu\/aTt37rT\/9J\/+kz3xxBP2zW9+s\/zlmaSCeKx\/paPxeS1sE893Sdv+7\/\/+741UFPNlxo0bZ0Rt0rbLQ5sAJg+jqD7IATkgB+RAph3gO5j27Nljv\/Vbv2Vf\/vKXy\/ByqwdNaon3Bm54L2CJ73UihcT3HHk76\/KmPnnrUK30R8cpB+SAHJAD9eMAc1NmzJhhP\/\/5z+3ZZ58tT8TtjAOeMvKUUmf2qYdtBDD1MMrqoxyQA3JADlTNAea8MDeFyMi\/\/\/f\/3p577jn7+te\/3iHEED1JSxl5SqlqHaruB5c\/XQBTtkIVOSAH5IAckAPd7wBpnSVLlpTfmOUvfOELnUojMbF30qRJAXaamppsyJAhtm\/fvvJ71XNFAFPPo6++ywE5IAfkQNcc6OWtea4Lt1GvW7cu3Km0ffv2ToFPLx9mVT5OAFMV2\/WhckAOyAE5IAfaOwCszJkz5zpAYZIuk3IR27Tfq36XBDD1O\/bquRyQA7XngI5YDsiBNgcEMG1GqJADckAOyAE5IAdqxwEBTO2MlY5UDlTfAR2BHJADciAjDghgMjIQOgw5IAfkgByQA3Kg8w4IYDrvlbasvgM6AjkgB+SAHJADwYH\/BQAA\/\/\/71U0AAAAABklEQVQDACGKFrZ81J4NAAAAAElFTkSuQmCC","height":337,"width":560}}
%---
%[output:693537c4]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('signal')\" style=\"font-weight:bold\">signal<\/a>\nToo many input arguments."}}
%---
