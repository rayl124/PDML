close all;
clear all;

%%
strnum_l = '007';
strnum_r = '008';

input_string_l = ['Input',strnum_l,'.txt'];
part_string_l = ['NumberPart',strnum_l,'.txt'];
electron_string_l = ['Electrons',strnum_l,'.txt'];

input_string_r = ['Input',strnum_r,'.txt'];
part_string_r = ['NumberPart',strnum_r,'.txt'];
electron_string_r = ['Electrons',strnum_r,'.txt'];

energy_vid = ['ElectronEnergy',strnum_r,'.avi'];

%%
fileID = fopen(part_string_l); num_col = 5;
field_string = [];
for i = 1:num_col
    field_string = [field_string, '%d '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_l = field_data{1}; e_np_inner_l = field_data{4}; 

fileID = fopen(part_string_r);
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_r = field_data{1}; e_np_inner_r = field_data{4}; 
close('all');

fileID_l = fopen(electron_string_l);
e_data_l = textscan(fileID_l, '%f %f', e_np_inner_l(1), 'HeaderLines',1);
e_data_l = e_data_l{2};

fileID_r = fopen(electron_string_r);
e_data_r = textscan(fileID_r, '%f %f', e_np_inner_r(1), 'HeaderLines',1);
e_data_r = e_data_r{2};

v = VideoWriter(energy_vid);
open(v);
histogram(e_data_r,'BinWidth',1)
set(gca,'yscale','log')
ylim([0,0.5*max(e_np_inner_r)])
xlim([1,200])
ylabel('# of Macroparticles')
xlabel('Energy [eV]')
title(['Iter: ',num2str(iters_r(1))])
frame = getframe(gcf);

for i = 1:3
    writeVideo(v, frame);
end

for i = 2:length(iters_r)
    e_data_r = textscan(fileID_r, '%f %f', e_np_inner_r(i));
    e_data_r = e_data_r{2};
    %histogram(e_data_r,'BinWidth',0.1,'BinLimits',[1,1000])
    histogram(e_data_r,'BinWidth',1)
    set(gca,'yscale','log')
    ylim([0,0.5*max(e_np_inner_r)])
    xlim([1,200])
    ylabel('# of Macroparticles')
    xlabel('Energy [eV]')
    title(['Iter: ',num2str(iters_r(i))])
    frame = getframe(gcf);
    for j = 1:3
        writeVideo(v, frame);
    end
end
close(v);

