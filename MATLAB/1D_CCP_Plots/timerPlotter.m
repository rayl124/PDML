close all;
clear all;

%%
strnum = '004';

input_string = ['Input',strnum,'.txt'];
timer_string = ['Timer',strnum,'.txt'];


%%
fileID = fopen(timer_string); num_col = 9;
field_string = [];
for i = 1:num_col
    field_string = [field_string, '%d '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters = field_data{1}; coll = field_data{2}; push1 = field_data{3};
fluid = field_data{4}; rho = field_data{5}; phi = field_data{6};
E = field_data{7}; push2 = field_data{8}; total = field_data{9};

coll_bar = mean(coll)
push1_bar = mean(push1)
fluid_bar = mean(fluid)
rho_bar = mean(rho)
phi_bar = mean(phi)
E_bar = mean(E)
push2_bar = mean(push2)
total(end)
