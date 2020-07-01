% Add path to ESFieldData.txt
clear; close all;

%% File info
strnum_l = '003';
strnum_r = '004';

input_string_l = ['Input',strnum_l,'.txt'];
CCfield_string_l = ['FieldCCData',strnum_l,'.txt'];
NCfield_string_l = ['FieldNCData',strnum_l,'.txt'];
part_string_l = ['NumberPart',strnum_l,'.txt'];
average_string_l = ['FieldAverageData',strnum_l,'.txt'];

input_string_r = ['Input',strnum_r,'.txt'];
CCfield_string_r = ['FieldCCData',strnum_r,'.txt'];
NCfield_string_r = ['FieldNCData',strnum_r,'.txt'];
part_string_r = ['NumberPart',strnum_r,'.txt'];
average_string_r = ['FieldAverageData',strnum_r,'.txt'];

EP_vid = ['ElectricPotential',strnum_l,'+',strnum_r,'.avi'];
EF_vid = ['ElectricField',strnum_l,'+',strnum_r,'.avi'];
Charge_vid = ['ChargeDensity',strnum_l,'+',strnum_r,'.avi'];
ne_vid = ['ElectronDensity',strnum_l,'+',strnum_r,'.avi'];
nn_vid = ['NeutralDensity',strnum_l,'+',strnum_r,'.avi'];
nstar_vid = ['ExcitedDensity',strnum_l,'+',strnum_r,'.avi'];
ni_vid = ['IonDensity',strnum_l,'+',strnum_r,'.avi'];
%% Input info
fileID = fopen(input_string_l);
header_string = [];
for i = 1:13
    header_string = [header_string, '%f '];
end
input_data = textscan(fileID, header_string, 'HeaderLines', 2);
e_np0_l = input_data{2}; i_np0_l = input_data{2};
ts = input_data{10}; nn = input_data{12};

fileID = fopen(input_string_r);
header_string = [];
for i = 1:13
    header_string = [header_string, '%f '];
end
input_data = textscan(fileID, header_string, 'HeaderLines', 2);
e_np0_r = input_data{2}; i_np0_r = input_data{2};
ts = input_data{10}; nn = input_data{12};

%% Field plotter
fileID = fopen(CCfield_string_l);
field_string = [];
for i = 1:9
    field_string = [field_string, '%f '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_l = field_data{1}; t = field_data{2};
x_l = field_data{3}; rho_l = field_data{4};
phi_l = field_data{5}; nn_l = field_data{6};
nstar_l = field_data{7}; ni_l = field_data{8};
ne_l = field_data{9};

fclose('all');

for i = 1:9
    field_string = [field_string, '%f '];
end
fileID = fopen(CCfield_string_r);
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_r = field_data{1};
x_r = field_data{3}; rho_r = field_data{4};
phi_r = field_data{5}; nn_r = field_data{6};
nstar_r = field_data{7}; ni_r = field_data{8};
ne_r = field_data{9};
fclose('all');

figure
pos = get(gcf, 'position');
set(gcf, 'Position', [pos(1), pos(2), pos(1) + 3*(pos(3)-pos(1)), pos(4)]);
v = VideoWriter(EP_vid);
open(v);
for i = 0:(ts/(iters_l(nn)))
   subplot(1,2,1)
   plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       phi_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(phi_r); min(phi_l)]), ...
        max([max(phi_l);max(phi_r)])]);
   line([0.085, 0.085], [min([min(phi_r); min(phi_l)]), ...
        max([max(phi_l);max(phi_r)])]);
   hold off
   axis([0,0.095, min([min(phi_r); min(phi_l)]), ...
       max([max(phi_l);max(phi_r)])])
   ylabel('\phi [V]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1) + 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       phi_r(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(phi_r); min(phi_l)]), ...
        max([max(phi_l);max(phi_r)])]);
   line([0.085, 0.085], [min([min(phi_r); min(phi_l)]), ...
        max([max(phi_l);max(phi_r)])]);
   hold off
   axis([0,0.095, min([min(phi_r); min(phi_l)]), ...
       max([max(phi_l);max(phi_r)])])
   ylabel('\phi [V]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end
end
close(v);

v = VideoWriter(Charge_vid);
open(v);
for i = 0:(ts/(iters_l(nn)))
   subplot(1,2,1)
   plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       rho_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(rho_r); min(rho_l)]), ...
        max([max(rho_l);max(rho_r)])]);
   line([0.085, 0.085], [min([min(rho_r);min(rho_l)]), ...
        max([max(rho_l);max(rho_r)])]);
   hold off
   axis([0,0.095, min([min(rho_r); min(rho_l)]), ...
       max([max(rho_l);max(rho_r)])])
   ylabel('\rho [C/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1)+ 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       rho_r(i*(nn-1) + 1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(rho_r);min(rho_l)]), ...
        max([max(rho_l);max(rho_r)])]);
   line([0.085, 0.085], [min([min(rho_r); min(rho_l)]), ...
        max([max(rho_l);max(rho_r)])]);
   hold off
   axis([0,0.095, min([min(rho_r);min(rho_l)]), ...
       max([max(rho_l);max(rho_r)])])
   ylabel('\rho [C/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end
   
end
close(v);

v = VideoWriter(ne_vid);
open(v);
for i = 0:(ts/(iters_l(nn)))
   subplot(1,2,1)
   plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       ne_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([ne_r; ne_l]), ...
        max([ne_r;ne_l])]);
   line([0.085, 0.085], [min([ne_r; ne_l]), ...
        max([ne_r;ne_l])]);
   hold off
   axis([0,0.095, min([ne_r; ne_l]), ...
        max([ne_r;ne_l])])
   ylabel('n_e [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1) + 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       ne_r(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([ne_r; ne_l]), ...
        max([ne_r;ne_l])]);
   line([0.085, 0.085], [min([ne_r; ne_l]), ...
        max([ne_r;ne_l])]);
   hold off
   axis([0,0.095, min([ne_r; ne_l]), ...
        max([ne_r;ne_l])])
   ylabel('n_e [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end
end
close(v);
% %*
% v = VideoWriter(nn_vid);
% open(v);
% for i = 0:(ts/(iters_l(nn)))
%    subplot(1,2,1)
%    plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
%        nn_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
%    hold on
%    line([0.01, 0.01],[min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])]);
%    line([0.085, 0.085], [min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])]);
%    hold off
%    axis([0,0.095, min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])])
%    ylabel('n_n [1/m^3]')
%    xlabel('x [m]')
%    title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1) + 1))]);
%    
%    subplot(1,2,2)
%    plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
%        nn_r(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
%    hold on
%    line([0.01, 0.01],[min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])]);
%    line([0.085, 0.085], [min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])]);
%    hold off
%    axis([0,0.095, min([nn_r; nn_l]), ...
%         max([nn_r;nn_l])])
%    ylabel('n_n [1/m^3]')
%    xlabel('x [m]')
%    title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
%    frame = getframe(gcf);
%    
%    for j = 1:3
%        writeVideo(v, frame);
%    end
% end
% close(v);

v = VideoWriter(nstar_vid);
open(v);
for i = 0:(ts/(iters_l(nn)))
   subplot(1,2,1)
   plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       nstar_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   line([0.085, 0.085], [min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   hold off
   axis([0,0.095, min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   ylabel('n_{AR^*} [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1) + 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       nstar_r(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   line([0.085, 0.085], [min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   hold off
   axis([0,0.095, min([nstar_r; nstar_l]), ...
        max([nstar_r;nstar_l])]);
   ylabel('n_{AR^*} [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end
end
close(v);

v = VideoWriter(ni_vid);
open(v);
for i = 0:(ts/(iters_l(nn)))
   subplot(1,2,1)
   plot(x_l(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       ni_l(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([ni_r; ni_l]), ...
        max([ni_r;ni_l])]);
   line([0.085, 0.085], [min([ni_r; ni_l]), ...
        max([ni_r;ni_l])]);
   hold off
   axis([0,0.095, min([ni_r; ni_l]), ...
        max([ni_r;ni_l])])
   ylabel('n_i [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*(nn-1) + 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn-1) + 1:(i+1)*(nn-1)),  ...
       ni_r(i*(nn-1) +1:(i+1)*(nn-1)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([ni_r; ni_l]), ...
        max([ni_r;ni_l])]);
   line([0.085, 0.085], [min([ni_r; ni_l]), ...
        max([ni_r;ni_l])]);
   hold off
   axis([0,0.095, min([ni_r; ni_l]), ...
        max([ni_r;ni_l])])
   ylabel('n_i [1/m^3]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*(nn-1) + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end
end
close(v);

%%
fileID = fopen(['FieldNCData',strnum_r,'.txt']);
field_data = textscan(fileID, '%f %f %f %f %f', 'HeaderLines',1);
iters_r = field_data{1};
x_r = field_data{3}; E_r = field_data{4};

fileID = fopen(['FieldNCData',strnum_l,'.txt']);
field_data = textscan(fileID, '%f %f %f %f %f', 'HeaderLines',1);
iters_l = field_data{1};
x_l = field_data{3}; E_l = field_data{4};
fclose('all');

v = VideoWriter(EF_vid);
open(v);
for i = 0:(ts/(iters_l(nn+1)))
   subplot(1,2,1)
   plot(x_l(i*nn + 1:(i+1)*nn),  ...
       E_l(i*nn +1:(i+1)*nn),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(E_r); min(E_l)]), ...
        max([max(E_l);max(E_r)])]);
   line([0.085, 0.085], [min([min(E_r);min(E_l)]), ...
        max([max(E_l);max(E_r)])]);
   hold off
   axis([0,0.095, min([min(E_r);min(E_l)]), ...
       max([max(E_l);max(E_r)])])
   ylabel('E [V/m]')
   xlabel('x [m]')
   title(['Test: ',strnum_l, ', Iter: ', num2str(iters_l(i*nn + 1))]);
   
   subplot(1,2,2)
   plot(x_r(i*(nn) + 1:(i+1)*(nn)),  ...
       E_r(i*(nn) +1:(i+1)*(nn)),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min([min(E_r);min(E_l)]), ...
        max([max(E_l);max(E_r)])]);
   line([0.085, 0.085], [min([min(E_r);min(E_l)]), ...
        max([max(E_l);max(E_r)])]);
   hold off
   axis([0,0.095, min([min(E_r);min(E_l)]), ...
       max([max(E_l);max(E_r)])])
   ylabel('E [V/m]')
   xlabel('x [m]')
   title(['Test: ',strnum_r, ', Iter: ', num2str(iters_r(i*nn + 1))]);
   frame = getframe(gcf);
   
   for j = 1:3
       writeVideo(v, frame);
   end

end
close(v);



%% Num particles plotter
fileID = fopen(part_string_l);
field_string = [];
for i = 1:5
    field_string = [field_string, '%d '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_l = field_data{1}; e_np_l = field_data{2}; i_np_l = field_data{3};
e_np_inner_l = field_data{4}; i_np_inner_l = field_data{5};

fileID = fopen(part_string_r);
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_r = field_data{1}; e_np_r = field_data{2}; i_np_r = field_data{3};
e_np_inner_r = field_data{4}; i_np_inner_r = field_data{5};

fclose('all');
figure;
pos = get(gcf, 'position');
set(gcf, 'Position', [pos(1), pos(2), pos(1) + 3*(pos(3)-pos(1)), pos(4)]);
subplot(1,2,1)
plot([iters_l],[e_np_l],'Linewidth',2)
hold on
plot([iters_l],[i_np_l],'Linewidth',2)
xlabel('Iter #')
ylabel('# of Macroparticles')
legend('Electrons','Ions')

subplot(1,2,2);
plot([iters_r],[e_np_r],'Linewidth',2)
hold on
plot([iters_r],[i_np_r],'Linewidth',2)
xlabel('Iter #')
ylabel('# of Macroparticles')
legend('Electrons','Ions')

figure;
pos = get(gcf, 'position');
set(gcf, 'Position', [pos(1), pos(2), pos(1) + 3*(pos(3)-pos(1)), pos(4)]);
subplot(1,2,1)
plot([iters_l],[e_np_inner_l],'Linewidth',2)
hold on
plot([iters_l],[i_np_inner_l],'Linewidth',2)
xlabel('Iter #')
ylabel('# of Macroparticles')
title('Inner Particles')
legend('Electrons','Ions')

subplot(1,2,2);
plot([iters_r],[e_np_inner_r],'Linewidth',2)
hold on
plot([iters_r],[i_np_inner_r],'Linewidth',2)
xlabel('Iter #')
ylabel('# of Macroparticles')
title('Inner Particles')
legend('Electrons','Ions')

%% Average Plotter
fileID = fopen(average_string_l);
field_string = [];
for i = 1:5
    field_string = [field_string, '%f '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_l = field_data{1}; n_bar_l = field_data{2}; ex_bar_l = field_data{3};
i_bar_l = field_data{4}; e_bar_l = field_data{5};

fileID = fopen(average_string_r);
field_string = [];
for i = 1:5
    field_string = [field_string, '%f '];
end
field_data = textscan(fileID, field_string, 'HeaderLines',1);
iters_r = field_data{1}; n_bar_r = field_data{2}; ex_bar_r = field_data{3};
i_bar_r = field_data{4}; e_bar_r = field_data{5};

figure;
pos = get(gcf, 'position');
set(gcf, 'Position', [pos(1), pos(2), pos(1) + 3*(pos(3)-pos(1)), pos(4)]);
ax1 = subplot(1,2,1);
plot(iters_l/200,e_bar_l,'Linewidth',2)
xlabel('LF Cycle #')
ylabel('$\bar{n_e}$','Interpreter','Latex')
title(['Test: ',strnum_l])

ax2 = subplot(1,2,2);
plot(iters_r/200, e_bar_r,'Linewidth',2)
xlabel('LF Cycle #')
ylabel('$\bar{n_e}$','Interpreter','Latex')
title(['Test: ',strnum_r])

linkaxes([ax1,ax2]);


%%
% 
% figure;
% semilogy(iters_l, residual_phi,'Linewidth',2);
% xlabel('Iteration')
% ylabel('Residual')
% title(['\phi Residual for ',strnum_l,'+',strnum_r])
% 
% figure;
% semilogy(iters_l, residual_E,'Linewidth',2);
% xlabel('Iteration')
% ylabel('Residual')
% title(['E Residual for ',strnum_l,'+',strnum_r])
% 
% figure;
% semilogy(iters_l, residual_rho,'Linewidth',2);
% xlabel('Iteration')
% ylabel('Residual')
% title(['\rho Residual for ',strnum_l,'+',strnum_r])

clear all;