% Add path to ESFieldData.txt
clear; close all;

%% File info
strnum = '006';
input_string = ['Input',strnum,'.txt'];
field_string = ['ESFieldData',strnum,'.txt'];
EP_vid = ['ElectricPotential',strnum,'.avi'];
EF_vid = ['ElectricField',strnum,'.avi'];
Charge_vid = ['ChargeDensity',strnum,'.avi'];
part_string = ['NumberPart',strnum,'.txt'];

%% Input info
fileID = fopen(input_string);
header_string = [];
for i = 1:12
    header_string = [header_string, '%f '];
end
input_data = textscan(fileID, header_string, 'HeaderLines', 2);
e_np0 = input_data{2}; i_np0 = input_data{2};
ts = input_data{10}; nn = input_data{12};

%% Field plotter
fileID = fopen(field_string);
field_data = textscan(fileID, '%f %f %f %f %f %f %f', 'HeaderLines',1);
iters = field_data{1}; t = field_data{2};
x = field_data{3}; rho = field_data{4};
phi = field_data{5}; E = field_data{6};
fclose('all');

figure
v = VideoWriter(EP_vid);
open(v);
for i = 0:(ts/(iters(1)+1)-1)
   plot(x(i*nn + 1:(i+1)*nn),  ...
       phi(i*nn +1:(i+1)*nn),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min(phi), max(phi)])
   line([0.085, 0.085], [min(phi), max(phi)])
   hold off
   axis([0,0.095, min(phi), max(phi)])
   ylabel('\phi [V]')
   xlabel('x [m]')
   title(['Iter: ', num2str(iters(i*nn + 1))]);
   frame = getframe(gcf);
   
   for j = 1:5
       writeVideo(v, frame);
   end
end
close(v);

v = VideoWriter(EF_vid);
open(v);
for i =0:(ts/(iters(1)+1)-1)
   plot(x(i*nn + 1:(i+1)*nn),  ...
       E(i*nn +1:(i+1)*nn),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min(E), max(E)])
   line([0.085, 0.085], [min(E), max(E)])
   hold off
   axis([0,0.095, min(E), max(E)])
   ylabel('E [V/m]')
   xlabel('x [m]')
   title(['Iter: ', num2str(iters(i*nn + 1))]);
   frame = getframe(gcf);
   
   for j = 1:5
       writeVideo(v, frame);
   end
end
close(v);

v = VideoWriter(Charge_vid);
open(v);
for i = 0:(ts/(iters(1)+1)-1)
   plot(x(i*nn + 1:(i+1)*nn),  ...
       rho(i*nn +1:(i+1)*nn),'r-','Linewidth',2)
   hold on
   line([0.01, 0.01],[min(rho), max(rho)])
   line([0.085, 0.085], [min(rho), max(rho)])
   hold off
   axis([0,0.095, min(rho), max(rho)])
   ylabel('Charge Density [C/m^3]')
   xlabel('x [m]')
   title(['Iter: ', num2str(iters(i*nn + 1))]);
   frame = getframe(gcf);
   
   for j = 1:5
       writeVideo(v, frame);
   end
end
close(v);

%% Num particles plotter
fileID = fopen(part_string);
field_data = textscan(fileID, '%d %d %d', 'HeaderLines',1);
iters = field_data{1}; e_np = field_data{2}; i_np = field_data{3};
fclose('all');
plot([0;iters],[e_np0;e_np],'Linewidth',2)
hold on
plot([0;iters],[i_np0;i_np],'Linewidth',2)
xlabel('Iter #')
ylabel('# of Ion Macroparticles')
legend('Electrons','Ions')