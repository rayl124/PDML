clear ; 
close all;
nx1 = 16; nx2 = 10;
nn = nx1*nx2;

fileID = fopen('ESFieldData.txt');
field_data = textscan(fileID, '%d %f %f %f %f', 'HeaderLines',1);
rho_i_all = field_data{2}; phi_all = field_data{3};
E1_all = field_data{4}; E2_all = field_data{5};

fileID = fopen('NumParticles.txt');
np_data = textscan(fileID, '%d %d', 'HeaderLines', 1);
iters = np_data{1}; np = np_data{2};

fileID = fopen('ParticleInfo.txt');
particle_data = textscan(fileID, '%f %f %f %f %f','HeaderLines', 1);
x1 = particle_data{2}; x2 = particle_data{3};

fclose('all');
v = VideoWriter('ElectricPotential.avi');
open(v);
figure(1)
for i = 1:length(iters)
    phi2 = reshape(phi_all(nn*(i-1)+1:nn*i),nx2,nx1);
    contourf(phi2);
    colorbar
    title(['Electric Potential: Time Step ',num2str(iters(i))])
    line([5,5,7,7,5],[1,5,5,1,1],'Color','black','LineWidth',2)
    frame = getframe(gcf);
    for j = 1:30
        writeVideo(v, frame);
    end
end
close(v);

v = VideoWriter('IonDensity.avi');
open(v);
figure(1)
for i = 1:length(iters)
    rho = reshape(rho_i_all(nn*(i-1)+1:nn*i),nx2,nx1);
    contourf(rho,1e11:1e11:1.1e12);
    colorbar
    caxis([1e11, 1.1e12]);
    title(['Ion Density: Time Step ',num2str(iters(i))])
    line([5,5,7,7,5],[1,5,5,1,1],'Color','black','LineWidth',2)
    frame = getframe(gcf);
    for j = 1:30
        writeVideo(v, frame);
    end
end
close(v);

v = VideoWriter('ParticlePositions.avi');
open(v);
figure(1)
interval_start = 1;
for i = 1:length(iters)
    x1_local = x1(interval_start:interval_start + np(i) - 1);
    x2_local = x2(interval_start:interval_start + np(i) - 1);
    interval_start = interval_start + np(i);
    plot(x1_local,x2_local, '.')
    xlim([0, max(x1)])
    ylim([0, max(x2)])
    lambdaD = max(x2)/9;
    line([lambdaD*4, lambdaD*4, lambdaD*6, lambdaD*6, lambdaD*4], ...
        [0,lambdaD*4,lambdaD*4,0,0],'Color','black','LineWidth',2)
    title(['Particle Positions: Time Step ',num2str(iters(i))])
%     line([5,5,7,7,5],[1,5,5,1,1],'Color','black','LineWidth',2)
    frame = getframe(gcf);
    for j = 1:30
        writeVideo(v, frame);
    end
end
close(v);
