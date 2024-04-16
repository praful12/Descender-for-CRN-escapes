%take delta_x - smoothen it
%input f0 for smoothening frequency. Used in line 29

function delta_x_s =  smooth_traj_butter(delta_x,f0) 
traj_d = delta_x;
[traj_pts,num_spec] = size(delta_x);


%------- concat traj
concat_traj = zeros((traj_pts)*2-2,num_spec);
flip_traj = flip(traj_d);

%front to front
concat_traj(1:traj_pts-1,:) = -flip_traj(2:end,:);
concat_traj(traj_pts-1:(traj_pts)*2-2,:) = traj_d(:,:);
conc_old = concat_traj;


con_concat = zeros(traj_pts*4-4,num_spec);
con_concat(1:(traj_pts)*2-2,:) = concat_traj;
con_concat((traj_pts)*2-1:4*(traj_pts)-4,:) = concat_traj;

con_concat_old = con_concat;



%----- matlab butterworth filter

fs = size(traj_d,1);
%f0 used here. This is the filter frequency.
[b,a] = butter(4,f0/(fs/2));
%freqz(b,a,[],fs)

%
%------- filt filt
for i = 1:num_spec
    con_concat(:,i) = filtfilt(b,a,con_concat(:,i));
end

%{
figure()
plot(con_concat)
hold on
plot(con_concat_old,'--')
hold off
%}
%traj_s = con_concat(1:traj_pts,:);
delta_x_s = con_concat(traj_pts-1:2*traj_pts-2,:);

delta_x_s(1,:) = zeros(1,num_spec);
delta_x_s(end,:) = zeros(1,num_spec);
%{
for i = 2:floor(traj_pts/4)
    for j = 1:num_spec
        if abs(delta_x(i,j))<abs(delta_x_s(i,j)) && abs(delta_x(i,j))<2*10^(-4) 
            delta_x_s(i,j) = delta_x(i,j);
        end
    end
end
%}
%}
end
