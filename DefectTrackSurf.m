function [charges] = DefectTrackSurf(r, v, dr)
%%%%% finds defects in a vector field v, sampled at points r.
% v has dimensions (number of sample points)x3x(time length)
% r has dimensions (number of sample points)x3
% dr is the distance at which to consider neghboring points. if not
% supplied an estimate is done for a quasi-uniformly sampled spherical
% surface.

[l_grid,~] = size(r);
[~,~,t_max] = size(v);

if nargin < 3
    dr = sin(2*pi/sqrt(l_grid));
end


t_arr = 1:t_max;


% find the average distance r between two neighboring points
%dr = sin(2*pi/sqrt(l_grid)); %% nearest-neighbors
%dr = 3*sin(2*pi/sqrt(l_grid)); %% larger


% represent vector field in the regular grid
%vx = zeros(l_grid-2, t_max);
vx = squeeze(v(:,1,:));
vy = squeeze(v(:,2,:));
vz = squeeze(v(:,3,:));


%  for each point: 
%   i.   find nearest neighbors
%   ii.  compute local basis vectors
%   iii. order points along contour
% at each time, for each point:
%   compute integral along Burger circuit: 
%           compute angles, then charge from angles.
                       
% find neighbors
[Idx,~] = rangesearch(r,r,dr);
% uncomment following line for alternative neighbor finding
%[Idx,~] = knnsearch(er,er,'K',7,'IncludeTies',true); 

% compute and store local basis at each point
for i=1:(l_grid)
    r0 = squeeze(r(i,:));
    hood = Idx{i};
    hood = hood(hood ~= i); % remove i from the list
    sh = length(hood);
    nh = zeros(sh, 3);
    th = zeros(sh, 3);
    ks = zeros(sh, 3);
    ts1 = zeros(sh, 3);
    ts2 = zeros(sh, 3);
    % for each neighbor: compute local basis
    for s=1:sh
        nh(s,:) = (r(hood(s),:) - dot(r(hood(s),:),r0) * r0);
        nh(s,:) = nh(s,:) / norm(nh(s,:));
        th(s,:) = cross(nh(s,:),r0);
        th(s,:) = th(s,:)/ norm(th(s,:));
    end

    % order points

    alphas = zeros(sh,1);
    for s=1:sh
        alphas(s) = atan2(dot(nh(s,:), th(1,:)), dot(nh(s,:), nh(1,:)));% + pi; 
    end
    [~,I] = sort(alphas);
    hood = hood(I);
    for s=1:sh
        if s== sh
            s1 = 1;
        else
            s1 = s+1;
        end
        ks(s,:) = cross(r(hood(s),:), r(hood(s1),:));
         ks(s,:) = ks(s,:) / norm(ks(s,:));
         ts1(s,:) = cross(r(hood(s),:), ks(s,:));
         ts1(s,:) = ts1(s,:)/norm(ts1(s,:));
         ts2(s,:) = cross(r(hood(s1),:), ks(s,:));
         ts2(s,:) = ts2(s,:)/norm(ts2(s,:));
    end
    basis(i).hood = hood; 
    basis(i).nh = nh(I,:);
    basis(i).th = th(I,:);
    basis(i).ks = ks;
    basis(i).ts1 = ts1;
    basis(i).ts2 = ts2;   
end


for t=t_arr 
    charges_raw(t).pos = [];
    charges_raw(t).q = [];
    for i=1:(l_grid)
        % compute angles, integral Q
        sh = length(basis(i).hood);
        dtheta = zeros(sh,1);
        for s=1:sh % compute angle of vector field with a set direction ( angle with vector from center to first neighbor on the list)
             if s == sh
                 s1 = 1;
             else
                 s1 = s + 1;
             end
             ks = basis(i).ks(s,:);
             ts1 = basis(i).ts1(s,:);
             ts2 = basis(i).ts2(s,:);
             
             
             vt1 = vx(basis(i).hood(s),t)* ts1(1) +...
                                vy(basis(i).hood(s),t)* ts1(2) +...
                                     vz(basis(i).hood(s),t) * ts1(3);
             vt2 = vx(basis(i).hood(s1),t)* ts2(1) +...
                                vy(basis(i).hood(s1),t)* ts2(2) +...
                                     vz(basis(i).hood(s1),t) * ts2(3);
             
             vk1 = vx(basis(i).hood(s),t)* ks(1) +...
                                vy(basis(i).hood(s),t)* ks(2) +...
                                     vz(basis(i).hood(s),t) * ks(3);
             vk2 = vx(basis(i).hood(s1),t)* ks(1) +...
                                vy(basis(i).hood(s1),t)* ks(2) +...
                                     vz(basis(i).hood(s1),t) * ks(3);
             
             theta1 = atan2(vk1/sqrt(vt1^2+vk1^2),...
                            vt1/sqrt(vt1^2+vk1^2));
                        
             theta2 = atan2(vk2/sqrt(vt2^2+vk2^2),...
                            vt2/sqrt(vt2^2+vk2^2));
             
            dtheta(s) = mod(theta2 - theta1 +pi, 2*pi) - pi;  % theta(s+1) - theta(s) \in [-pi, pi]
        end
        
        Q = sum(dtheta)/(2*pi); % topological charge

        if abs(Q) > .05 % threshold above which a point is recorded as a defect
            
            charges_raw(t).q = [charges_raw(t).q round(Q)];
            charges_raw(t).pos = [charges_raw(t).pos i];
        end
    end
end

% %%% debug plot
% numcharge = zeros(1, length(t_arr));
% for t=t_arr
%     numcharge(t) = length(charges_raw(t).q);
% end
% figure
% plot(2*t_arr, numcharge)
% xlabel("Time (min)")
% ylabel("number of charges")

%%%% Reduce clustered sets of defect to single points (important if defects
%%%% not aligned with the sampling grid)

for t=t_arr
    charges(t).q = [];
    charges(t).pos = [];
    
    % find clusters of points
    idx_raw = charges_raw(t).pos;
    pos = r(idx_raw,:);
    T = clusterdata(pos,'Criterion','distance','Cutoff',dr);
    
    n = length(unique(T));
    for i = 1:n
        % find first index k such as T(k) = n
        k = find(T == i,1);
        charges(t).pos = [charges(t).pos idx_raw(k)];
        charges(t).q = [charges(t).q charges_raw(t).q(k)];
    end
end
toc
    



end