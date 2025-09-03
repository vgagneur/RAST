function  [s100,s110,s112,s123]=RAST_euler_to_slip(euler,crystal)
%Plane types + Euler to slip trace caculation, original version was written in python by Johann Pauli Magnussen (co-author). 
% V. Gagneur translated it to matlab and merged it with the other functions of the code.

if crystal==1

    %Slip Planes
    %The crystal structures of interest are A2 (BCC) and B2 (BCC with the central atom replaced with different element). The slip directions and planes for the A2 and B2 systems can be different. Half the entries have been commented out because of symmetry.
    %
    % Define set of slip planes for A2 (BCC) and B2 (BCC with central atom replaced with different element).
    % A2 may slip on {110}, {112}, and {123} planes. NB some literature
    % evidence of <100> {110} in bcc chromium has also been reported
    % https://doi.org/10.1016/j.actamat.2021.117485.
    
    planes_110={[1 1 0], [-1 1 0], [1 0 1], [1 0 -1], [0 1 1], [0 -1 1]};
    planes_112={[1 1 2], [-1 1 2], [1 -1 2], [1 1 -2], [1 2 1], [-1 2 1], [1 -2 1], [1 2 -1], [2 1 1], [-2 1 1], [2 -1 1], [ 2 1 -1]};
    
    planes_123={[1 2 3], [-1 2 3], [1 -2 3], [1 2 -3], [1 3 2], [-1 3 2], [1 -3 2], [1 3 -2]...
        , [2 1 3], [-2 1 3], [2 -1 3], [2 1 -3], [2 3 1], [-2 3 1], [2 -3 1], [2 3 -1]...
        , [3 2 1], [-3 2 1], [3 -2 1], [3 2 -1], [3 1 2], [-3 1 2], [3 -1 2], [3 1 -2]};
    
    %B2 can slip on all above, as well as {100} planes.
    
    planes_100={[ 1 0 0], [ 0 1 0], [ 0 0 1]};

elseif crystal==2

    % each plane type will be named with the corresponding bcc plane type in
    % the output (e.g. (111) match will be called "100" match.
    % fcc plane types {[1 1 1],[-1 1 1],[1 -1 1],[1 1 -1]}
    % planes are entered twitce  because code crashes if less than two
    % planes defined, but doesnt change anything for the analysis and
    % counting procedures

    planes_100={[1 1 1],[1 1 1]};
    planes_110={[-1 1 1],[-1 1 1]}; 
    planes_112={[1 -1 1],[1 -1 1]};
    planes_123={[1 1 -1],[1 1 -1]};

end

euler=str2num(euler{1}); %converts Euler format

%calculates intersections for each plane type
s100 = calc_slip_traces(planes_100,euler); 

s110 = calc_slip_traces(planes_110,euler);

s112 = calc_slip_traces(planes_112,euler);

s123 = calc_slip_traces(planes_123,euler);

return

%Rotating Planes
%Rotations are done by Bunge convention.
%NB convention may vary depending on software used for ebsd acquisition, MTEX functions used... be careful, adapt to your system, and make sure to test a "simple" sample of known slip activity first to check conventions are correct for you.

%Rotation matrix about x-axis, theta in degrees.
function M=rot_x(theta)
    theta = deg2rad(theta);
    M=[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
    return 
end

function M=rot_y(theta)
    theta = deg2rad(theta);
    M=[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
    return 
end

function M=rot_z(theta)
    theta = deg2rad(theta);
    M=[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
    return 
end

function M=rotate_zxz(v,euler)
    v = v;
    om = rot_z(euler(3)) * rot_x(euler(2)) * rot_z(euler(1));
    M=squeeze(v * om);
    return
end
function M=rotate_zxy(v,euler)
    v = v;
    om = rot_y(euler(3)) * rot_x(euler(2)) * rot_z(euler(1)); %fail ????? from the original coede ?
    M=squeeze(v * om);
    return
end


% Calculate slip trace orientation.
function M=calc_slip_trace(slip_plane, euler)
    % x-y plane is the view plane.
    xy_plane = [0 0 1];
    % Slip plane is rotated using orientation matrix (Bunge convention).
    slip_plane_rotated = rotate_zxz(slip_plane,euler);
    % Slip direction is cross product.
    M=cross(slip_plane_rotated,xy_plane);
    return
end

% Calculate orientation for multiple slip traces.
function M=calc_slip_traces(slip_planes, euler)
    slip_dirs = [];
    for x=1:length(slip_planes)
        plane = slip_planes{x};
        new_slip_dir = calc_slip_trace(plane,euler);
        slip_dirs=[slip_dirs;new_slip_dir];
    end
    slip_dirs(:,3)=[];
    
    M=slip_dirs;
    return 
end

end