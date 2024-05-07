function [ext_angle] = seg_dist(d1, d2, d3)
    ba = d1 - d2;
    bc = d3 - d2;

    % Initialize array of length d1
    ext_angle = zeros(size(d1, 1), 1);

    % Loop over height of d1
    for i = 1:size(d1, 1)
 
        % Calculate dot product
        dot_product = dot(ba(i,:), bc(i,:), 2);

        % Calculate angle in radians
        angle_rads = dot_product / (norm(ba(i,:)) * norm(bc(i,:)));

        % Calculate inside angle in degrees
        inside_angle_degs = (180 / pi) * acos(angle_rads);

        % Calculate exterior angle
        ext_angle(i, 1) = 180 - inside_angle_degs;
    end
end