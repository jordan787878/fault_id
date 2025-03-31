function angle_deg = get_angle_deg(v1, v2)
    % Normalize the vectors
    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);
    
    % Calculate the dot product
    dot_product = dot(v1, v2);
    
    % Calculate the angle in radians
    angle_rad = acos(dot_product);
    
    % Convert angle to degrees
    angle_deg = rad2deg(angle_rad);
end