function [ in_box, far_field ] = FindPointsInBox( data, box_center, box_width, ws_index )
    
    if (nargin == 3)
       ws_index = 2; 
    end

    in_box = [];
    far_field = [];
    
    % so I don't have to take square roots later 
    in_box_dist = box_width * box_width;
    far_dist = ws_index * ws_index * in_box_dist;
    
    box_center_dists = sum(bsxfun(@minus, data, box_center).^2,1);
    
    for i = 1:size(data,2)
       
        if (box_center_dists(i) < in_box_dist)
            
            in_box(end+1) = i;
            
        elseif(box_center_dists(i) > far_dist) 
            
            far_field(end+1) = i;
            
        end
        
        
    end
    
end

