function [ in_box, far_field ] = FindPointsInBox( data, box_center, box_width )
    
    in_box = [];
    far_field = [];
    
    % so I don't have to take square roots later 
    in_box_dist = box_width * box_width;
    far_dist = 4 * in_box_dist;
    
    for i = 1:size(data,2)
       
        point = data(:,i);
        
        box_center_dist = (box_center - point)' * (box_center - point);
        
        if (box_center_dist < in_box_dist)
            
            in_box(end+1) = i;
            
        elseif(box_center_dist > far_dist) 
            
            far_field(end+1) = i;
            
        end
        
        
    end
    
end

