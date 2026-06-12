function update_position(obj,currsub)
        
    pos = get(currsub, 'Position');
    new_pos = pos + [-0.05 -0.05 0.07 0.07];
    set(currsub, 'Position', new_pos)

end