function [annulus_positions] = multi_annulus_positions_generator(annulus_radius, annuli_number, destination_folder)
    
    global screen_width
    global screen_height
    
    annulus_positions = zeros(100, 4);
    for a = 1:100
        annulus_center_x = randi(screen_width - annulus_radius * 2) + annulus_radius;
        annulus_center_y = randi(screen_height - annulus_radius * 2) + annulus_radius;
        annulus1_x = annulus_center_x - annulus_radius;
        annulus1_y = annulus_center_y - annulus_radius;
        annulus2_x = annulus_center_x + annulus_radius - 1;
        annulus2_y = annulus_center_y + annulus_radius - 1;
        annulus_positions(a,:) = [annulus1_y, annulus2_y, annulus1_x, annulus2_x];
    end
    
% % % %     sca
% % % %     keyboard
    
% % % %     table_header = ['annulus_p1_y ' 'annulus_p2_y ' 'annulus_p1_x ' 'annulus_p2_x '];
% % % %     dlmwrite([destination_folder 'annulus_positions.txt'], table_header, 'delimiter', '');
% % % %     dlmwrite([destination_folder 'annulus_positions.txt'], annulus_positions, '-append', 'roffset', 1, 'delimiter','\t');

end
