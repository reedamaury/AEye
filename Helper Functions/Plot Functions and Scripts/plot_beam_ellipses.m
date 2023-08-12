function plot_beam_ellipses(beam_ellipses)   
    for i = 1:size(beam_ellipses, 1)
        for j = 1:size(beam_ellipses, 2)
            % Retrieve the ellipse parameters
            a = beam_ellipses{i, j}(1);
            b = beam_ellipses{i, j}(2);
            angle = beam_ellipses{i, j}(3);
            x0 = beam_ellipses{i, j}(4);
            y0 = beam_ellipses{i, j}(5);
            
            % Plot the ellipse
            % plot_ellipse(a, b, angle, x0, y0);
            figure;
            t = linspace(0, 2*pi, 100);
            x = a*cos(t);
            y = b*sin(t);
            r = [cos(angle) -sin(angle); sin(angle) cos(angle)];
            xy = r * [x; y];
            plot(xy(1,:) + x0, xy(2,:) + y0);
        end
    end  
end


