function [frame, xx] = obtainframedata(initialFrame, finalFrame, filenameFormat)

    % initial read for import options
    a = readtable(sprintf(filenameFormat, initialFrame));
    opts = detectImportOptions(sprintf(filenameFormat, initialFrame));
    hummaq = ["range_m_", "intensity", "laserRow", "laserCol"];
    opts.SelectedVariableNames = hummaq;

    RAD = int16(40);
    thresh = int16(80);
    
    for i = initialFrame:finalFrame
        % taking selected variables from each excel file and putting it into the cell array X
        X{i-initialFrame+1} = table2array(readtable(sprintf(filenameFormat, i), opts)); 
        % creating a matrix from current cell array
        Xint = X{i-initialFrame+1}(:, [3,4,1,2]);
        kk = find(diff(Xint(:, 1)) > 0, 1, 'first'); % flicker bug removal
        if (~isempty(kk)); 
            Xint(kk+1:end, :) = [];
        end 
        Xint(:, 3) = (min(256, Xint(:, 3)) * 64);
        X{i-initialFrame+1} = Xint;
    end

    frame = cell(0);
    for j = 1:length(X)
        xx = cell(0);
        k = (unique(X{j}(:, 1)));
        kaz = (unique(X{j}(:, 2)));
        XX = (X{j});
        shot_number = 0;
        sizes = zeros(size(k, 1), size(kaz, 1));
        for iq = 1:length(k)
            for iqq = 1:length(kaz)
                temp = XX(XX(:, 1) == k(iq) & XX(:, 2) == kaz(iqq), [3, 4]);
                [~, iii] = sort(temp(:, 1));
                temp = temp(iii, :);
                if(~isempty(temp))
                    hummaQ = [temp, (temp(:, 1) * 0 + 1) * [(k(iq) - 16 * 40) / 40, (kaz(iqq) - 32 * 40) / 40]];
                    hummaQ(:, 1) = hummaQ(:, 1) / 64;
                    hummaQ(:, 2) = min(1, hummaQ(:, 2) / 6000);
                    shot_number = shot_number + 1;
                    hummaQ(:, 5) = shot_number;
                    xx{iq, iqq} = hummaQ;
                    sizes(iq, iqq) = size(xx{iq, iqq}, 1);
                end
            end
        end
        disp(['on frame number', num2str(j)])
        frame{j} = cell2mat(xx(:));
    end
end
