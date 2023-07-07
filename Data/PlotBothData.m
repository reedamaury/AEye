%chdir("C:\Users\asteinhardt\Desktop\upscaling_Data_Highway_Bloom\Both");%% aos 6/14/23
a=readtable("Both_10.10.10.66_"+num2str(1716)+".csv");opts = detectImportOptions("Both_10.10.10.66_"+num2str(1716)+".csv");
hummaq=["range_m_","intensity","laserRow","laserCol"] ;opts.SelectedVariableNames =hummaq;    

RAD=int16(40);thresh=int16(80);
for i=1716:1724;
    % taking selected variables from each excel file and putting it into the cell array X
    X{i-1715}=table2array(readtable("Both_10.10.10.66_"+num2str(i)+".csv",opts)); %#ok<*SAGROW>convert to integ 0-65535
    
    % creating a matrix from current cell array (note that the columns have
    % been rearranged)
    Xint=X{i-1715}(:,[3,4,1,2]);

    % if a laser in a successive row fires while current row is firing,
    % remove it. (i.e. [880 880 880 874 880 880 880] remove 874) corrects
    % laser misfires. diff does backwards differentiation, starts from last
    % element of vector and continues to first element. [0 0 -6 6 0 0]
    kk=find(diff(Xint(:,1))>0,1,'first'); % flicker bug removal
    if(~isempty(kk)); 
        Xint(kk+1:end,:)=[]; % when you get to new row of lasers, remove the data from the first laser in that row
    end % mapped to uniform grid
    Xint(:,3)= (min(256,Xint(:,3))*64);
    X{i-1715}=Xint; %range now unit16, int unit l, 6 echo, 
end %% we now have for each return up to 4 entries reading data,  curate first tick is higher it is building plumbing done input read
%%***********************************************************************
frame=cell(0);
for j=1:length(X); % for each frame
    xx= cell(0) ;
    k= (unique(X{j}(:,1))); % laser rows (no repetitions)
    kaz=(unique(X{j}(:,2))); % laser columns (no repetitions)
    XX= (X{j}); % XX gets the current frame of data
    sizes=zeros(size(k,1),size(kaz,1)); % zeros for all the rows and columns of the lasers %% we can use upscale delta in future to determine which entries to upscale for now we do all
    for iq=1:length(k) % length of laser rows  % baseline shots on grid redundant ops can accelerate obvious
        for iqq=1:length(kaz) % length of laser columns 
            temp=XX(XX(:,1)==k(iq) & XX(:,2)==kaz(iqq),[3,4]); % the indices of XX where all the rows equal the current row k(iq) (first value will be 484) and columns equal the current column k(iqq) first value will be 29,
                                                                            % get the range and intensity values of the current row k(iq) and column k(iqq) and assign to temp
            [~,iii]=sort(temp(:,1)); % sort range values in ascending (increasing) order and assign indices of ascending order temp to iii
            temp=temp(iii,:);  % sorts temp so all range values (with corresponding intensities) are in ascending order for the  % example init    23754 prededup   23694 after dedup 2% not much
            % so for each laser row and laser column (i.e. each laser shot)
            % sort the range in ascending order with corresponding
            % intensity values...why might there be multiple
            % range-intensity values for each laser? because a single
            % shot/pulse can reflect off of multiple surfaces/objects at
            % varying distances. These are echoes

            if(~isempty(temp))
               try
                    hummaQ =     [temp,(temp(:,1)*0+1)*[(k(iq)-16*40)/40,(kaz(iqq)-32*40)/40]]; % columns 1 and 2 are range and intensity, 3 and 4 are elevation and azimuth. 1,3, and 4 give us spherical coordinates.
               catch
                    keyboard
               end
               hummaQ(:,1)=hummaQ(:,1)/64; % divide range values by 64 and assign to all the rows of first column of hummaQ
               hummaQ(:,2)=min(1,hummaQ(:,2)/6000); % divide intensity values by 6000, if any resulting intensities are greater than 1, then intensity just gets 1. assign this to all the rows of second column of hummaQ %hummaQ(1)=hummaQ(1)/64; hummaQ(2)=min(1,hummaQ(2)/6000);
               xx{iq,iqq}=hummaQ; % assign range, intensity, row, column to cell array
               sizes(iq,iqq)=size(xx{iq,iqq},1); % lets us know how many echos each shot has sense the # of rows will be equal to number of returns for that shot
            end
        end
     end
    disp(['on frame number',num2str(j)])
    frame{j}=cell2mat(xx(:));
end
%% frame{j} is frame with all echos, in range intensity el az (deg)
%% convert to x y z:
for i=1:length(frame)
    Q=frame{j};
    %%curate kill all >200m
    Q(Q(:,1)>200,:)=[]; % if range is greater than 200 meters, remove data (all columns)
    r=Q(:,1);
    fee=Q(:,3);
    thet=Q(:,4);
    x=r.*cosd(fee).*cosd(thet);
    y=r.*cosd(fee).*sind(thet);
    z=r.*sind(fee);
    I=Q(:,2);I=min(1000,I)/100;
    scatter3(x,y,-z,10,I,'filled'); colorbar; ylim([-10 10]); zlim([-2 5]); view(-90,35); xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)")
    keyboard
end