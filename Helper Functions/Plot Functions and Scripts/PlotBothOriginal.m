% chdir("C:\Users\asteinhardt\Desktop\upscaling_Data_Highway_Bloom\Both");%% aos 6/14/23
a=readtable("Both_10.10.10.66_"+num2str(1716)+".csv");opts = detectImportOptions("Both_10.10.10.66_"+num2str(1716)+".csv");
    hummaq=["range_m_","intensity","laserRow","laserCol"] ;opts.SelectedVariableNames =hummaq;    
RAD=int16(40);thresh=int16(80);
for i=1716:1724;X{i-1715}=table2array(readtable("Both_10.10.10.66_"+num2str(i)+".csv",opts)); %#ok<*SAGROW>convert to integ 0-65535
 Xint=X{i-1715}(:,[3,4,1,2]);k=find(diff(Xint(:,1))>0,1,'first');%flicker bug removal
 if(~isempty(k)); Xint(k+1:end,:)=[];end%mapped to uniform grid
 Xint(:,3)= (min(256,Xint(:,3))*64);X{i-1715}=Xint;%range now unit16, int unit l, 6 echo, 
end%% we now have for each return up to 4 entries reading data,  curat first tick is higher it is building plumbingdone input read
%%***********************************************************************
frame=cell(0);
for j=1:length(X);x= cell(0) ;k= (unique(X{j}(:,1)));
    kaz=(unique(X{j}(:,2)));XX= (X{j});
  sizes=zeros(size(k,1),size(kaz,1));%% we can use upscale delta in future to determine which entries to upscale for now we do all
   for iq=1:length(k)%baseline shots on   grid redundant ops can accelerate obvious
   for iqq=1:length(kaz)
temp=XX(XX(:,1)==k(iq)&XX(:,2)==kaz(iqq),[3,4]);
       [~,iii]=sort(temp(:,1));temp=temp(iii,:);  %example init    23754 prededup   23694 after dedup 2% not much
   if(~isempty(temp))
       try
      hummaQ=     [temp,(temp(:,1)*0+1)*[(k(iq)-16*40)/40,(kaz(iqq)-32*40)/40]];
       catch
           keyboard
       end
       hummaQ(1)=hummaQ(1)/64;hummaQ(2)=min(1,hummaQ(2)/6000);
       x{iq,iqq}=hummaQ;
     sizes(iq,iqq)=size(x{iq,iqq},1);
   end
   end
   end
disp(['on frame number',num2str(j)])
   frame{j}=cell2mat(x(:));
end
%% frame{j} is frame with all echos, in range intensity el az (deg)
%% convert to x y z:
for i=1:length(frame)
Q=frame{j};
%%curate kill all >200m
Q(Q(:,1)>200,:)=[];
r=Q(:,1);fee=Q(:,3);thet=Q(:,4);
x=r.*cosd(fee).*cosd(thet);y=r.*cosd(fee).*sind(thet);z=r.*sind(fee);
I=Q(:,2);I=min(1000,I)/100;
scatter3(x,(z),-(y),10,I,'filled'); colorbar; ylim([-10 10]); zlim([-10 5])
keyboard
end