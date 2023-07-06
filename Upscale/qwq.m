chdir("C:\Users\asteinhardt\Desktop\upscaling_Data_Highway_Bloom\Both");%% aos 6/14/23
a=readtable("Both_10.10.10.66_"+num2str(1716)+".csv");opts = detectImportOptions("Both_10.10.10.66_"+num2str(1716)+".csv");
    hummaq=["range_m_","intensity","laserRow","laserCol"] ;opts.SelectedVariableNames =hummaq;    
RAD=int16(40);thresh=int16(80);
for i=1716:1724;X{i-1715}=table2array(readtable("Both_10.10.10.66_"+num2str(i)+".csv",opts)); %#ok<*SAGROW>convert to integ 0-65535
 Xint=X{i-1715}(:,[3,4,1,2]);k=find(diff(Xint(:,1))>0,1,'first');%flicker bug removal
 if(~isempty(k)); Xint(k+1:end,:)=[];end%mapped to uniform grid
 Xint(:,3)=int16(min(256,Xint(:,3))*64);X{i-1715}=Xint;%range now unit16, int unit l, 6 echo, 
end%% we now have for each return up to 4 entries reading data,  curat first tick is higher it is building plumbingdone input read
%%***********************************************************************

for j=1:length(X);x= cell(0) ;k=int16(unique(X{j}(:,1)));kaz=int16(unique(X{j}(:,2)));XX=int16(X{j});
    humma=uint8(0);clumps=cell(0);sizes=zeros(size(k,1),size(kaz,1));%% we can use upscale delta in future to determine which entries to upscale for now we do all
   for iq=1:length(k)%baseline shots on   grid redundant ops can accelerate obvious
   for iqq=1:length(kaz)
temp=XX(XX(:,1)==k(iq)&XX(:,2)==kaz(iqq),[3,4]);
       [~,iii]=sort(temp(:,1));temp=temp(iii,:);
     %example init    23754 prededup   23694 after dedup 2% not much
     %% dedup!
      if(size(temp,1)==2&&abs(diff(temp(1:2,1)))<RAD);temp=mean(temp,1);
      elseif(size(temp,1)==3); kkk=find(abs(diff(temp(:,1)))<RAD);
          if(~isempty(kkk));kkk=unique([kkk,kkk+1]);%% options 1,2,3  1,2  2,3
           if(length(kkk)==2);temp=[temp(setdiff([1:3],kkk),:);mean(temp(kkk,:),1)];
           elseif(length(kkk)==3);temp=mean(temp,1);
           end
          end
      elseif(size(temp,1)==4)%% options size kk 2,3,4 
          kkk=find(abs(diff(temp(:,1)))<RAD);
          if(~isempty(kkk));kkk=unique([kkk,kkk+1]);
              if(length(kkk)==4);temp=mean(temp,1);%done length4
              elseif(length(kkk)==3)
                   temp=[temp(setdiff([1:4],kkk),:);mean(temp(kkk,:),1)];
              elseif(length(kkk)==2)% options  1,2 | 2,3 | 3,4 
  temp=[mean(temp(setdiff([1:4],kkk),:),1);mean(temp(kkk,:),1)];
              end
          end
      end
              %% end dedup 6202023
     x{iq,iqq}=temp;
     % if(~isempty(temp));keyboard;end
     sizes(iq,iqq)=size(x{iq,iqq},1);
 if(sizes(iq,iqq)>1);[~,iii]=sort(-x{iq,iqq}(:,1));
     x{iq,iqq}=x{iq,iqq}(iii,:)';
     x{iq,iqq}=reshape(x{iq,iqq},1,[]); 
 end
    end
   end
 
 %% end build of x this is all in fpga except perhaps dedup (this could be there also)

 [shotlist,sizes,rowses]=upscale(RAD,x,k,kaz,thresh);
  [shotlist2,sizes2]=            upscale2(x,k,kaz,RAD,thresh);
end
%%%%%%%%%%%%%%%%% 23694 baseline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shotlist,sizes,rowses]=upscale(RAD,x,k,kaz,thresh)
%core upscale prep max output 2iq iqq = 2lengthk-2 etc unless+1 for diag
shotlist=cell(length(k)*2-1,2*length(kaz)-1);
sizes=zeros(length(k)*2-1,2*length(kaz)-1);rowses=sizes;%see below sizes
for iq=1:length(k)-1%uses all input data 1:k 1:kaz , output up to 2k-1 2kaz-1
   for iqq=1:length(kaz)-1
                       if(~isempty(x{iq,iqq}))
                           c=x{iq,iqq};
                           xtempy=c(1,:);
                           c(1,:)=[];
          while(~isempty(c));xtempy=[xtempy,c(1,:)];c(1,:)=[];end %#ok<AGROW>
          % if(length(xtempy)>4);xtempy=xtempy(:,[1:4]);
          % end
          shotlist{iq*2-1,-1+2*iqq}=int16(xtempy);
          n=length(shotlist{iq*2-1,-1+2*iqq});
          S=shotlist{iq*2-1,-1+2*iqq};
          sizes(iq*2-1,-1+2*iqq)=n;
                       else;S=[];
                       end%gen base shot above

                 if(~isempty(x{iq,iqq+1})&&~isempty(S));c=x{iq,iqq+1};xtempy=c(1,:);c(1,:)=[];
          while(~isempty(c));xtempy=[xtempy,c(1,:)];c(1,:)=[];end %#ok<AGROW>
        m=length(xtempy);
           if(n<=m);qaz=[S,[xtempy(:,1:n)]];%kills closer echos
                 elseif(m<n);qaz=[[S(:,1:m)],[xtempy]];
           end;shotlist{iq*2-1,2*iqq}=int16(qaz);sizes(iq*2-1,2*iqq)=length(qaz);
                 end%end az shot

            if(~isempty(x{iq+1,iqq})&&~isempty(S));c=x{iq+1,iqq};xtempy=c(1,:);c(1,:)=[]; %remove closest if 4 shots etc usually ase
                        while(~isempty(c));xtempy=[xtempy,c(1,:)];c(1,:)=[];end%#ok<AGROW>
        m=length(xtempy);
                    if(n<=m);qaz=[S,[xtempy(:,1:n)]];%kills closer echos
                 elseif(m<n);qaz=[[S(:,1:m)],[xtempy]];
                     end;shotlist{iq*2,2*iqq-1}=int16(qaz);sizes(iq*2,2*iqq-1)=length(qaz);
            end%end el
   end
end; RANK=int16(zeros(length(k)*2-1,2*length(kaz)-1));
for i=1:length(k)-1
    for jj=1:length(kaz)-1
        RANK(i*2,2*jj-1)=int16(length(cell2mat(shotlist(2*i,2*jj-1))));
        RANK(2*i-1,2*jj)=int16(length(cell2mat(shotlist(2*i-1,2*jj))));
        % if(isempty(shotlist{i,jj}));shotlist{i,jj}=int16(zeros(0,0));end
    end
end
%%assumes 4 or 8 az or el now timing in 16bit 
t4=(RANK==4);t8=(RANK==8);
SHOT4=cell2mat(shotlist(t4));
SHOT8=cell2mat(shotlist(t8)); 
tic
t1=abs(SHOT4(:,1)-SHOT4(:,3))<RAD;SHOT4(t1,[1,2])=(SHOT4(t1,[1,2])+SHOT4(t1,[3,4]))/2;
t2=abs(SHOT8(:,1)-SHOT8(:,5))<RAD;t3=abs(SHOT8(:,3)-SHOT8(:,7))<RAD; %2 3 on 8
tboth=t2&t3;t2loner=t2&~tboth;t3loner=t3&~tboth;
SHOT8(t2loner,[1,2])=(SHOT8(t2loner,[1,2])+SHOT8(t2loner,[5,6]))/2;
SHOT8(t3loner,[1,2])=(SHOT8(t3loner,[7,8])+SHOT8(t3loner,[3,4]))/2;
SHOT8(tboth,1:4)=(SHOT8(tboth,[1:4])+SHOT8(tboth,[5:8]))/2;
toc
SHOT4=SHOT4(:,1:2);SHOT8=SHOT8(:,1:4);
m1=length(SHOT4(t1,:));
t4(t4==1)=t1;t8copy=t8;t8copy(t8copy==1)=tboth;t8(t8==1)=t2loner|t3loner;%perpsall!
shotlist(t4)=mat2cell(SHOT4(t1,:),ones(m1,1),2);
m2=size(SHOT8(tboth,1),1);m3=length(SHOT8(t2loner|t3loner,:));
shotlist(t8copy)=mat2cell(SHOT8(tboth,1:4),ones(m2,1),4);
shotlist(t8)=mat2cell(SHOT8(t2loner|t3loner,1:2),ones(m3,1),2);%% below pointers
z=zeros(length(k)*2-1,2*length(kaz)-1);z([2:2:length(k)*2-1]',[2:2:2*length(kaz)-1]')=1;zdi=logical(z);
  z=zeros(length(k)*2-1,2*length(kaz)-1);z([1:2:length(k)*2-1]',[1:2:length(kaz)*2-1]')=1;zbase=logical(z);
s=logical(~(t8|t8copy|t4|z|zdi));tau=sum(s(:)); hummay=cell(tau,1);
hummay(:)=mat2cell(zeros(tau,0),ones(tau,1),0);  
  shotlist(logical(~(t8|t8copy|t4|z|zdi)))=hummay;
  %% sizes 

for i=1:length(k)*2-1
    for jj=1:length(kaz)*2-1
sizes(i,jj)=length(shotlist{i,jj});
if sizes(i,jj)>8;keyboard;end
    end
end
 
%% ready for diag note +1 so up to 2k-2+1=2k-1
for iq=1:length(k)-1
  for iqq=1:length(kaz)-1
      if(~isempty(x{iq+1,iqq+1})||~isempty(shotlist{iq*2+1,2*iqq})||~isempty(shotlist{2*iq,2*iqq+1}))
        qazyy=[f(x{iq+1,iqq+1}) ;f(shotlist{iq*2+1,2*iqq});f(shotlist{2*iq,2*iqq+1})];
        else;qazyy=[];
      end
        if(~isempty(qazyy))
        [~,iii]=sort(-qazyy(:,1));qazyy=qazyy(iii,:);kkk=find(abs(diff(qazyy(:,1)))<RAD);
        kkkk=find(diff(kkk)>1,2,'first');
        kkk=unique([kkk,kkk+1]);
            if(~isempty(kkk)&&isempty(kkkk))
                qazyy=int16(mean(qazyy(kkk,:),1));
        elseif(~isempty(kkk)&&~isempty(kkkk)&&length(kkkk)==1)
            qazyy=[int16(mean(qazyy(kkk(1:kkkk+1),:),1));...
                int16(mean(qazyy(kkk(kkkk+2:end),:),1))];
                elseif(~isempty(kkk)&&~isempty(kkkk)&&length(kkkk)>=2) 
            keyboard
        elseif(isempty(kkk))
            qazyy=[];
            end
            shotlist{2*iq,2*iqq}=qazyy;
        end
 end
end%now thresh
      for i=1:length(k)
          for jj=1:length(kaz)
Q=shotlist{i*2-1,2*jj-1};
          if(~isempty(Q))
          kkk=find(Q(2:2:end)>thresh);
          pointer=[2:2:length(Q)];
          shotlist{i*2-1,2*jj-1}=[Q(pointer(kkk)-1),Q(pointer(kkk))];
          end
          end
       end
   %regenerate size from scratch to ensure error free
   sizes=sizes*0;
      for i=1:size(sizes,1)
       for jj=1:size(sizes,2)
           if(iscell(shotlist{i,jj}));shotlist{i,jj}=cell2mat(shotlist{i,jj});end
           sizes(i,jj)=numel(shotlist{i,jj});
       end
   % di is 0,3 6 9 or 12 shots
      end%dunndun
 
      disp([' shots size net diag then base then ~either ',num2str(sum(sizes(zdi(:)))/2), ' ',num2str(sum(sizes(zbase(:)))/2), ' ',num2str(sum(sizes(~zbase(:)&~zdi(:)))/2)]);
 t4samples=sizes==4;
 t2samples=sizes==2;
 disp([' 4 and 2 upscale ',num2str(sum(t2samples(:))),' ',num2str(2*sum(t4samples(:)))])
end       
%%%% 76979 shots end started with 23064 above thresh 656>2echo = .85%

function [shotlist2,sizes]=upscale2(x,k,kaz,RAD,thresh)
shotlist2=cell(length(k)*2-1,length(kaz)*2-1);sizes=zeros(length(k)*2-1,length(kaz)*2-1);
%% done with prep now see if we can improve az el or if further tuning is needed this is optimized module
z=zeros(length(k)*2-1,2*length(kaz)-1);z([2:2:length(k)*2-1]',[2:2:2*length(kaz)-1]')=1;
zdi=logical(z);
  z=zeros(length(k)*2-1,2*length(kaz)-1);z([1:2:length(k)*2-1]',[1:2:length(kaz)*2-1]')=1;zbase=logical(z);
%% main loop
for iq=1:length(k)-1
 for iqq=1:length(kaz)-1
   if(~isempty(x{iq,iqq}))
                           c=x{iq,iqq};
                           xtempy=c(1,:);
                           c(1,:)=[];
          while(~isempty(c));xtempy=[xtempy,c(1,:)];c(1,:)=[];end %#ok<AGROW>
          shotlist2{iq*2-1,-1+2*iqq}=int16(xtempy);n=length(shotlist2{iq*2-1,-1+2*iqq});
          S=shotlist2{iq*2-1,-1+2*iqq};sizes(iq*2-1,-1+2*iqq)=n;
                       else;S=[];
   end%gen base shot above
   i=iq;jj=iqq;%isomorphic
       if(~isempty(x{i,jj+1})&&~isempty(S))
      c=x{i,jj+1};shotlist2{i*2-1,2*jj}=sordid(S,c);
          range=shotlist2{i*2-1,2*jj}(1:2:end);
          % if(length(range)>2);keyboard;end
          int=shotlist2{i*2-1,2*jj}(2:2:end);
          kk=find(abs(diff(range))<RAD);kkk=find(diff(kk)>1); 
          if(~isempty(kk))
                      if(isempty(kkk))
                          shotlist2{i*2-1,2*jj}=int16([sum(range([kk,kk(end)+1]))/(1+length(kk)),sum(int([kk,kk(end)+1]))/(1+length(kk))]);
                      elseif(length(kkk)>=1)
              shotlist2{i*2-1,2*jj}=int16([sum(range([kk(1:kkk),kk(kkk)+1]))/(1+length(kk(1:kkk))),...
                  sum(int([kk,kk(end)+1]))/(1+length(kk(1:kkk))),...
             sum(range(kk(kkk)+1:kk(end)+1))/(1+length(kk(kkk)+1:kk(end)+1)),...
             sum(int([kk(kkk)+1:kk(end)+1]))/(1+length(kk(kkk)+1:kk(end)+1))]);
                      end
               KKK=length(shotlist2{i*2-1,2*jj});sizes(i*2-1,2*jj)=KKK;
          end
          end%2below is el optimal see upscale for fast
      %endif
    
      if(~isempty(x{i+1,jj})&&~isempty(S));c=x{i+1,jj};shotlist2{i*2,2*jj-1}=[sordid(S,c)];
          range=shotlist2{i*2,2*jj-1}(1:2:end);int=shotlist2{i*2,2*jj-1}(2:2:end);
          kk=find(abs(diff(range))<RAD);kkk=find(diff(kk)>1); 
          if(~isempty(kk))
                      if(isempty(kkk))
 shotlist2{i*2,2*jj-1}=int16([sum(range([kk,kk(end)+1]))/(1+length(kk)),sum(int([kk,kk(end)+1]))/(1+length(kk))]);
                      elseif(length(kkk)>=1)
              shotlist2{i*2,2*jj-1}=int16([sum(range([kk(1:kkk),kk(kkk)+1]))/(1+length(kk(1:kkk))),...
                  sum(int([kk,kk(end)+1]))/(1+length(kk(1:kkk))),...
             sum(range(kk(kkk)+1:kk(end)+1))/(1+length(kk(kkk)+1:kk(end)+1)),...
             sum(int([kk(kkk)+1:kk(end)+1]))/(1+length(kk(kkk)+1:kk(end)+1))]);
                      end
              sizes(i*2,2*jj-1)=length(shotlist2{i*2,2*jj-1});
          end
  
      end%doneel    %azeldun               
 end
end

%% now we do the diagonal prep and range coincidence diagonal completion
for iq=1:length(k)-1
  for iqq=1:length(kaz)-1
      if(~isempty(x{iq+1,iqq+1})||~isempty(shotlist2{iq*2,2*iqq+1})||~isempty(shotlist2{2*iq+1,2*iqq}))
        qazyy=[f(x{iq+1,iqq+1}) ;f(shotlist2{iq*2,2*iqq+1});f(shotlist2{2*iq+1,2*iqq})];
        else;qazyy=[];
      end
        if(~isempty(qazyy))
        [~,iii]=sort(-qazyy(:,1));qazyy=qazyy(iii,:);kkk=find(abs(diff(qazyy(:,1)))<RAD);
        kkkk=find(diff(kkk)>1,1,'first');kkk=unique([kkk,kkk+1]);
            if(~isempty(kkk)&&isempty(kkkk));qazyy=int16(mean(qazyy(kkk,:),1));
        elseif(~isempty(kkk)&&~isempty(kkkk));qazyy=int16(mean(qazyy(kkk(1:kkkk+1),:),1));  
        elseif(isempty(kkk));qazyy=[];
            end;shotlist2{2*iq,2*iqq}=qazyy;sizes(2*iq,2*iqq)=length(qazyy);
        end
 end
end
 %% now we do the baseline thresholding on intensity
for iqq=1:length(k)
    for jj=1:length(kaz)
        Q=shotlist2{2*iqq-1,2*jj-1};
        if(~isempty(Q))
             kkk=find(Q(2:2:end)>thresh);pointer=[2:2:length(Q)];
          shotlist2{iqq*2-1,2*jj-1}=[Q(pointer(kkk)-1),Q(pointer(kkk))];
          sizes(iqq*2-1,-1+2*jj)=length(Q);
        end
    end
end
%% done
      disp([' shots size net diag then base then ~either ',num2str(sum(sizes(zdi(:)))/2), ' ',num2str(sum(sizes(zbase(:)))/2), ' ',num2str(sum(sizes(~zbase(:)&~zdi(:)))/2)]);
 t4samples=sizes==4;
 t2samples=sizes==2;
 disp([' 4 and 2 upscale ',num2str(sum(t2samples(:))),' ',num2str(2*sum(t4samples(:)))])
end
%%%%% 261 >2 echos
function y=sordid(xdat,ydat)
top=[]; %#ok<NASGU>
if(~isempty(xdat));n=size(xdat,2);else;n=0;end
if(~isempty(ydat));nn=size(ydat,2);else;nn=0;end
 if(n==4);xdat=[xdat(1:2);xdat(3:4)];end
  if(n>=8);xdat=[xdat(1:2);xdat(3:4);xdat(5:6);xdat(7:8)];end
     if(nn==4);ydat=[ydat(1:2);ydat(3:4)];end
      if(nn==6);ydat=[ydat(1:2);ydat(3:4);ydat(5:6)];end
           if(n==6);xdat=[xdat(1:2);xdat(3:4);xdat(5:6)];end
  if(nn>=8);ydat=[ydat(1:2);ydat(3:4);ydat(5:6);ydat(7:8)];end
y=[xdat;ydat]'; 
% if(length(y)>2);keyboard;end
[~,iii]=sort(y(1,:));y=y(:,iii);
y=y(:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%
function y=f(x)
if(~isempty(x))
[~,m]=size(x);if(m>8);x=x(:,1:8);m=8;end
if(m==2)
      y=x(:,1:2);
elseif(m==4)
    y=[x(:,1:2);x(:,3:4)];
elseif(m==6)
    y=[x(:,1:2);x(:,3:4);x(5:6)];
elseif(m==8)
        y=[x(:,1:2);x(:,3:4);x(5:6);x(7:8)];
elseif(m==1)
    y=x;
end
else
    y=[];
end
end