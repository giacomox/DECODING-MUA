% ------- to re-obtain psths_3t -------------------------------------------
%clear all
close all
%load(['/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data/order/20111005 ch5_R_dt.mat']) %matrix of rasters for the 80 cells obtained by "inport_seq2013_2.m" (to inport from Elphy files)
addpath '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/SCRIPTS'
clear R
load('R')
load('index')

%
disp=1; % display the psths and the cos fits?
fit_psths=0; % do u want to find the psth baricentrum fitting them with a gaussian?
bin=10; % psth bin
ord=1:12; ; % in Elphy condition1 of dt exp is the bar comeing from 90 deg
cv=(0:30:330); %conversion vector
 ff1=(0:pi/6:2*pi-pi/6)*-1;
%
% inport index and R
ct2=hist(index,1:13);
R_dt=R;


for cell=1:96
    clf
    clear fitt ff limm pik
    
   
       
         %----- find the thrashold -------------
         psths_dt(:,13,cell)=zeros(250,1);
 

             for ep=1:ct2(13)
                 psths_dt2(:,ep,13,cell)=hist(R_dt(:,ep,13,cell),1:bin:2500);
                 psths_dt(:,13,cell)=psths_dt(:,13,cell)+psths_dt2(:,ep,13,cell)*100;
             end
             
                  psths_dtS(:,13,cell)=smoothfredo(psths_dt(:,13,cell)/ct2(13),5); % cond 13 is the control condition, blanc
                  mm= mean(psths_dtS(40:150,13,cell));
                  ss= std(psths_dtS(40:150,13,cell));
                  thr=mm+4*ss; % upper thrashold for significant response in the psths
                  thrd=mm-1*ss; % lower thrashold - inhibition
%                   
      %------ dt psths ----------------------------
         for cond=1:12
         
             psths_dt(:,cond,cell)=zeros(250,1);


             for ep=1:ct2(cond)
                 psths_dt2(:,ep,cond,cell)=hist(R_dt(:,ep,cond,cell),1:bin:2500);
                 psths_dt(:,cond,cell)=psths_dt(:,cond,cell)+psths_dt2(:,ep,cond,cell)*100;
             end
            
             psths_dt(1,cond,cell)=0;
             psths_dt(:,cond,cell)=psths_dt(:,cond,cell)/ct2(cond);
             prs(:,cond,cell)=smoothfredo(psths_dt(:,cond,cell),10); % I want to smooth a lot the psths to find the baricentrum of the response
             psths_dtS(:,cond,cell)=smoothfredo(psths_dt(:,cond,cell),5);
    


      %------- find the psths picks ----------------------------      
        
             [W1 pick2]=max(prs(40:150,cond,cell));  % I get the picks only between 400 and 1500 ms (400 is when the bar appear)
             
             [WS1 S]=min(prs(40:150,cond,cell));
            
           
             pick(cond,cell)=pick2*10+400;   
             W(cond,cell)=W1;
                
             pickS(cond,cell)=S*10+400;
             WS(cond,cell)=WS1;
            
%            
       %---- fit psths with a gaussian --------------------------------------------
       x=1:10:2500;
           if fit_psths
                
                if W(cond,cell)>thr
                   aff=1;
                  
                   pp0=[pick(cond,cell), 1000,exp(W(cond,cell))];

                   LB=[400 1 W(cond,cell)-2];
                   UB=[1500 200 W(cond,cell)+1];
                   
                   options=optimset('MaxFunEvals',50000000);
                   
                    prs(1:30,cond,cell)=prs(1:30,cond,cell)./10;
                    prs(150:250,cond,cell)=prs(150:250,cond,cell)./100;
                   
                    psths_dtdd=psths_dt;
                    psths_dtdd(150:250,cond,cell)=0;%psths_dt(150:250,cond,cell)./100;
                    
                   % first fir to find the std
                   pp=fminsearchbnd(@(pp) norm_fitG(pp,prs(:,cond,cell),prs(:,cond,cell),x,aff),pp0,LB,UB);
                   %psths_dt
                   % second fit usign the exp of the psth to find the mean
                   pp1=fminsearchbnd(@(pp) norm_fitG2(pp,prs(:,cond,cell),prs(:,cond,cell),x,aff),pp,LB,UB,options);
                   
                   
%                    pp20=[pp(2) pp(3)];
%                     LB=[10 pp(3)-5];
%                     UB=[1000 pp(3)+5];
%                     pp2=fminsearchbnd(@(pp2) norm_fitG2(pp2,psths_dt(:,cond,cell),pp',x,aff),pp20,LB,UB);
                   
                    picks_fit(:,cond,cell)=pp;cd '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data'
                    picks_fit(1,cond,cell)=pp1(1);
                  else
                    picks_fit(:,cond,cell)=ones(1,3)*nan;
                end
                
           end
%         %---------------------------------------------------------------------------  
           
           
         end
         
         
         for cond=1:12
         %****** DISPLAY *********
                 if disp   
                     figure(1)
                     
%                       pps=[5 1 2 3 4 8 12 16 15 14 13 9 6]; % order in the figure
                       pps=[5 9 13 14 15 16 12 8 4 3 2 1];
                      subplot(4,4,pps(ord(cond)))
                      plot(1:10:2500,prs(:,cond,cell));
                      hold on
                      plot(x, psths_dtS(:,cond,cell),'color',[0.8 0.7 0.8]);
                    
                      
                     % ----- DISPLAY RASTER      
                      for i=1:ct2(cond)
                          scatter(R_dt(1:100,i,cond,cell),ones(100,1)+i*5,32,'k.')
                          hold on
                      end
                      %------------
                      
                     title(num2str(cond))
                    if W(cond,cell)>thr
                     plot(pick(cond,cell),W(cond,cell),'g*')
                     
                    if fit_psths
                          m5=picks_fit(:,cond,cell);
                    
                       plot(x,m5(3)*(m5(2)*sqrt(2*pi))*normpdf(x,m5(1),m5(2)),'k')
                     end
%                      elseif W(cond,cell)<thrd
%                      plot(pick(cond,cell)*10,W(ord(cond),cell),'r*')
%                      line([pick(ord(cond))*10-100 pick(ord(cond))*10-100],yl); % this are the same param I used in DECODING_2011.m
%                      line([pick(ord(cond))*10+300 pick(ord(cond))*10+300] ,yl);    
                     end
                     
                    % plot(fit1,prs(fit1(cond),cond,cell))
                     xlim([0 1000])
                     ylim([0 max(W(:,cell))+5])
                     j40=line([0 2000] , [thr thr]);
                     %set(j40,'--')
                     
                     
                     
                     grid
               end
               
         clear pick2  
           end
           
           
           
          
   % ----------- FIT -------------------------------------------------------         
   clear limm1 good pik ord2
   good=find(W(:,cell)>thr);
    
    if numel(good)>3% if there are at least 3 points
       
        if fit_psths
            
            limm1=picks_fit(1,:,cell)';
        else
                    
        limm1=pick(:,cell);
        end
        
        
         limm=limm1(good);
         
         ff=ff1(good);
         
         pik=W(good,cell)';
         
         ord2=ord(good);
         
         options=optimset('MaxFunEvals',50000000);%,'display','iter'
         coss0=[1 1 1];
        % cosa=fminsearch(@(cosa) fitcosW(cosa,ff',limm,pik'),coss0,options);% weighted for the amount of response (pick)!!!
         cosa=fminsearch(@(cosa) fitcosW(cosa,ff',limm,pik'),coss0,options);

         clear fitt
         fitt=(cosa(1)+cosa(2)*cos(ff1+cosa(3)));
      %  fitt(find(fitt)<400)=400;
   %-------------------------------------------------------------------------

    
    % MEAN VECTOR
%         
         clear uus vvs mmx mmy th2 R2 uusM vvsM
         [uus,vvs]=pol2cart(ff1(ord)',fitt');
         [mmx mmy]=pol2cart(pi,(fitt(10)-fitt(4))/2);

         [th2,R2]=cart2pol(nanmean(uus),nanmean(vvs)); % CONVERT TO POLAR COORDINATE THE AVERAGE OVER THE 360 dimension vectors        
         [uusM,vvsM]=pol2cart(th2,R2);
      
        if disp  
            for cond=1:12
             subplot(4,4,pps(ord(cond)))  
                yl=ylim;
                line([fitt(cond)-200 fitt(cond)-200],yl); % this are the same param I used in DECODING_2011.m
                line([fitt(cond)+200 fitt(cond)+200] ,yl);      
            end
         
         
         yu=subplot(4,4,6)
         h=polarFC(ff1(ord2)',limm,'*g', max(limm),1);
         hold on
        
         h=polarFC(ff1(ord)',fitt','*k', max(fitt),1);
                   
         set(gca,'xdir','reverse')
        
         header  (['DIRECTION TUNING  horiz offset= ' num2str((fitt(10)-fitt(4))/2) ' ms  ((fitt(10)-fitt(4))/2)']);
         
         quiver(0,0,uusM,vvsM,3,'k','linewidth',3)
        
         quiver(0,0,mmx,mmy,3,'b','linewidth',3)
          
         grid
        
         quiver(zeros(12,1),zeros(12,1),uus,vvs,'r','ShowArrowHead','off')
        % legend('psths picks', 'cos fit', 'mean vector', 'Location' ,'NorthEastOutside') %'horiz component'
        
         ax=get(yu,'Position');
         ax([3,4])=ax([3,4])+0.15; %or wathever
         ax([2])=ax([2])-0.15;
          ax([1])=ax([1])-0.1;
         set(yu,'Position',ax);
         axis square
%          
%          ax=get(figure(1),'Position');
%          ax([3])=ax([3])+100; %or wathever
%          ax([4])=ax([4])+500;
%          set(figure(1),'Position',ax);
%          set(gcf, 'Units','centimeters', 'Position',[0 0 50 35])
%          set(gcf, 'PaperPositionMode','auto')
%          footer('from -200 to + 200 respect to the pick - fitcosW - fit_psths')
%          %saveas(gca,['/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/FIGURES/RECENTER/' num2str(cell) '.png'] )
        end   
         clear re reas
         
         
         %TUNING PLOT
         for cond=1:12
             i=ord(cond);  % here I calcule the integral of the response for the different conditions
             winP=fitt(cond)-200;
             winP2=winP+400;
             re(i)=numel(find((R_dt(:,:,cond,cell)>winP) & (R_dt(:,:,cond,cell)<winP2)))/ct2(cond)/400*1000
             reas(i)=numel(find(R_dt(:,:,13,cell)>0))/ct2(13) ;
             
         end
         
         re(13)=re(1); % to "close" the polar plot
         reas(13)=reas(1);
         
         re=re-reas; % BASELINE SUBTRACTION
         
   
     if disp     
         gg2=subplot(4,4,7)
        
         h=polarFC(pi/6:pi/6:2*pi+pi/6,re,'k', max(re),1);
         hold on
         polarFC(pi/6:pi/6:2*pi+pi/6,reas,'r', max(reas),1);
         
         set(h,'linewidth',3)
       %  header(['DIRECTION TUNING cell ' num2str(bb) ' ch' num2str(oo) ])
         hold on
        % quiver(0,0,cos(thetaD*pi/180), sin(thetaD*pi/180),max(re));
         set(gca,'xdir','reverse')
        
         ax3=get(gg2,'Position');
         ax3([3,4])=ax3([3,4])+0.15; %or wathever
       %  ax([1,2])=ax([1,2])-0.1;
          %ax3([2])=ax3([2])-0.15;
           ax([1])=ax([1])+0.25;
         set(gg2,'Position',ax);
         axis square
         
         saveas(figure(1),['/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/FIGURES/MEA/BM_MEA_20130927_dt1/RECENTER/' num2str(cell) '.png'] )
        end
   %offsetRFx(cell)=bx;
   offset(:,cell)=fitt;
   
   tuning(:,cell)=re;
   else
       
           %offsetRFx(cell)=0;
           offset(:,cell)=ones(12,1)*400;
    end

end


% save('tuning','tuning')
% save('offset','offset')

%    
 %  
 %   end
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------


 
 %% TEST ON SINGLE CELL with correction
 
 addpath '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/SCRIPTS'
 
 %clear R index
 %exp='BM_MEA_20130927_3T1';
 %pathl=(['W:\invibe\ANALYSIS\EXTRACELLULAR\USERS\GIACOMO\Mat_data\' exp '\']);
 pathl=(['C:\DATABASE\MATLAB\' exp '\']);
%  cd([pathl])
  load ([pathl 'R']);
  load ([pathl 'index']);
%  load index
%  clear R_3t
 R_3t=R;
 
 cond_n=max(index) ;

 ep_n = hist(index,1:cond_n);

 ct2=hist(index,cond_n)
 bin=10;
 
 for cell=1:96
     for cond=1:cond_n
         %ct2(cond)=max(find(R_3t(1,:,cond,cell)>0));
         psths_3t(:,cond,cell)=zeros(250,1);
         
         for ep=1:ct2(cond)
             psths_3t2(:,ep,cond,cell)=hist(R_3t(:,ep,cond,cell),1:bin:2500);
             psths_3t(:,cond,cell)=psths_3t(:,cond,cell)+psths_3t2(:,ep,cond,cell)*100;
         end
         
         psths_3t(1,cond,cell)=0;
         psths_3t(:,cond,cell)=psths_3t(:,cond,cell)/ct2(cond);
         
     end
 end
 
 
%%
 
 clf
 clear pd DEL
 %--------- MANUAL SETTING --------------------------
 
 exp='BM_MEA_20130927_dt1';
 
 disp=0; % DISPLAY IMAGES?

 correct=1

 col=['k' 'b']
 
 speed=13.3;
 
 dist=[1.5 1.5 1.5 6 6 6 6];
 
 
 
 ori=1:7;
 
 cond_n=numel(ori)*numel(dist);
 %-----------------------------------------------------
 
 
 pathl=(['W:\invibe\ANALYSIS\EXTRACELLULAR\USERS\GIACOMO\Mat_data\'  exp '\']);
 
 load([pathl 'offset']) % exp dt
 
 
 %**************************************************************************
  del=round((offset(1,:)-offset(6,:)./2./(offset(1,:)+offset(6,:))*900)./10);
 %**************************************************************************
 
 
  % del=round((offset(1,:)-offset(6,:)./2)./10);
 
 %trajj=[0 0 0 0 0 227 227 227 227 227 682 682 682 682 682];
 %uu=[2 5];
 
 
 
 
  cond_n=numel(dist);
   for cell=1:96
       clf
        ct=1;   
        for cond=1: cond_n
          
          if cond<cond_n;
          step=round((dist(end)-dist(cond))/speed*100)
          else
          step=1;
          end
     
    
       if correct
           del(cell)
         % del(cell)=40  
        %------------------------------------------------------  
          if del(cell)==0
             elseif del(cell)<0
             
               
                pd(del(cell)*-1+1:length(psths_3t),ori(cond),cell)=psths_3t(1:length(psths_3t)+del(cell),ori(cond),cell); % here I add time
                               
             else
                 
%                 if cond==dist(1)  % 1.5 traj
%                    pd(:,cond,cell)=nan;%psths_3t(:,cond,cell);
%                    pd(:,cond,cell)=psths_3t(:,cond,cell);
%                  else
                        
                 pd(1:length(psths_3t)-del(cell),ori(cond),cell)=psths_3t(del(cell)+1:length(psths_3t),ori(cond),cell);% here I remove time
               
                %end
           end
         %------------------------------------------------------   
           
           if disp
              hold on
              plot(step:length(pd(:,ori(cond),cell)),pd(1:length(pd(:,ori(cond),cell))-step+1,ori(cond),cell),col(cond))
     
           end
      
       
       else
            if disp
            hold on
             plot(step:250,psths_3t(1:250-step+1,ori(cond),cell),col(cond))
             
            end
       end       
            
         % plot(1:250,psths_3t(1:250,,40),'k')
           if disp  
              
              yl=ylim;
              line([40+del(cell)+step 40+del(cell)+step],yl)
              line([80+del(cell)+step 80+del(cell)+step],yl)
              
              l= line([40+step 40+step],yl);
               set(l,'Color',[1 0 0])
              l=line([80+step 80+step],yl);
               set(l,'Color',[1 0 0])
%            
               grid
          
             
             xlim([0 100])
             
             ct=ct+1;
             
              cd '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data'
            if exist ([exp])<1
                mkdir ([exp])

            else
                cd([exp]);
            end

           
         saveas(gcf,[num2str(cell) '_3T.jpg'])
           end  
         end
        
     end
 
  % header (['pick.m re-centering test - cell ' num2str(cell)])  
  % footer('in blue you have the resonse after recentering, in red before - this cell in dt had a RF shifted on the left')
%% **************************************************************************************************
%*****************************************************************************************************

% 3T MEAN
clf
col=['k' 'b' 'r' 'g' 'y'];

%pdM=mean(pd,3)

pdM=mean(psths_3t,3)

pdM=squeeze(psths_3t(:,:,40))

%%
disp=0;

col=['k' 'b' 'r' 'g' 'y'];

clf
for cell=1:96;

pdM=squeeze(psths_3t(:,:,cell));

oriS=[1 2; 3 4; 5 6]

for i=1:3
ori=oriS(i,:)
subplot(3,1,i)
speed=[ 6.6 6.6 sp(1) sp(1) sp(2) sp(2) ] 
for cond=1:numel(ori)
    
    pdM(:,ori(cond))=smoothfredo(pdM(:,ori(cond)),10);
%if ori(cond)==3
 if cond==2
    step=4.5/speed(ori(cond))*100
else
    step=1
 end
        if disp 
         h=plot(step:length(pdM(:,ori(cond)))+step-1,pdM(:,ori(cond)),col(cond))
         hold on
         set(h,'linewidth',3)
         title(['speed = ' num2str(speed(ori(cond)))] )
         grid
         xlim([30 150])
        end
        
  % ------ anticipation ----------------------------      
   st=round(step)+40
   bef(ori(cond),cell)=sum(pdM(st-10:st,ori(cond)));     
   plot(bef(ori(cond),:),[col(cond) '+']);
   hold on
  %------------------------------------------------- 
end


end
%header(['exp 20130923\_3T3\_Speed.DAT'])
%saveas(figure(1), ['C:\DATABASE\ANALYSIS\MEA\BM_MEA_20130923_3T3_speed\ch' num2str(cell) '.jpg'])
end
%%




%% display TUNINGS

%load('tuning.m')

for ch=1:96
re=tuning(:,ch);
if sum(re>0)
subplot(10,10,ch)
h=polarFC(pi/6:pi/6:2*pi+pi/6,re','k', max(re),0);
ax=get(gca,'Position');
ax(3)=ax(3)+0.015;
ax(4)=ax(4)+0.015;

set(h,'linewidth',2);
set(gca,'Position',ax);
end

end

%%
clf
for ep=1:15
cond=1
ch=40
subplot(5,3,ep)
plot(L(:,ep,cond,ch))
hold on
%subplot(2,1,2)
scatter(R(1:100,ep,cond,cell),ones(100,1),32,'r*')
 xlim([0 800])    
 grid
end


%%
clf
lfp=squeeze(mean(L,2));
cond=6
ch=40

plot(-lfp(:,cond,ch),'r')
hold on
plot(1:10:2500,psths_3t(:,cond,ch))
grid


xlim([0 1000])
%%
plot(-lfp(:,2,ch))
hold on
plot(-lfp(:,5,ch),'r')

