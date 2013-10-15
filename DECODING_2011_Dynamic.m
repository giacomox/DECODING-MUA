%% DECODING_2011 INHERITED BY DIRECTION TUNING EXPERIMENTS 2011 - population decoding analysis
% by Giacomo Benvenuti - 14/9/2013  - Invibe - INT - CNRS

% this script is herited from "dir_tun_2011_PLOODEC.m "






%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% TUNING FIT

clear all

%------------ SET MANUAL PARAMS !!!!!!!!! ------------
display_image=0; %DO YOU WANT TO DISPLAY THE FUGURE?

spkN=0; % do u want to use the number of spks or the spk rate?

N=70; % how much do you want to smooth the psth?

stim_duration=2000;

lag=200 ; % size of the response window considered

bin=10;


%------------------------------------------------------

addpath '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/SCRIPTS'
giacstyle; % function with fontsizes and some other trivial stuff
paths= '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data/order/'

%----- LOAD global---------------------------------------------------------
    load([paths 'fnl']); % LOAD CELLS LIST
    load([paths '20111005_ch5_R_3t.mat']) % LOAD POPULATION RASTERS
    load([paths '20111005 ch5_R_dt.mat'])
    load([ '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data/22-sep-2013/offset']) % offset of the response when the RF is not well centered
%--------------------------------------------------------------------------    

nb_cond=size(R_dt,3);% 12 ori + blanc

theta = (pi/6:pi/6:2*pi); % ori of the stimuli in radiants

num_c= size(fnl,1); % cells number

hzz=0:1:200;

for cell=1:num_c
    
    
    %-------- LOAD 3t  --------
    bb=fnl(cell,69:76); % recover exp date
    oo=fnl(cell,87:87);% here I recover the channel numb
    
    load([paths num2str(bb) '_ch' num2str(oo) '_dir_dt.mat']) % index of condition for the exp dt
    load([paths num2str(bb) '_ch' num2str(oo) '_cond_3t.mat']) % " " for 3t
   
    %-------------------------------
    
    RASTERS=permute(R_3t,[2 1 3 4]); % xxxx INPUT xxxxx
    RASTERSD=permute(R_dt,[2 1 3 4]);
    
    nbrep(:,cell)=hist(cond_3t,1:16);% assays number in 3t raster
    
    nbrepDt(:,cell)=hist(dir_dt,1:13);
   

          

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% CALCULE SPIKE RATE DIRECTION TUNING  ------------------------------------------------------

clear sspkAS

for th=1:12
    for as=1:nbrepDt(th,cell)
        sspkAS(as,th)=numel(find((RASTERSD(as,:,th,cell)<400) & (RASTERSD(as,:,th,cell)>1))); % baseline before stimulus appear - Spontaneous activity
    end
end

if spkN
    spk_rateAS(cell)=mean(sspkAS(:,th))
else
    
    spk_rateAS(cell)=mean(sspkAS(:,th))./400*1000;
end

for th=1:12
    
   
    %------ RESPNSE WINDOW SELECTION --------------------------------------
    if (offset(th,cell)-lag/2)<400;
         win1(th)=400;
    else
         win1(th)=offset(th,cell)-lag/2; % offset come from the fit of DT smooth psth max
    end
    
    %win1(th)=400; % 26/9
     win2(th)=win1(th)+lag;  % lag is set by the user
    %----------------------------------------------------------------------
    
  
    
    num_spk(th,cell)=numel(find((RASTERSD(:,:,th,cell)>win1(th)) & (RASTERSD(:,:,th,cell)<win2(th))));
   
    if spkN
        spk_rate(th,cell)= num_spk(th,cell)/nbrep(th);
    else
        spk_rate(th,cell)= num_spk(th,cell)/nbrep(th)/(win2(th)-win1(th))*1000; % skp/sec
    end
    
    for ep=1:nbrepDt(th,cell)
        qq=numel(find((RASTERSD(ep,:,th,cell)>win1(th)) & (RASTERSD(ep,:,th,cell)<win2(th))));
        if spkN
            spk_epD(ep)= qq;
        else
            spk_epD(ep)= qq./lag*1000 ;
        end
        %   qq(ep)= spk_ep(ep,th,cell); % <<<<<<<< Poisson likes integers - after meeting with Lolo 10/9/2013
        
       
    end
    
    
    P(:,th,cell)=hist(spk_epD,hzz); % the input fort the fit is the distrib of rate in a raster
    
    clear spk_ep_rate qq spk_epD
    
    
end

tuning=spk_rate(:,cell)' ;%- spk_rateAS(cell)' ;


%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
% FITTING -multiplicative von mises fitting -

[maxx, i] = max(tuning);
if maxx>100
    maxx=50;
end

% 1) 
%--------------------------------------------------------------------------
% a0=[theta(i),spk_rateAS(cell),maxx,2,1,1];
% options=optimset('MaxFunEvals',500000000) ;
% LB=[theta(i)-pi/2,spk_rateAS(cell)/2,0,0,0,1];%lower bound vector (same order than a0)for iteration
% UB=[theta(i)+pi/2, spk_rateAS(cell)*2,100,10,10,5];
% aff=0;
% 
% a1=fminsearchbnd(@(a1) SEKLG(a1,theta, P(:,:,cell),hzz,aff),a0,LB,UB,options) ;  % try nbinpdf negative binomial Goris - Movshon
%--------------------------------------------------------------------------



% 2) Theta Dir free to change - 6 params
%--------------------------------------------------------------------------
a0=[theta(i),spk_rateAS(cell),maxx,2,1,1,theta(i)];
options=optimset('MaxFunEvals',500000000) ;
LB=[theta(i)-pi/2,spk_rateAS(cell)/2,0,0,0,1,-1000];%lower bound vector (same order than a0)for iteration
UB=[theta(i)+pi/2, spk_rateAS(cell)*2,100,10,10,5,+1000];
aff=0;

a1=fminsearchbnd(@(a1) SEKLG6(a1,theta, P(:,:,cell),hzz,aff),a0,LB,UB,options) ; 
%--------------------------------------------------------------------------



fits_params(:,cell)=a1;% OUTPUT - array of param of the fitted function
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

end









%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% MAXIMUM LIKELYHOOD EXTIMATIONS


step=30;

eee=jet;
clf
thetafit= (0:pi/60:2*pi)';
for cell=1:80
        [fitx(:,cell) dirrr(:,cell) ori(:,cell) Rm]=VMm(fits_params(:,cell),thetafit); %thetafix
end

clear re ru reB ruB MLE MLED

for tm=1:30 % time steps
  
  for ep=1:10
    clear spk3t spk3tB spkDt
    
     for cell=1:80
        clear win1 win2
        
        % CALCULE SPIKE RATE 3t ------------------------------------------------------
        
        tth=[1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12];
        
        trajj=[0 0 0 0 0 227 227 227 227 227 682 682 682 682 682];
        
        
        for th=1:15
%           %*****************  3t    ****************************
          
            %-----------------------------------------------------------------------
             del=(offset(10,cell)-offset(4,cell))./2);% re-center responses based on the exp dt
           
            %-----------------------------------------------------------------------
           
            if (400-del+trajj(th) < 400;
                win1(th)=400+trajj(th)   ;%????
            else
                win1(th)=400+del+ trajj(th); 
            end
            
            %Dynamic
            win1=win1-400+(step*tm);
            
            win2(th)=win1(th)+lag   ;
            
          
            %-------------- DT ---------------------
            if th<13
            if (offset(th,cell)-lag/2)<400;
                winD1(th)=400;
            else
                winD1(th)=offset(th,cell)-lag/2; % offset come from the fit of DT smooth psth max
            end 
            
            %Dynamic
            winD1=winD1-400+(step*tm);
            winD2(th)=winD1(th)+lag   ;
            end
          %*********************************************   
           
            
         %------------  random assay selection -------------  
            sel=round(rand(1)*nbrep(th,cell));
            if th<13
            selD=round(rand(1)*nbrepDt(th,cell));% we want to select random mix of assaies from the population for each repetition
            if selD==0
                selD=1;
            end
            end
            if sel==0
                sel=1;
            end
        %--------------------------------------------------- 
            
            clear spk_ep spk_epB
            spk_ep=numel(find((RASTERS(sel,:,th,cell)>win1(th)) & (RASTERS(sel,:,th,cell)<win2(th))));
            
            spk_epB=numel(find((RASTERS(sel,:,th,cell)>win1(th)-lag) & (RASTERS(sel,:,th,cell)<win2(th)-lag)));
            
            if th<13
            spk_epD=numel(find((RASTERSD(selD,:,th,cell)>winD1(th)) & (RASTERSD(selD,:,th,cell)<winD2(th))));
            end
            
            % end
            
            if spkN
                spk3t(th,cell)=spk_ep; % single assay response  to paramtrize cell tuning
                %spk3tB(th,cell)=spk_epB;
               if th<13 
                spkDt(th,cell)=spk_epD;
               end
             else
                spk3t(th,cell)=spk_ep./lag*1000;
               % spk3tB(th,cell)=spk_epB./lag*1000;
                 if th<13 
                spkDt(th,cell)=spk_epD./lag*1000;
               end
            end
        end
    end
    
    
  
    %---------------------------------------------------------------------
    %  DECODING ----------------------------------------------------------
    %---------------------------------------------------------------------
    
    clear ML dirrr MLS
    
   
    
    disp=0;
    
    for cell=1:80
                
        for th=1:15                
                                    
                                      %JAZAYERI
         %*****************************************************************   
            ML(:,th,cell) = log(fitx(:,cell))*(spk3t(th,cell)); %    MAXIMUM LIKELYHOOD EXTIMATION SINGLE TRIAL
         %*****************************************************************      
            % the Von MIses is just the exp of a cosine function for this a
            % I normilize the tuning curves using the log
            
         
          if th<13  
            MLD(:,th,cell) = log(fitx(:,cell))*(spkDt(th,cell));
          end
            
        end
        
    end
    
    MLS=squeeze(sum(ML,3));
   
    MLSD=squeeze(sum(MLD,3));
    
    kern=sum(fitx,2); % KERNEL 
    
   
    for i=1:15
        
        
                                       %JAZAYERI
        %*****************************************************************   
        MLE(:,i,ep,tm)=MLS(:,i) - kern - sum(log(gamma(spk3t(i,cell)+1)),2); % REGULARIZATION 
       %*****************************************************************  
        
        [re(i,ep) ru(i,ep)]=max(MLE(:,i,ep));
        [min_y(i,ep) min_x(i,ep)]=max(MLE(:,i,ep));
        
        
        if i<13
        MLED(:,i,ep,tm)=MLSD(:,i) - kern - sum(log(gamma(spkDt(i,cell)+1)),2);  
        [reD(i,ep) ruD(i,ep)]=max(MLED(:,i,ep));
        
        end
       
        
        if disp
            subplot(3,5,i);
            hold on
            plot(0:3:360,MLE(:,i,ep),'Color',eee(ep*4,:));
            
            
            %  ylim([500 3000])
            xlim([1 360]);
            mylin=line([90 90],ylim);
            xlabel('ori')
            ylabel('sum(response_cell_orix*f(tuning)_i)')
            grid
            % ylim([-3 6])
        end
    end
    
    %    end
    %header('LEAVE ONE OUT EXTIMATION - 10 REPET')
    %footer('MLE(:,i,rep)=MLS(:,i) - kern - sum(log(gamma(looSR(i,:,rep)+1)),2);')
    
end
end




%% DISPLAY DECODING : PICKS OF THE MEAN MLE


close all
figure(1) 

cd '/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/FIGURES/movie/'
filename = 'DT_DECOD_dynamic.gif';
set(gcf,'Position',[500 500 2000 1000])

MLEM=(mean(MLE,3));
MLEDM=(mean(MLED,3));

C=[63:3:360, 3:3:60]

n=1;
clear ss aa mov frame
for tm=1:30
    
    
    clear ss aa
    
    
    GG=MLEDM(:,:,tm); % set here if u want to show DT or 3t
   
    
    [aa ss]=max(GG);
    %ss(find(ss<10) )=120;                             %
    
    
    clf
    
    aa(1)=2000
    
    if dt
        nn=1;
    else
        nn=1:3
    
    for tr=nn
      
      if dt==0  
        subplot(1,3,tr)
      end   
        g1=1+5*(tr-1)
        g2=g1+4;
        col=['b' 'c' 'r' 'g' 'y' 'b' 'c' 'r' 'g' 'y' 'b' 'c' 'r' 'g' 'y' ]
        ct=1;
        
        if dt 
            nn2=1:12
        else
            nn2=g1:g2
        end
        
        for i=nn2
            % polar(ss(g1:g2).*3,aa(g1:g2),'*')
            h=polar([0,C(ss(i))*pi/180],[0,aa(i)],[col(ct)])
            set(h,'linewidth',3)
            hold on
            h2=polar([0,(C(ss(i))+180)*pi/180],[0,aa(i)],[col(ct)])
            ct=ct+1;
            %  polar(mod(ss(g1:g2)+180,360),aa(g1:g2),'r*')
        end
        hold off
       if dt==0
        title(['TRAJ = ' num2str(tr)])
       end
       
    end
    header(['TIME: ' num2str(30*tm) ':' num2str(30*tm+lag)])
    footer('10 REP - WIN 200 ms  - MLE MEANS PICKS ')
    
 figure(1) 
 

drawnow
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
  
 if n == 1;
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
else
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
end
 n=n+1;
end
    % close(mov);
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX














%% ************
figure(2)
clf
MLEM=(mean(MLE,3));
MLEBM=(mean(MLEB,3));

clear aa ss 
ct=1;
for i=1:5
  
    subplot(1,5,ct)
    plot(MLEM(:,i),'b')
    yl=ylim;
    [aa ss]=max(MLEM(60:120,i))
    line([ss+60 ss+60] ,[yl(1) aa])
    xlim([60 120])
    set(gca,'xtick',[0:10:120])
    grid
    ct=ct+1;
    hold on
end



%set(gca,'xtick',[0:10:120])
%legend('-60', '-30', '0', '30','60')
%footer('in number of spks - NO OFFSET CORRECTION - cond 11:15 MLEB BEFORE the RF - MEAN 100 REPETITIONS')
%grid
%saveas(figure(2),'/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/FIGURES/MLE11.jpg')




%%   DISPLAY DECODING 

% set(gca,'xtick',[0:30:360])

% HELP
% 'ru( 3tcond, repet)' are the max of the MLE inside the RF
% 're( 3tcond, repet)' are the hights of the max
% 'ruB( 3tcond, repet)' is the same than ru before the RF


% NOTE:in 'r' you find values from 0 to 120. it depend onf the fact we have
% 12 conditions and to generate the vonmises we used 120 point to have a
% better resolution. ru=10 is equals to 10*3, 30 deg



%********************* LOAD FOR TESTS *****************************************
cd /riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/Mat_data/26-sep-2013/
load re 
load ru
load reB
load ruB
%******************************************************************************

%%
clear vx vy 
clf
col=['r' 'b' 'k' 'g' 'c' 'r' 'b' 'k' 'g' 'c' 'r' 'b' 'k' 'g' 'c']; 

for tt=1:2                                  % before and after the RF
  clear rux rex  
    if tt==1
        rux=ru*3;
        rex=re;
        
       
    else
        rux=ruB*3;
        rex=reB;
    end
    
 %   ruxC=complex(rux,rex);
   % ruxC=mod(ruxC+60,360);
    
    
    rux=mod(rux+60,360);
    
    figure(3)
       
     ct=1;     
        for cond=1:15  
            hold on
      %  pospos=cond+(trajj-1)*5+(tt-1)*15;
        subplot(6,5,ct+15*(tt-1))
        clear hh fg
        hh1=hist(rux(cond,:)',0:3:360); 
        hh(1,:)=smoothfredo(hh1,10);
        hh(1,end)=0;
        hh(1,1)=0;
         
%          [rrx rry]=pol2cart(-(2*pi/3)+pi/6*cond,max(hh))
%         quiver(0,0,rrx,rry,0.5)
         
        fg=polarFC(0:pi/60:2*pi,hh,'r',max(hh),'off');
        set(fg,'Color',col(cond),'linewidth',3)
     
      % polarFC(-(2*pi/3)-pi/6cond,max(hh),'*b',max(hh),'off');
    % [ba1 ba2 ] =get(gca,'Position')
       ct=ct+1;
      %  end
    end
end
%%
for i=1:5
   subplot(6,5,i)
   title(['ORI ' num2str(90-i*30)])
    
end

for i=1:6
   subplot(6,5,1+5*(i-1))
   if i<4
   ylabel(['TRAJ ' num2str(i)])
   
   else
   ylabel(['TRAJ ' num2str(i-3)])
   end
    
end

%%












%%  OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD 








ORI=0;  % DIRECTIONS OR ORIENTATIONS DECODING?
clear vx vy 
clf
col=['r' 'b' 'k' 'g' 'c' 'r' 'b' 'k' 'g' 'c' 'r' 'b' 'k' 'g' 'c']; 

for tt=1:2                                  % before and after the RF
  clear rux rex  
    if tt==1
        rux=ru*3;
        rex=re;
    else
        rux=ruB*3;
        rex=re;
    end
    
    
    rux=mod(rux+60,360);
    if ORI
 
        aa=(find(rux>90 & rux<270));
       
        rux(aa)=mod(rux(aa)+180,360);
        rux2= rux(aa);
      %  [kk, bb]=(find(rux>180 & rux<210));
       
       % rux(bb)=rux(bb)-180;
        clear aa bb
    
        
   
       
    end
    
    
     % condition 10 (100) in dt is ori 90deg dir 180 deg --> is the same that cond 3 in 3t
    
    
    rum=mode(rux,2); % MEAN OR MODE FOR THE ORIENTATIONS VECTORS
    
        
    
    for ww=1:3                                % 3 trajectories
        
        figure(1)
        subplot(2,3,ww+(tt-1)*3)
        set(gca,'xdir','reverse')
         ct=1;
        for rep=1:100                         % repetitions
            ui=1+(ww-1)*5:5+(ww-1)*5;
            
            for cond=ui                          % orientations
               if rux(cond,rep)<90 | rux(cond,rep)>270 
                polarFC(rux(cond,rep)*pi/180,rex(cond,rep),[col(cond) '.'],max(re(:)),1);
                hold on
                  
                [vx(cond,ct) vy(cond,ct)]=pol2cart(rux(cond,rep)*pi/180,rex(cond,rep));
                   
                
                  %  quiver(0,0,vx(cond,rep),vy(cond,rep))
                ct=ct+1;
               end              
            end
            
        end
  
   
        
        
        if ww+tt==2
            legend('-60', '-30' ,'0' ,'30', '60','Location','NorthEastOutside')
        end
        
        text(max(re(:))+500,0,num2str(mod(rum(ui).*3,360)))
        
        for cond=ui                               % DISPLAY MODE VECTORS
           % [uus,vvs]=pol2cart((rum(i)*pi/180),mean(rex(i,:),2)); % directions
             hold on
            %[MMx(cond) MMy(cond)]=cart2pol(mean(vx(cond,:),2), mean(vy(cond,:),2)); 
            vx(find(vx<1))=nan;
            vy(find(vy<1))=nan;
            h=quiver(0,0,nanmedian(vx(cond,:),2) ,nanmedian(vy(cond,:),2));
             set(h,'Color',col(cond),'linewidth',3)
        end
        
       % set(gca,'fontsize',15)
    end
end

for yu=1:3
    subplot(2,3,yu)
    title(['TRAJECTORY ' num2str(yu)])
end


subplot(2,3,1); ylabel('INSIDE THE RF')
subplot(2,3,4); ylabel('BEFORE THE RF')

ax=get(figure(1),'Position');
ax([3])=ax([3])+300; %or wathever
ax([4])=ax([4])+100;
set(figure(1),'Position',ax);
set(gcf, 'Units','centimeters', 'Position',[0 0 50 35])
set(gcf, 'PaperPositionMode','auto')

if ORI
    header('ORIENTATIONS  DECODING')
else
    header('DIRECTIONS DECODING')
end
footer('100 random repetiotions - time win of 400 ms - script:GIACOMO/backup/DECODING_2011_4.m')

%%
clf
clear like
MLEM=mean(MLEB,3);
for i=1:12
    
    subplot(3,5,i)
    like(:,i)=(MLEM(:,i)-mean(MLEM(:,i))).*3;
    polarFC(0:pi/60:2*pi,like(:,i)','r',max(like(:,i)),'off');

end



%%  Angle
clf
dt=0;  % SELECT HERE DT OR 3T

if dt
    oris=1:12;
else
    oris=1:15;
end


for ori=oris
    
    subplot(3,5,ori)
    clear pp pp_x pp_y pp2 d1 xx l1 l2
    for tt=1:30 %time
        clear pp1
        
        if dt
            [pp1 pp2]=max(MLED(:,ori,:,tt));
            l1=ori*10-10;
            l2=ori*10+10;
            
        else
            ori
            [pp1 pp2]=max(MLE(:,ori,:,tt));
            l1=100-10;
            l2=100+10;
       
        end
        
        pp_y(tt,:)=squeeze(pp1);
        pp_x(tt,:)=squeeze(pp2);
        
       
        
        hits(tt)= numel(find(pp_x(tt,:)>l1 & pp_x(tt,:)<l2));
        fa(tt)=10-hits(tt);
        
        sm=mean(hits(tt),2);
        nm=mean(fa(tt),2);
        
        ss=std(permute(hits,[2 1]));
        ns=std(permute(fa,[2 1]));
        
        d1(tt)= (sm-nm)/sqrt(1/2*(ss^2+ns^2));
        
    end
    if dt
        xx =130:30:30*30+100;
    else
        xx=30:30:30*30;
    end
    %plot(xx,(pp_x-10)*3,'*')  %plot MLED
    
    %---------- FIG 1 -------------------
    
    %    figure(1)
    %     plot(xx,(pp_x-10)*3,'*') %plot MLE
    %    figure(2)
    %     plot(xx,pp_y,'*') % plot MLE absolute hight
    %     clear sm nm ss ns
    %     ylim([-50 360])
    %     xlabel('TIME (ms)')
    %     ylabel('directions (alpha)')
    %     header(['condition ' num2str(ori) ' 300'])
    %     title([num2str(ori)])
    %-------------------------------------
    
   %-------- FIG 3 -------------------------- 
    figure(3)
    plot(xx,d1)
    
    ylim([-10 5])
    yl=ylim;
    line([300 300],yl)
    grid
end





