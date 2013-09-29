%**************************************************************************
%                   IMPORT MEA DATA RECORED BY ELPHY
%
%        by Giacomo Benvenuti - Invibe - INT - CNRS Marseille
%        29-Sep-2013     <giacomox@gmail.com>
%
% HELP 
% This script allow you to inport in matlab DATA exported from the ELPHY
% format .DAT by the ELPHY script EXPORT_MEA.pg2
%
% It works on a single expteriment 
% -------------------------------------------------------------------------
%   INPUTS:
% * MUA_[channel nuber]_[episode] = Multi unit activity, timestamps
%                                   representig spks
% * LFP_[channel nuber]_[episode]
%
% * FIX_[channel nuber]_[episode] = eye moviments. 98 and 99 from the
%                                   channels in which the vertical and horizontal 
%                                    eye posistion was recorded
% * stim_index= vector to correlate any episode number with a stimulus number 
%               (which kind of stimulus was presented)
%--------------------------------------------------------------------------
% * OUTPUTS:
%
% * R(episode,time,condition,channel) = ME raster for spks
%
% * L(episode,time,condition,channel) = ME raster for LFPs
%
%--------------------------------------------------------------------------


clear all

%XXX set by the USER!! XXXXXXXXXXXXXX

exp='BM_MEA_20130927_dt';

cond_n= 13 ; % conditions number

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

pathx= (['/riou/work/invibe/ANALYSIS/EXTRACELLULAR/USERS/GIACOMO/fromELPHYMEA/' exp '/']);


% Inport the stimuli index and recover the repetiton number for every condition

index=LoadElphyVector([ pathx  'stim_index'],1,'Vector');

ep_n = hist(index,1:cond_n);  

% Inport MUA

R=ones(500,max(ep_n),cond_n,96)*nan;

for ep=1:ep_n
   
    for ch=1:96
      
        cond=index(ep);
      
        R(:,rep(cond),cond,ch)=LoadElphyVector([pathx 'MUA_' num2str(ch) '_' num2str(ep) ],1,'Vector');
    
        L(:,rep(cond),cond,ch)=LoadElphyVector([pathx 'LFP_' num2str(ch) '_' num2str(ep) ],1,'Vector');
    
    end
end








