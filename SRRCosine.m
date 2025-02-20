function [filter_coeff,TotalTap,t]=SRRCosine(NISI,Alpha,Tc,Beta)
%%
% A routine for Square Root Raised Cosine pulse shape.
% The Tap value are normalized in energy.
% - "NISI" is the number of pre-cursor and post-cursor
% - "Alfa" is teh roll-off factor
% - "Tc" is the sampling time normalized to 1
% - "Beta" is the oversampling factor
% - "filter_coeff" is the array containing the filter coefficients
% - "Totaltap" is the number of filter coefficients, i.e.,
% length(filter_coeff)
% - "t" is the time vector at the symbol rate
% Example: [filter_coeff,TotalTap,t]=SRRCosine(3,0.3,1,4);



%% Time axis computation
% Creation of an auxiliary array for the routine to compute SRRC coefficients
% "start" = Initial Point to calculate TapValue, depend on the interferent symbols
TotalTap = (2*NISI*Beta)+1;	% Number of Tap for SRRC Filter
t=[-NISI:1/Beta:NISI]; %sampling points
t=t-0.000001; % to avoid computation of  undetermined ratios 

%% Generation Tap Value
tnorm=t./Tc;
Num=sin((1-Alpha)*pi.*tnorm)+4*Alpha*cos((1+Alpha)*pi.*tnorm).*tnorm;
Den=pi.*tnorm.*(1-16*Alpha^2.*tnorm.^2);

filter_coeff_not_norm=Num./Den;
%Energy computation
Energy_SRRC=sum(filter_coeff_not_norm.^2);
% Energy Normalizzation
filter_coeff=filter_coeff_not_norm/sqrt(Energy_SRRC);



