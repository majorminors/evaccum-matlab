function [Bw,Bl,enaw,enal,sim_act,t,t0] = expected_activityLBA(mu,sigma,theta,t0,c0,rt,n_units,finger,avoid,fsample,b01,shape)
%estimated neuronal activation
%NOTE: model rt ->sec ; 



%envelope sample length = 0.002 sec
if ~exist('fsample','var'); fsample = 0.002;end
if ~exist('b01','var'); b01 = t0/2;end
if isempty(b01);b01 = t0/2;end
if ~exist('shape','var');shape = 1;end

%rounding issues with fmincom...
rt =      round(rt*1000)/1000;
t0 =      round(t0*1000)/1000;
c0 =      round(c0.*1000)/1000;
fsample = round(fsample*1000)/1000;



%Inverse Mill's Ratio
IMR = @(beta,mu,sigma) mu-sigma.*(normpdf((beta-mu)/sigma)./normcdf((beta-mu)/sigma));


%time from 0 to end of ramp
t = 0:fsample:(rt-t0);

%expected values
Bw = (theta(finger)-c0(finger)/2)./(rt-t0);%accumulating rate for the winner unit
enaw = Bw.*t;
tmp_act = enaw;
enal = {};
Bl  = 0;

if n_units >1
    %let's take into account the distributions of each finger (IMR calculated
    %separately, enal estimated individually and then summed)
    units = 1:4; units([finger,avoid]) =[];
    
    for iu = units

        %expected value for bl for each loser unit
        Bl(iu) = IMR(Bw,mu(iu),sigma(iu));%#ok
        
        %expected accumulated activity from loser unit
        enal{iu} = Bl(iu).*t;%#ok
        
        %total expected activity
        tmp_act  = tmp_act + enal{iu};
    end
    
end


%%center the activity and pad t0 with zeros
t = 0:fsample:rt;

sim_act = zeros(length(t),1);
%onoff = dsearchn(t',[b01;rt-(t0-b01)]);
onoff = dsearchn(t',b01);
onoff(2) = onoff(1)+length(tmp_act)-1;

sim_act(onoff(1):onoff(2)) = tmp_act;
if shape ==1
sim_act(onoff(2)+1:end) = tmp_act(1);
elseif shape == 2
sim_act(onoff(2)+1:end) = tmp_act(end);    
elseif shape == 3
sim_act = tmp_act;    
end
