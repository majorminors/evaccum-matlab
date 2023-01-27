function at_apply_ICA(inputf,outputf,ICAfile)

O = spm_eeg_load(inputf);
load(ICAfile);

%create a copy of the object with empty channels
new = clone(O,outputf,size(O),0);


modalities = {'MEGMAG', 'MEGPLANAR', 'EEG'};
newdata = zeros(size(new));

badchans = O.badchannels;

for m = 1:length(modalities)
    
    indCh = O.indchantype(modalities{m});
    indCh = indCh(~ismember(indCh,badchans));%if bad channels
    
    data = O(indCh,:,:);orig_size = size(data);
    data = reshape(data,orig_size(1),[],1);%checked OK
    
    W = all_ica.weights{m};remove = all_ica.remove{m};
    iW = pinv(W);
    
    IC = W*data;%unmix 
    
    if ~isempty(remove)
        IC(remove,:)=0;%remove bad ICs
    end
    
    data = iW * IC;%back-project
    
    newdata(indCh,:,:) = reshape(data,orig_size);%restore original dimensions
    
end
if ~isempty(badchans); newdata(badchans,:,:) = O(badchans,:,:);end %put back bad channels
new(:,:,:)=newdata;%populate channels
new.save;

%%check
% ch = 174;clf;
% subplot(2,2,1);plot(O(ch,1:1000));
% subplot(2,2,2);plot(O(ch,1:1000));hold on;plot(newdata((ch),1:1000),'r');
