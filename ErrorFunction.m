%INPUTS FROM YOU
vGoal1 =[]; % The voltage time series of goal-neuron 1 that you want to emulate
vGoal2 = []; % The voltage time series of goal-neuron 2 that you want to emulate
vModel1 = []; % The voltage time series of model-neuron 1 that is trying to mimic goal-neuron 1
vModel2 = []; % The voltage time series of model-neuron 2 that is trying to mimic goal-neuron 2
sModel1 = []; % The s time series of model-neuron 1 which should be one of the variables of the model-neuron
sModel2 = []; % The s time series of model-neuron 2 which should be one of the variables of the model-neuron

%PREALLOCATING THE BLENDED SYNAPTIC MODEL OF GOAL NEURONS
sGoal1_Array = zeros(length(time1),1); 
sGoal2_Array = zeros(length(time1),1); 
%FIND THE PEAKS IN THE GOAL-NEURON VOLTAGE TIME SERIES
peaks_vGoal1 = findpeaks(vGoal1,'MinPeakHeight',10, 'MinPeakDistance',5);
peaks_vGoal2 = findpeaks(vGoal2,'MinPeakHeight',10, 'MinPeakDistance',5);

%CALCULATE THE BLENDED SYNAPTIC MODEL OF GOAL NEURONS
for i=1:length(time1) 
    % CALCULATE sGoals USING THE BLENDED SYNAPTIC MODEL EQUATION
    sGoal1=sGoal1+step*(alpha*(1-sGoal1)/(1+exp(-10*(vGoal1(i)+20)))-beta*sGoal1);
    sGoal2=sGoal2+step*(alpha*(1-sGoal2)/(1+exp(-10*(vGoal2(i)+20)))-beta*sGoal2);
    % STORE sGoals IN THEIR RESPECTIVE ARRAYS
    sGoal1_Array(i)=sGoal1;
    sGoal2_Array(i)=sGoal2;
end  

%STORE POST-PROCESSED VARIABLES
post_sGoal1_Array=sGoal1_Array;
post_sGoal2_Array=sGoal2_Array;
smooth_post_vGoal1= movmean(vGoal1,3000);
smooth_post_vGoal2= movmean(vGoal2,3000);

post_vModel1=squeeze(vModel1); 
post_vModel2=squeeze(vModel2); 
sModel1=squeeze(sModel1); 
sModel2=squeeze(sModel2); 
smooth_vModel1= movmean(post_vModel1,3000); 
smooth_vModel2= movmean(post_vModel2,3000);
smooth_sModel1= movmean(sModel1,3000); 
smooth_sModel2= movmean(sModel2,3000); 

        
%CALCULATE ENVELOPES
if isfinite(post_vModel1)
    [voltsvGoal1upper,voltsvGoal1lower] = envelope(vGoal1,150,'peak');
    [voltvModel1upper,voltvModel1lower] = envelope(post_vModel1,150,'peak');
    envNormUpper1 = norm(voltvModel1upper-voltsvGoal1upper);
    envNormLower1 = norm(voltvModel1lower-voltsvGoal1lower);
    envNormTotal1 = (envNormUpper1+envNormLower1)/2;

    [voltsvGoal2upper,voltsvGoal2lower] = envelope(vGoal2,150,'peak');
    [voltvModel2upper,voltvModel2lower] = envelope(post_vModel2,150,'peak');
    envNormUpper2 = norm(voltvModel2upper-voltsvGoal2upper);
    envNormLower2 = norm(voltvModel2lower-voltsvGoal2lower);
    envNormTotal2 = (envNormUpper2+envNormLower2)/2;


else
    [voltsvGoal1upper,voltsvGoal1lower] = envelope(vGoal1,150,'peak');
    voltvModel1upper = voltsvGoal1upper;
    voltvModel1lower = voltsvGoal1lower;
    envNormUpper1 = norm(voltvModel1upper-voltsvGoal1upper);
    envNormLower1 = norm(voltvModel1lower-voltsvGoal1lower);
    envNormTotal1 = (envNormUpper1+envNormLower1)/2;

    [voltsvGoal2upper,voltsvGoal2lower] = envelope(vGoal2,150,'peak');
    voltvModel2upper = -1*voltsvGoal2upper;
    voltvModel2lower = -1*voltsvGoal2lower;
    envNormUpper2 = norm(voltvModel2upper-voltsvGoal2upper);
    envNormLower2 = norm(voltvModel2lower-voltsvGoal2lower);
    envNormTotal2 = (envNormUpper2+envNormLower2)/2;               
end
        

%CALCULATE ERRORS
error1 = (norm(post_sGoal1_Array - post_sModel1) + norm(post_sGoal2_Array - post_sModel2))/2;
peaks_post_vModel1 = findpeaks(post_vModel1,'MinPeakHeight',0, 'MinPeakDistance',5);
peaks_post_vModel2 = findpeaks(post_vModel2,'MinPeakHeight',0, 'MinPeakDistance',5);
error2 = abs(abs(length(peaks_vGoal2)-length(peaks_post_vModel2))+abs(length(peaks_vGoal1)-length(peaks_post_vModel1)))./2;
error3 = norm(((smooth_post_vGoal2-smooth_vModel2)+(smooth_post_vGoal1-smooth_vModel1))./2);
error4 = var(((abs(smooth_post_vGoal2-smooth_vModel2)+abs(smooth_post_vGoal1-smooth_vModel1))./2));
error5 = sigmf(((envNormTotal1 + envNormTotal2)/2),[0.001 9.2039e+03]);

%RESCALE COSTS
error1_max = 104.7806 + 0.0001;
error1_min = 0 - 0.0001;
error1_rescaled = rescale(error1, 'InputMin', error1_min, 'InputMax', error1_max);

error2_max = ((length(peaks_si3L)+length(peaks_si3R))/2) + 0.0001;
error2_min = 0 - 0.0001;
error2_rescaled = rescale(error2, 'InputMin', error2_min, 'InputMax', error2_max);

error3_max = 9.2767e+03 + 0.0001;
error3_min = 0 - 0.0001;
error3_rescaled = rescale(error3, 'InputMin', error3_min, 'InputMax', error3_max);

error4_max = 75.9981 + 0.0001;
error4_min = 0 - 0.0001;
error4_rescaled = rescale(error4, 'InputMin', error4_min, 'InputMax', error4_max);

error5_max = 1 + 0.0001;
error5_min = 0 - 0.0001;
error5_rescaled = rescale(error5, 'InputMin', error5_min, 'InputMax', error5_max);

%COMBINE RESCALED ERRORS
w1 = 30; % synapse norm
w2 = 5; % spikes
w3 = 5; % smooth voltage norm
w4 = 10; % variance of smooth voltage difference
w5 = 2; % envelope norm

weighted_error_combo = (w1*error1_rescaled + w2*error2_rescaled + w3*error3_rescaled + w4*error4_rescaled)/(w1 + w2 + w3 + w4);
weighted_error_combo_rescaled = (w1*error1_rescaled + w2*error2_rescaled + w3*error3_rescaled + w4*error4_rescaled + w5*error5_rescaled)/(w1 + w2 + w3 + w4 + w5); 

%SET ERROR Set cost as the rescaled weighted combination
error = weighted_error_combo_rescaled; 
error_max = 1;