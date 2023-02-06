'''
https://github.com/jungky10/RBD_actigraphy

KyoungeunPark
miniholic98@gmail.com

Ki-Young Jung
jungky@snu.ac.kr
'''

import pyActigraphy
import os
import csv
import pandas as pd
import numpy as np
import sys
import pickle
from matplotlib import pyplot as plt
import copy
import math
from scipy.signal import find_peaks
from FLM_func import FLM2
import statsmodels.api as sm

fpath = os.path.join(os.getcwd(), 'data_example/')
filename='AAA(00000000).csv'
raw = pyActigraphy.io.read_raw_rpx(fpath+filename)
raw_state = pyActigraphy.io.read_raw_rpx(fpath+filename)

instance=filename[0:13]
start=raw.start_time

'########################### reading file ###################################################'
f= open(fpath+filename, 'rt', encoding='UTF8')
reader = csv.reader(f)

csv_list=[]
for l in reader:
    csv_list.append(l)
f.close()
raw_sleepinfo = pd.DataFrame(csv_list)
raw_sleepinfo = pd.DataFrame(csv_list, columns=raw_sleepinfo.loc[64, :])

'####################################### masking #######################################'
mask = np.ones_like(raw.data)
raw.mask=pd.Series(mask, index=raw.data.index)
raw_state.mask=pd.Series(mask, index=raw.data.index)
raw.mask_inactivity=True
raw_state.mask_inactivity=True
nonwear =raw.data[raw.data.isna()].index
raw.mask[nonwear]=0
raw_state.mask[nonwear]=0

'####################### extract sleep time from actiware ###################################'
sleep = raw_sleepinfo.loc[raw_sleepinfo["Interval Type"] == "SLEEP", ["Start Date","Start Time", "End Date","End Time"]]
sleep.reset_index(drop=True, inplace=True)
sleepon=[]
sleepoff=[]

'''drop the invalid sleep'''
sleep.drop(index=set(np.where(sleep.values=='NaN')[0]), inplace=True) # 22/09/20 modified

'''sleep onset, offset resampling'''
if raw.frequency == pd.Timedelta('30S'):
    for i in range(len(sleep.values)):
        if sleep.values[i][1][-2] == '3':
            sleep.values[i][1] = sleep.values[i][1][0:-2] + '00'
    for j in range(len(sleep.values)):
        if sleep.values[j][3][-2] == '3':
            sleep.values[j][3] = sleep.values[j][3][0:-2] + '00'

if raw.frequency == pd.Timedelta('15S'):
    for i in range(len(sleep.values)):
        if sleep.values[i][1][-2] != '0':
            sleep.values[i][1] = sleep.values[i][1][0:-2] + '00'
    for j in range(len(sleep.values)):
        if sleep.values[j][3][-2] != '0':
            sleep.values[j][3] = sleep.values[j][3][0:-2] + '00'

for i in range(len(sleep.index)):
    sleep_onset = pd.to_datetime(' '.join([sleep.values[i][0], sleep.values[i][1]]))
    sleep_offset = pd.to_datetime(' '.join([sleep.values[i][2], sleep.values[i][3]]))
    sleepon.append(sleep_onset)
    sleepoff.append(sleep_offset)

for i in range(len(sleepon)):
    print(i, 'sleep onset:', sleepon[i], '   sleep offset:', sleepoff[i], '    sleep time:', sleepoff[i] - sleepon[i])


'######################################### resampling to 1 min #######################################'
i=0
if raw.frequency==pd.Timedelta('30S') or raw.frequency==pd.Timedelta('15S'):
    raw_state1 = raw_state.resampled_data(freq='1min')
    raw1 = raw.resampled_data(freq='1min')

'''######################################### sleep analysis #######################################'''
'''###### exclude nap and analyze only sleep time ########'''

sleep_act=[]
sleep_light=[]

for i in range(len(sleepon)):
    # The analysis only included sleep that began between 8:00 PM and 4:00 AM
    if sleep.values[i][1][0:2]>= '20' or sleep.values[i][1][0:2]<='04':
        temp_single_night = pyActigraphy.io.read_raw_rpx(fpath + filename, start_time=sleepon[i],
                                                         period=sleepoff[i] - sleepon[i])
        # masking
        mask = np.ones_like(temp_single_night.data)
        temp_single_night.mask = pd.Series(mask, index=temp_single_night.data.index)
        temp_single_night.mask_inactivity = True
        nonwear = temp_single_night.data[temp_single_night.data.isna()].index
        temp_single_night.mask[nonwear] = 0

        single_night_act = temp_single_night.resampled_data(freq='1min')
        single_night_light= temp_single_night.resampled_light(freq='1min')
        sleep_act.append(single_night_act)
        sleep_light.append(single_night_light)

'############ remove activity with sudden light in sleep ############'
# You can adjust these parameters
light_th = 1.5
light_baseline = 0.3
act_min_th= 200
act_min_th_2= 100
act_max_th=500

excluded_index=[]
sleep_act_mod_ver1=copy.deepcopy(sleep_act)
for i in range(len(sleep_light)):
    for j in range(len(sleep_light[i])-7):
        if j >=7:
            if (sleep_light[i][j-1] + sleep_light[i][j] + sleep_light[i][j+1] > light_th) & \
                    ((sleep_act[i][j - 1] + sleep_act[i][j] + sleep_act[i][j + 1] > act_min_th) or (sleep_act[i][j] > act_min_th_2)):
                if ((sleep_light[i][j - 7] < light_baseline) & (sleep_light[i][j + 7] < light_baseline)) or \
                        ((sleep_light[i][j - 3] < 0.2 * sleep_light[i][j]) & (sleep_light[i][j + 3] < 0.2 * sleep_light[i][j])) or\
                        ((sleep_light[i][j - 3] > light_th) & (sleep_light[i][j + 3] < max(0.3 * sleep_light[i][j], light_baseline))) or\
                        ((sleep_light[i][j + 3] > light_th) & (sleep_light[i][j - 3] < max(0.3 * sleep_light[i][j], light_baseline))):
                    excluded_index.append(sleep_act[i].index[j])
                    sleep_act_mod_ver1[i][j]=0
                    sleep_act_mod_ver1[i][j - 1] = 0
                    sleep_act_mod_ver1[i][j + 1] = 0

excluded_index_2=[]

sleep_act_mod_ver2=copy.deepcopy(sleep_act_mod_ver1)
for i in range(len(sleep_act_mod_ver1)):
    for j in range(len(sleep_act_mod_ver1[i]) - 1):
        if j >=1:
            if sleep_act_mod_ver1[i][j - 1]+sleep_act_mod_ver1[i][j]+sleep_act_mod_ver1[i][j + 1]>act_max_th:
                excluded_index_2.append(sleep_act_mod_ver1[i].index[j])
                sleep_act_mod_ver2[i][j]=0
                sleep_act_mod_ver2[i][j - 1] = 0
                sleep_act_mod_ver2[i][j + 1] = 0

'########################### active block detector ###########################'

print('-------------------------------------')
window_len=20  # moving average window
block_th=0.2  # block threshold, set to 20%
data=sleep_act_mod_ver2  # data
binarized_ratio_resized_list=[]
valid_act_count_min=20
valid_act_count_max=100

for s in range(len(sleep_act)):
    ratio_of_epochs_with_activity = []
    # calculate the ratio of epochs with the corresponding activity count value (20-100) in a 20-min window
    for i in range(len(data[s]) - window_len):
        counts_of_unique_values = data[s][i:i + window_len - 1].value_counts(normalize=True)

        if counts_of_unique_values.sort_index().index[0]<=valid_act_count_min or counts_of_unique_values.sort_index().index[-1]>=valid_act_count_max:
            # if there are epochs with activity counts outside the range of 20-100 activity counts,
            ratio_of_epochs_with_activity.append(1 - sum(counts_of_unique_values[counts_of_unique_values.index < valid_act_count_min]) - sum(counts_of_unique_values[counts_of_unique_values.index > valid_act_count_max]))
        else:  # if all epochs in the window have an activity count between 20-100
            ratio_of_epochs_with_activity.append(1)

    # smoothing ratio_of_epochs_with_activity using FLM
    freq = '10min'  # dummy value
    max_order = 15 # FLM order
    flm = FLM2(basis='fourier', sampling_freq=freq, max_order=max_order)
    flm._FLM2__nsamples = len(ratio_of_epochs_with_activity)
    X = np.stack(flm.basis_functions, axis=1)
    y = ratio_of_epochs_with_activity
    model = sm.OLS(y, X)
    results = model.fit()
    flm.beta['name'] = results.params
    y_est = np.dot(X, flm.beta['name'])  # smoothed ratio_of_epochs_with_activity

    # binarize the ratio (1:sleep, 0:wake)
    binarized_ratio = copy.deepcopy(list(y_est))
    for b in range(len(binarized_ratio)):
        if binarized_ratio[b] > block_th:  # if ration exceeds 20%
            binarized_ratio[b] = 1
        else:
            binarized_ratio[b] = 0

    # align the block and data sync by attaching dummies to the front and back of the shortened data
    front_dummy = [binarized_ratio[0]] * (window_len // 2)
    back_dummy = [binarized_ratio[-1]] * (window_len // 2)

    # reshape the list to match the size of sleep data
    binarized_ratio_resized = front_dummy + binarized_ratio + back_dummy

    # remove too short (<5min) block
    start_idx = 0
    for i in range(len(binarized_ratio_resized) - 1):
        if binarized_ratio_resized[i] != binarized_ratio_resized[i + 1]:
            if binarized_ratio_resized[i + 1] == 1:
                start_idx = i
            else:
                if (i - start_idx) < 5:
                    binarized_ratio_resized[start_idx + 1:i + 1] = [0 for x in range(i - start_idx)]

    # scaling the binarized_ratio_resized for plotting
    binarized_ratio_resized_for_plotting = [binarized_ratio_resized[i] * max(data[s]) for i in range(len(binarized_ratio_resized))]

    #find local maxima
    arr_ratio=np.array(ratio_of_epochs_with_activity)
    peaks, _ = find_peaks(arr_ratio, height=0.2, distance=30)

    # Attach dummies to the front and back of the shortened data
    dummy = int(window_len / 2) * [0]
    ratio = dummy + ratio_of_epochs_with_activity + dummy
    smoothed_ratio = dummy + list(y_est) + dummy

    # Replace the negative value of the smoothed ratio to 0
    for i in range(len(smoothed_ratio)):
        if smoothed_ratio[i] < 0:
            smoothed_ratio[i] = 0

    peaks_for_plot = [x + len(dummy) for x in peaks]

    # Process figure for each sleep
    ratio_percentage = [100 * i for i in ratio]
    smoothed_ratio_percentage = [100 * i for i in smoothed_ratio]
    fig, ax1 = plt.subplots()
    ax1.set_ylim(-4.7, 100)
    ax1.plot(ratio_percentage, color='green', alpha=0.5)
    ax1.plot(smoothed_ratio_percentage, color='slategrey')
    ax1.plot(peaks_for_plot, np.array(ratio_percentage)[peaks_for_plot], '*', color='darkgreen')
    plt.title('Sleep' + str(s) + ': Processed and smoothed signal')
    plt.ylabel('Ratio of epochs with activity (%)', fontsize=12)
    plt.xlabel('Time (min)', fontsize=12)
    # plt.axhline(y=20, linewidth=1, color='darkgreen', linestyle='--')

    # Block figure for each sleep
    fig, ax1 = plt.subplots()
    ax1.plot(data[s].values, alpha=0.7)
    ax1.plot(binarized_ratio_resized_for_plotting, color='C1')
    ax1.fill(binarized_ratio_resized_for_plotting, color='peachpuff', alpha=0.8)
    plt.title('Sleep' + str(s) + ': Raw signal + block')
    plt.ylabel('Activity counts', fontsize=12)
    plt.xlabel('Time (min)', fontsize=12)

    binarized_ratio_resized_list.append(binarized_ratio_resized)

total_inactive_percentage=[]  # inactive %
total_num_block=[]  # number of block
total_interval_list=[]  # interval of block
total_duration_list=[]  # duration of block
total_onset_list=[]  # block onset time
total_offset_list=[]  # block offset time

for k in range(len(binarized_ratio_resized_list)):
    binarized_ratio_resized_series = pd.Series(data=binarized_ratio_resized_list[k], index=sleep_act_mod_ver2[k].index)

    interval_list = []
    duration_list = []
    onset_list = []
    offset_list = []
    num_block=0

    print('Sleep', k+1)
    j = binarized_ratio_resized_series.index[0]  # start point
    for i in range(len(binarized_ratio_resized_series) - 1):
        if binarized_ratio_resized_series[i] != binarized_ratio_resized_series[i + 1]:
            if binarized_ratio_resized_series[i + 1] == 1:
                if len(interval_list)==0:
                    print('Sleep-to-1st block interval: ', binarized_ratio_resized_series.index[i] - j)
                else:
                    print('Block-to-block interval: ', binarized_ratio_resized_series.index[i] - j)
                print()
                interval_list.append(binarized_ratio_resized_series.index[i] - j)
                total_interval_list.append(binarized_ratio_resized_series.index[i] - j)
            elif binarized_ratio_resized_series[i + 1] == 0:
                num_block += 1
                print('Block', num_block)
                print('Start point: ', j)
                print('Stop point: ', binarized_ratio_resized_series.index[i])
                print('Duration: ', binarized_ratio_resized_series.index[i] - j)
                print()

                onset_list.append(j)
                offset_list.append(binarized_ratio_resized_series.index[i])
                duration_list.append(binarized_ratio_resized_series.index[i] - j)
                total_onset_list.append(j)
                total_offset_list.append(binarized_ratio_resized_series.index[i])
                total_duration_list.append(binarized_ratio_resized_series.index[i] - j)
            j = binarized_ratio_resized_series.index[i] + pd.Timedelta('1min')

        elif i== (len(binarized_ratio_resized_series) - 2):
            if binarized_ratio_resized_series[i]==0:
                print('Block-to-wake interval: ', binarized_ratio_resized_series.index[i + 1] - j)
                print()
                interval_list.append(binarized_ratio_resized_series.index[i + 1] - j)
                total_interval_list.append(binarized_ratio_resized_series.index[i + 1] - j)
            elif binarized_ratio_resized_series[i]==1:
                num_block += 1
                print('Block', num_block)
                print('Start point: ', j)
                print('Stop point: ', binarized_ratio_resized_series.index[i + 1])
                print('Duration: ', binarized_ratio_resized_series.index[i + 1] - j)
                print()

                onset_list.append(j)
                offset_list.append(binarized_ratio_resized_series.index[i + 1])
                duration_list.append(binarized_ratio_resized_series.index[i + 1] - j)
                total_onset_list.append(j)
                total_offset_list.append(binarized_ratio_resized_series.index[i + 1])
                total_duration_list.append(binarized_ratio_resized_series.index[i + 1] - j)


    inactive_percentage = binarized_ratio_resized_list[k].count(0) / len(binarized_ratio_resized_list[k])*100
    total_inactive_percentage.append(inactive_percentage)
    total_num_block.append(num_block)

    print('### Motor Activity Block summary for sleep', k+1)
    print('Inactive %: ', round(inactive_percentage), '%')
    print('Mean duration: ', np.mean(duration_list))
    print('Mean interval: ', np.mean(interval_list))
    print('Mean number of blocks: ', num_block)
    print('-------------------------------------')
    print()


print('Average MAB indices over', len(sleep_act), 'days')
print('Mean inactive %: ', round(np.mean(total_inactive_percentage)), '%')
print('Mean duration: ', np.mean(total_duration_list).floor('s'))
print('Mean interval: ', np.mean(total_interval_list).floor('s'))
print('Mean number of blocks: ', round(np.mean(total_num_block), 2))

print('--------------------------------------')
print('Inactive% Duration Interval Number')
print(round(np.mean(total_inactive_percentage)), round(np.mean(total_duration_list).seconds / 60, 2), round(np.mean(total_interval_list).seconds / 60, 2), round(np.mean(total_num_block), 2))

for i in range(len(total_interval_list)):
    total_interval_list[i]=total_interval_list[i].seconds/60
for i in range(len(total_duration_list)):
    total_duration_list[i] = total_duration_list[i].seconds / 60

'Histograms for the 4 MAB indices'
plt.figure()
plt.title('Inactive % histogram')
plt.hist(total_inactive_percentage, bins=10, range=[60, 100])

plt.figure()
plt.title('Block duration histogram')
plt.hist(total_duration_list, bins=10, range=[0,180])

plt.figure()
plt.title('Block interval histogram')
plt.hist(total_interval_list, bins=10, range=[0,180])

plt.figure()
plt.title('Number of block histogram')
plt.hist(total_num_block, bins=10, range=[0, 15])


plt.show()