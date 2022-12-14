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

sleep_num=[]
sleep_light=[]

for i in range(len(sleepon)):
    if sleep.values[i][1][0:2]>= '20' or sleep.values[i][1][0:2]<='04':
        test = pyActigraphy.io.read_raw_rpx(fpath + filename, start_time=sleepon[i],
                                            period=sleepoff[i] - sleepon[i])
        mask = np.ones_like(test.data)
        test.mask = pd.Series(mask, index=test.data.index)
        test.mask_inactivity = True
        nonwear = test.data[test.data.isna()].index
        test.mask[nonwear] = 0

        test1 = test.resampled_data(freq='1min')
        light1= test.resampled_light(freq='1min')
        sleep_num.append(test1)
        sleep_light.append(light1)

'############ remove activity with sudden light in sleep ############'

light_th = 1.5
light_min = 0.3
act_min= 200

exc=[]
sleep_num_exc1=copy.deepcopy(sleep_num)
for i in range(len(sleep_light)):
    for j in range(len(sleep_light[i])-7):
        if j >=7:
            if (sleep_light[i][j-1] + sleep_light[i][j] + sleep_light[i][j+1] > light_th) & \
                    ((sleep_num[i][j-1]+sleep_num[i][j]+sleep_num[i][j+1]>act_min) or (sleep_num[i][j]>100)):
                if ((sleep_light[i][j - 7] < light_min) & (sleep_light[i][j + 7] < light_min)) or \
                        ((sleep_light[i][j - 3] < 0.2 * sleep_light[i][j]) & (sleep_light[i][j + 3] < 0.2 * sleep_light[i][j])) or\
                        ((sleep_light[i][j - 3] > light_th) & (sleep_light[i][j + 3] < max(0.3 * sleep_light[i][j],light_min))) or\
                        ((sleep_light[i][j + 3] > light_th) & (sleep_light[i][j - 3] < max(0.3 * sleep_light[i][j],light_min))):
                    exc.append(sleep_num[i].index[j])
                    sleep_num_exc1[i][j]=0
                    sleep_num_exc1[i][j -1] = 0
                    sleep_num_exc1[i][j +1] = 0

exc2=[]
maxth=500

sleep_num_exc=copy.deepcopy(sleep_num_exc1)
for i in range(len(sleep_num_exc1)):
    for j in range(len(sleep_num_exc1[i])-1):
        if j >=1:
            if sleep_num_exc1[i][j-1]+sleep_num_exc1[i][j]+sleep_num_exc1[i][j+1]>maxth:
                exc2.append(sleep_num_exc1[i].index[j])
                sleep_num_exc[i][j]=0
                sleep_num_exc[i][j -1] = 0
                sleep_num_exc[i][j +1] = 0

'############################# active block detector #############################'

print('-------------------------------------')
window=20
th=0.2
data=sleep_num_exc
rebinper0_list=[]
ref=20
refmax=100

for s in range(len(sleep_num)):
    per0 = []
    # calculate the percent of 0 count within the window
    for i in range(len(data[s]) - window):
        a = data[s][i:i + window - 1].value_counts(normalize=True)
        # print(a[0])
        if a.sort_index().index[0]<=ref or a.sort_index().index[-1]>=refmax:
            per0.append(1-sum(a[a.index<ref])-sum(a[a.index>refmax]))
            #per0.append(a[0])
        else:
            #per0.append(1)
            per0.append(1)

    #smoothing per0 using FLM
    freq = '10min'  # dummy
    max_order = 15
    flm = FLM2(basis='fourier', sampling_freq=freq, max_order=max_order)
    flm._FLM2__nsamples = len(per0)
    X = np.stack(flm.basis_functions, axis=1)
    y = per0
    model = sm.OLS(y, X)
    results = model.fit()
    flm.beta['name'] = results.params
    y_est = np.dot(X, flm.beta['name'])

    # binarize the percent of 0 count, 1:sleep, 0:wake
    binper0 = copy.deepcopy(per0)
    binper0 = copy.deepcopy(list(y_est))
    for b in range(len(binper0)):
        if binper0[b] > th:
            binper0[b] = 1
            # 1
        else:
            binper0[b] = 0
            # 0


    temp1 = [binper0[0]] * (window // 2)
    temp2 = [binper0[-1]] * (window // 2)

    # reshape the list to match the size of sleep data
    rebinper0 = temp1 + binper0 + temp2

    # remove <5min active block
    j = 0
    for i in range(len(rebinper0) - 1):
        if rebinper0[i] != rebinper0[i + 1]:
            if rebinper0[i + 1] == 1:
                j = i  # ??????????????? 0?????? ?????????
            else:
                if (i - j) < 5:
                    rebinper0[j + 1:i + 1] = [0 for x in range(i - j)]

    # scaling the list for plotting
    scbinper0 = [rebinper0[i] * max(data[s]) for i in range(len(rebinper0))]

    #find local maxima
    nper0=np.array(per0)
    peaks, _ = find_peaks(nper0, height=0.2, distance=30)

    # FLM ?????? (Y_est)?????? ?????? ??? 0?????? ??????
    # ????????? ??????????????? ????????? ????????? ????????? ?????? ?????????(0) ?????????(Padding ??? ??????)
    scaled_per0 = per0
    scaled_yest = list(y_est)
    dummy = int(window / 2) * [0]

    for i in range(len(scaled_yest)):
        if scaled_yest[i] < 0:
            scaled_yest[i] = 0
    scaled_per0 = dummy + scaled_per0 + dummy
    scaled_yest = dummy + scaled_yest + dummy
    scaled_peaks = [x+len(dummy) for x in peaks]

    # plotting
    # processed by sliding window & smoothing (FLM)
    scaled_per0_perc = [100 * i for i in scaled_per0]
    scaled_yest_perc = [100 * i for i in scaled_yest]
    fig, ax1 = plt.subplots()
    ax1.set_ylim(-4.7, 100)
    ax1.plot(scaled_per0_perc, color='green', alpha=0.5)
    ax1.plot(scaled_yest_perc, color='slategrey')
    ax1.plot(scaled_peaks, np.array(scaled_per0_perc)[scaled_peaks], '*', color='darkgreen')
    plt.title('Sleep' + str(s) + ': Processed and smoothed signal')
    plt.ylabel('Ratio of epochs with activity (%)', fontsize=12)
    plt.xlabel('Time (min)', fontsize=12)
    # plt.axhline(y=20, linewidth=1, color='darkgreen', linestyle='--')

    # Block detection
    fig, ax1 = plt.subplots()
    ax1.plot(data[s].values, alpha=0.7)
    ax1.plot(scbinper0, color='C1')
    ax1.fill(scbinper0, color='peachpuff', alpha=0.8)
    plt.title('Sleep' + str(s) + ': Raw signal + block')
    plt.ylabel('Activity counts', fontsize=12)
    plt.xlabel('Time (min)', fontsize=12)

    rebinper0_list.append(rebinper0)

total_activeper=[]
total_activeblock=[]
total_interval_list=[]
total_duration_list=[]
total_onset_list=[]
total_offset_list=[]

for k in range(len(rebinper0_list)):
    rebinper0_ser = pd.Series(data=rebinper0_list[k], index=sleep_num_exc[k].index)
    timevoting = pd.Series(data=[0] * len(rebinper0_list[k]), index=sleep_num_exc[k].index)

    interval_list = []
    duration_list = []
    onset_list = []
    offset_list = []
    activeblock=0

    print('Sleep', k+1)
    j = rebinper0_ser.index[0]  # start point
    for i in range(len(rebinper0_ser) - 1):
        if rebinper0_ser[i] != rebinper0_ser[i + 1]:
            if rebinper0_ser[i + 1] == 1:
                if len(interval_list)==0:
                    print('Sleep-to-1st block interval: ', rebinper0_ser.index[i] - j)
                else:
                    print('Block-to-block interval: ', rebinper0_ser.index[i] - j)
                print()
                interval_list.append(rebinper0_ser.index[i] - j)
                total_interval_list.append(rebinper0_ser.index[i] - j)
            elif rebinper0_ser[i + 1] == 0:
                activeblock += 1
                print('Block',activeblock)
                print('Start point: ', j)
                print('Stop point: ', rebinper0_ser.index[i])
                print('Duration: ', rebinper0_ser.index[i] - j)
                print()

                onset_list.append(j)
                offset_list.append(rebinper0_ser.index[i])
                duration_list.append(rebinper0_ser.index[i] - j)
                total_onset_list.append(j)
                total_offset_list.append(rebinper0_ser.index[i])
                total_duration_list.append(rebinper0_ser.index[i] - j)
            j = rebinper0_ser.index[i] + pd.Timedelta('1min')

        elif i== (len(rebinper0_ser) - 2):
            if rebinper0_ser[i]==0:
                print('Block-to-wake interval: ', rebinper0_ser.index[i+1] - j)
                print()
                interval_list.append(rebinper0_ser.index[i+1] - j)
                total_interval_list.append(rebinper0_ser.index[i+1] - j)
            elif rebinper0_ser[i]==1:
                activeblock += 1
                print('Block',activeblock)
                print('Start point: ', j)
                print('Stop point: ', rebinper0_ser.index[i+1])
                print('Duration: ', rebinper0_ser.index[i+1] - j)
                print()

                onset_list.append(j)
                offset_list.append(rebinper0_ser.index[i+1])
                duration_list.append(rebinper0_ser.index[i+1] - j)
                total_onset_list.append(j)
                total_offset_list.append(rebinper0_ser.index[i+1])
                total_duration_list.append(rebinper0_ser.index[i+1] - j)

    # for i in range(len(onset_list)):
    #    timevoting[onset_list[0]:offset_list[0]]+=1

    activeper = rebinper0_list[k].count(0) / len(rebinper0_list[k])
    total_activeper.append(activeper)
    total_activeblock.append(activeblock)

    print('### Motor Activity Block summary for one night ###')
    print('Inactive %: ', activeper)
    print('Mean duration: ', np.mean(duration_list))
    print('Mean interval: ', np.mean(interval_list))
    print('Mean number of block: ', activeblock)
    print('-------------------------------------')
    print()


print('Total characteristics')
print('Total inactive %: ', round(np.mean(total_activeper),2))
print('Total mean duration: ', np.mean(total_duration_list))
print('Total mean interval: ', np.mean(total_interval_list))
print('Total mean number of active block: ', round(np.mean(total_activeblock),2))

print('--------------------------------------')
print('Inactive% Duration Interval Number')
print(round(np.mean(total_activeper),2), round(np.mean(total_duration_list).seconds/60,2), round(np.mean(total_interval_list).seconds/60,2), round(np.mean(total_activeblock),2))
#print(round(np.mean(total_duration_list).seconds/60,2))
#print(round(np.mean(total_interval_list).seconds/60,2))
#print(round(np.mean(total_activeblock),2))

for i in range(len(total_interval_list)):
    total_interval_list[i]=total_interval_list[i].seconds/60
for i in range(len(total_duration_list)):
    total_duration_list[i] = total_duration_list[i].seconds / 60


plt.figure()
plt.title('inactive % histogram')
plt.hist(total_activeper,bins=10, range=[0,1])

plt.figure()
plt.title('block duration histogram')
plt.hist(total_duration_list, bins=10, range=[0,180])

plt.figure()
plt.title('block interval histogram')
plt.hist(total_interval_list, bins=10, range=[0,180])

plt.figure()
plt.title('block number histogram')
plt.hist(total_activeblock, bins=10, range=[0,15])


plt.show()