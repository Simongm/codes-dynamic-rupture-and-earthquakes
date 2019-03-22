# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:46:33 2017

@author: Snhh26
"""
import os,sys
import matplotlib.pyplot as plt
import numpy as np


## ADD THE PATH TO FIND AND IMPORT THE READ FUNCTION ####
sys.path.append("C:\\Users\\Snhh26\\Desktop\\DU965_9June2017\\functions")
from read_function import read_function

##SET THE WORKING DIRECTORY: ######################################
os.chdir('C:\\Users\\Snhh26\\Desktop\\DU965_9June2017')
print(os.getcwd()) # to check working directory where the files are

################################parameters to change#################################
#richter='richter-m' #name of richter 
first=5 #number of 1st richter file
last=5 #number of last richter-file
### if we want to display  only a few datasets:
dataset1=1 # 1st dataset
dataset2='all' #last dataset, write dataset2='all' if you want to display everything

## if you want to take only few points in the dataset, otherwise leave default values for 1MHz
firstindex=0
lastindex=4980736
cuto=10000000 # low pass cutoff frequency

ratio=1 ## decimating ration, chose a multiple of 2 : [1,2,4,7,8,14,16,28,32,56,64,112,128,224,256,448,512,896,1024,1792,2048,3584,4096,7168,8192,14336,16384,28672,32768,57344,65536,114688,229376,458752]
filtercondition=0 # 0= no filter; 1=lowpass at cuto freq
indochannels=[0] #indices of channels to display, =[0,1,2,3] if 4 channels



[data_m,time_m]=read_function('richter-m',first,last,dataset1,dataset2,firstindex,lastindex,cuto,ratio,filtercondition,indochannels)
[data_s1,time_s1]=read_function('richter-s1',first,last,dataset1,dataset2,firstindex,lastindex,cuto,ratio,filtercondition,indochannels)
#
#example of plot
plt.plot(time_m,data_m)
plt.xlabel('time in s')
plt.ylabel('voltage in V')

## example of event counter on 1 channel..........................
threshold=0.5
Window=0.01

cond=data_m>0.5 ## find all the points with amplitude > 0.5V
indall=(np.where(cond))[0]

event=data_m[indall]
WindowL=200 # window length to adjust, not to pick the same event

## loop to remove points in a same event
t0=0
for i in range (0,len(indall)):
    t1=time_m[indall[i]]
    if t1-t0<0.005: 
        indall[i]=0
    t0=t1
    
plt.plot(time_m[indall],data_m[indall],'r.')


#time_event=[]
#i=0
#for time_ev in time_m:
#    data_ev=data_m[i]
#    if time_ev>
#    time_event=np.append(time_event,time_ev)
#    
#    i=i+1



### part to translate:
    
############################################################################
# % Hit count from Richter data files at 1kHz sampling rate
############################################################################
#tic
#clear
#
#% import wve file and get length of SRM file
#datafile_master = 'G:\Stephan\COSx-1-59_MHF\COSx-1-59_MHF_Richter\raw data\\rct-uop-090414.data.00006'; % file path and file name (excluding file extension) of richter time file containing failure
#datafile_slave1 = 'G:\Stephan\COSx-1-59_MHF\COSx-1-59_MHF_Richter\raw data\\rct-uop-180113.data.00006'; % file path and file name (excluding file extension) of richter time file containing failure
#datafile_slave2 = 'G:\Stephan\COSx-1-59_MHF\COSx-1-59_MHF_Richter\raw data\\rct-uop-021115.data.00006'; % file path and file name (excluding file extension) of richter time file containing failure
#
#data = importdata (char([datafile_master '.wve']));
#time_start = data.data(10);
#time_finish = data.data(11);
#hhmmss_start = num2str(time_start);
#hhmmss_finish = num2str(time_finish);
#start_time_ss = str2double(hhmmss_start(1:2))*3600 + str2double(hhmmss_start(3:4))*60 + str2double(hhmmss_start(5:end));
#length_srm = (str2double(hhmmss_finish(1:2))*3600 + str2double(hhmmss_finish(3:4))*60 + str2double(hhmmss_finish(5:end)))-...
#    (str2double(hhmmss_start(1:2))*3600 + str2double(hhmmss_start(3:4))*60 + str2double(hhmmss_start(5:end)));
#last_step = (floor(length_srm*10) -1.000001)/10;
#rest = ((length_srm*10 - floor(length_srm*10))/10)+0.000001;
#
#% Read 16 bit, 4 channel SRM data
#
#%import Master Richter first 0.1s
#
#%namefile=sprintf('%s-%s-%s',datafile(16:18),datafile(20:21),datafile(23:26));
#namefile = 'data';
#fname_bin1=char([datafile_master '.srm']);
#fname_inx=char([datafile_master '.wve']);
#Fs=10E6; % sampling frequency
#L = 1E4; # time window
#% open binary data file, read data
#fid = fopen(fname_bin1,'r');
#secs=0; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 10
#HPfilter = 4e4;
#LPfilter = 1e6;
#d = fdesign.bandpass('N,F3dB1,F3dB2',10, HPfilter,LPfilter,Fs); 
#Hd = design(d,'butter');
#data1 = filter (Hd,data);
#
#% Import Slave 1 Richter first 0.1s
#
#%eval(sprintf('datafile=''G:\\Stephan\\%s%s-%s-%.0f_HF\\%s%s-%s-%.0f_HF_Richter\\raw data\\rct-uop-180113.data.00%s'';',rock,dir,block,k,rock,dir,block,k,a))
#%namefile=sprintf('%s-%s-%s',datafile(16:18),datafile(20:21),datafile(23:26));
#namefile = 'data';
#fname_bin2=char([datafile_slave1 '.srm']);
#fname_inx=char([datafile_slave1 '.wve']);
#Fs=10E6; % sampling frequency
#L = 1E4; % time window
#% open binary data file, read data
#fid = fopen(fname_bin2,'r');
#secs=0; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 10
#data2 = filter (Hd,data);
#
#% Import second Slave 2 Richter first 0.1s
#
#%eval(sprintf('datafile=''G:\\Stephan\\%s%s-%s-%.0f_HF\\%s%s-%s-%.0f_HF_Richter\\raw data\\rct-uop-021115.data.00%s'';',rock,dir,block,k,rock,dir,block,k,a))
#%namefile=sprintf('%s-%s-%s',datafile(16:18),datafile(20:21),datafile(23:26));
#namefile = 'data';
#fname_bin3=char([datafile_slave2 '.srm']);
#fname_inx=char([datafile_slave2 '.wve']);
#Fs=10E6; % sampling frequency
#L = 1E4; % time window
#% open binary data file, read data
#fid = fopen(fname_bin3,'r');
#secs=0; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data3 = filter (Hd,data);
#%get absolute values of the waveform
#data = [abs(data1),abs(data2),abs(data3)]; 
#%set threshold for each channel
#threshold = [4 0.01 0.01 0.01 0.01 0.01 0.01 0.0115 0.0315 0.01 0.01 0.01]; 
#
#% find points of relative maxima on each channel
#
##################################################################################################
#% Alternate Version superfast by Eduardo Rossi
#% Pay attention to these indices, there is I was wrong. Recheck.
#% The idea is to take two subsets of the data matrix
#% Initial and make the differences between these sub-matrices. A eye
#% They should have a size SIZE-2 to perform what
#% first did in for loops. Recheck making a silly example with
#% an initial array of, for example, 5 x 8
#
#derivative_1 = data(2:end-1,:)-data(1:end-2,:);
#derivative_2 = data(3:end,:)-data(2:end-1,:);
#
#dev_1_pos = derivative_1>0;
#dev_2_pos = derivative_2<0;
#
#clear derivative_1
#clear derivative_2
#
#matrix_signs = dev_1_pos+dev_2_pos; 
#% If you write well wherever it appears a "2" is equivalent to having derivata1 greater than zero and less than zero derivata2 
#clear dev_1_pos
#clear dev_2_pos
#
#% Find the matrix elements that have the value 2, which
#% Means that they have the condition sought
#[elementi_i, elementi_j] = find(matrix_signs==2);
#clear matrix_signs
#% Now you have to be careful, because © indexes found actually they are not
#% Those of the initial matrix, but are those of the sub-matrices, so you have to correct one or two indices
#% With the silly example of an initial matrix 5 x 8 you see it as easy as
#% Correct. This correction I call for now "CORR", see you
#% quant'Ã¨
#
#% If the found items are not many, you can use a for,
#% Remembering to preallocate ALWAYS the object to make for (that
#% Greatly increases the speed of for loops in Matlab
#hits = zeros(1, 12);  % 12 channels
#ampl = linspace(0.01,2.3,20); % CREATE A 3D MATRIX, WITH SECONDS, CLUSTERS AND CHANNELS AS DIMENSIONS
#for i=1:length(elementi_j)  
#    if data(elementi_i(i) + 1, elementi_j(i)) >= threshold(elementi_j(i));
#        hits(1,elementi_j(i)) = hits(1,elementi_j(i))+1;  
#    end
#end
#
#time_step = 0.0999999:0.001:last_step;
#for secs = time_step; % add 0.1 minus 0.1us each time
#%import Master Richter second 0.1s
#L = 1E4+1; % time window (plus last point of the previous window time)
#% open binary data file, read data
#fid = fopen(fname_bin1,'r');
#secs=secs; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data1 = filter (Hd,data);
#
#% Import Slave 1 Richter second 0.1s 
#L = 1E4+1; % time window
#% open binary data file, read data
#fid = fopen(fname_bin2,'r');
#secs=secs; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data2 = filter (Hd,data);
#
#% Import second Slave 2 Richter second 0.1s 
#L = 1E4+1; % time window
#% open binary data file, read data
#fid = fopen(fname_bin3,'r');
#secs=secs; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data3 = filter (Hd,data);
#data = [abs(data1),abs(data2),abs(data3)]; %get absolute values of the waveform
#derivative_1 = data(2:end-1,:)-data(1:end-2,:);
#derivative_2 = data(3:end,:)-data(2:end-1,:);
#
#dev_1_pos = derivative_1>0;
#dev_2_pos = derivative_2<0;
#
#clear derivative_1
#clear derivative_2
#
#matrix_signs = dev_1_pos+dev_2_pos;
#clear dev_1_pos
#clear dev_2_pos
#[elementi_i, elementi_j] = find(matrix_signs==2);
#clear matrix_signs
#hits1 = zeros(1, 12);
#
#for i=1:length(elementi_j)  
#    if data(elementi_i(i) + 1, elementi_j(i)) >= threshold(elementi_j(i));
#        hits1(1,elementi_j(i)) = hits1(1,elementi_j(i))+1;  
#    end
#end
#hits = [hits;hits1];
#end
#
#% get hits of the incomplete 0.1 second
#%import Master Richter rest time
#L = rest*10E6; % time window (plus last point of the previous window time)
#% open binary data file, read data
#fid = fopen(fname_bin1,'r');
#secs=length_srm - rest; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data1 = filter (Hd,data);
#
#% Import Slave 1 Richter rest time
#L = rest*10E6; % time window
#% open binary data file, read data
#fid = fopen(fname_bin2,'r');
#secs=length_srm - rest; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data2 = filter (Hd,data);
#
#% Import second Slave 2 Richter rest time
#L = rest*10E6; % time window
#% open binary data file, read data
#fid = fopen(fname_bin3,'r');
#secs=length_srm - rest; %number of seconds to skip
#skipNbytes=secs*Fs*2*4; %4 bytes per int32, 2 bytes per int16, 4samples per time step (multiplexed)
#fseek(fid, skipNbytes, 'bof'); % skip secs seconds
#data = fread(fid,[4,L],'uint16'); %1E6=0.1s
#data = ((((data/2^16)*5)-2.5)*2)';
#fclose(fid);
#%apply bandpass Butterworth filter of order 100
#data3 = filter (Hd,data);
#data = [abs(data1),abs(data2),abs(data3)]; %get absolute values of the waveform
#derivative_1 = data(2:end-1,:)-data(1:end-2,:);
#derivative_2 = data(3:end,:)-data(2:end-1,:);
#
#dev_1_pos = derivative_1>0;
#dev_2_pos = derivative_2<0;
#
#clear derivative_1
#clear derivative_2
#
#matrix_signs = dev_1_pos+dev_2_pos;
#clear dev_1_pos
#clear dev_2_pos
#[elementi_i, elementi_j] = find(matrix_signs==2);
#clear matrix_signs
#hits1 = zeros(1, 12);
#
#
#for i=1:length(elementi_j)  
#    if data(elementi_i(i) + 1, elementi_j(i)) >= threshold(elementi_j(i));
#        hits1(1,elementi_j(i)) = hits1(1,elementi_j(i))+1;  
#    end
#end
#hits = [hits;hits1];
#hits = [zeros(length(hits(:,1)),1) hits];
#for i = 1:length (hits(:,1))
#    hits(i,1) = (start_time_ss*1000+i)/1000;
#end
#% get hit rate count per second
#hits_second =[];
#for i=1:floor(length_srm) % complete seconds in the file
#    hits_second = [hits_second;sum(hits(i*10-9:i*10,2:end))];
#end
#hits_second = [hits_second; sum(hits(floor(length_srm)*10+1:end,2:end))];
#hits_second = [zeros(floor(length_srm)+1,1),hits_second];
#for i=1:floor(length_srm)+1
#    hits_second(i,1) = i+start_time_ss;
#end
#
#eval(sprintf('save(''MHF_hitsx1kHz.txt'',''hits'',''-ascii'')'));
#eval(sprintf('save(''MHF_hitsx1Hz.txt'',''hits_second'',''-ascii'')'));
#toc
