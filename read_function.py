#!/usr/bin/python
import h5py
from scipy import signal
import numpy as np

def cutoff_filt(data,order,freqs,Samplingrate,ratio):
    CutOffFreq_l=freqs[0]/ratio
    SamplingDec=Samplingrate/ratio
    Wn_l=np.float64(CutOffFreq_l)/(SamplingDec/2);
    CutOffFreq_h=freqs[1]/ratio
    SamplingDec=Samplingrate/ratio
    Wn_h=np.float64(CutOffFreq_h)/(SamplingDec/2);
    b, a = signal.butter(order, [Wn_l,Wn_h],'bandstop')
    data = signal.filtfilt(b, a, np.transpose(data))
    data=np.transpose(data)
    return data

def read_function(richter,first,last,dataset1,dataset2,firstindex,lastindex,cuto,ratio,filtercondition,indochannels):

    channels=np.array(indochannels)

    ##############################################################################
    Files=np.arange(first,last+1,1) #syntax np.arange(2nd file,last file+1,1)
    Number_files=np.size(Files, axis=0)
    #initialize the array with the 1st file
    filnum='{0:05}'.format(Files[0])
    file1 = h5py.File('.'.join([richter,'data',filnum,'h5']))
    ##other files:
    print('Files',Files)
    print(file1)
    print('Number_files',Number_files)
    print('channels',file1.attrs['Channels Acquired']+1)
    Samplingrate=file1.attrs['Sampling Rate (Hz)']
    print('Sampling rate (Hz)',Samplingrate)   
    #initialize the datasets                      
    datasetlist1=list(file1.keys())
#    print('datasets',datasetlist1)
    nbdatasets_file1=np.size(datasetlist1, axis=0)
    if dataset2=="all":
        dataset2=len(datasetlist1)
        print(dataset2)
    #nbdatasets_file1=10 ######just to check->comment it otherwise
    datasetname='{0:08}'.format(dataset1)
    ## initialize the concatenated data
#    firstindex=0
#    lastindex=4980736
    index=np.array(np.arange(firstindex,lastindex,ratio))
    
    alldata=(np.array(file1[datasetname])[:,channels])[index,:]
    print('1st file','dataset number',datasetname,'of',nbdatasets_file1)
    
    # concatenate data from 1st file
    #for i in range(1,nbdatasets_file1+1):
    #     if only few datasets to select in the 1st file:
    for i in range(dataset1+1,dataset2+1):
       datasetname='{0:08}'.format(i)
#       dataset=file1[datasetname]
       alldata=np.concatenate((alldata,(np.array(file1[datasetname])[:,channels])[index,:]))
       
       print('1st file','dataset number',datasetname,'of',nbdatasets_file1)
    
       
    ## concatenate data from other files if more than one.........................
    if Number_files>1:
    
        filnum='{0:05}'.format(Files[1])
        first_file2 = h5py.File('.'.join([richter,'data',filnum,'h5']),'r')
        list_of_files='.'.join([richter,'data',filnum,'h5'])
        ## initialize the concatenated data
        datasetlist=first_file2.keys() #nb of datasets in the 1st file
    
        print('list of files',list_of_files) 
    
        ## create a list of remaining files to concatenate
        for n in range(1,Number_files):
            filnum='{0:05}'.format(Files[n])
            currentfilename='.'.join([richter,'data',filnum,'h5'])
            list_of_files=np.append(list_of_files,currentfilename)
                
       ## from the list of files, concatenate all the datasets
        for i in range(1,Number_files):
            print('list_of_files[i]',list_of_files[i])
            currentfile=h5py.File(list_of_files[i])
            datasetlist=list(currentfile.keys())
            nbdatasets=np.size(datasetlist, axis=0)
            
            for j in range(0,nbdatasets):
                datasetname='{0:08}'.format(j+1)
#                dataset=currentfile[datasetname]
                print('file number',i+1,'of',Number_files,'dataset number',j+1,'of',nbdatasets)
                alldata=np.concatenate((alldata,(np.array(currentfile[datasetname])[:,channels])[index,:]))
            
        
    npts=np.size(alldata, axis=0)
    print('length data', npts)
    print(alldata.shape)
    dt=np.float(ratio)/Samplingrate;
    timevector=np.arange(0,npts*dt,dt)
    alldata=(alldata/6553.6-5) #as read
    print('dt',dt)
    print('last time',timevector[-1])
#    plt.plot(timevector,alldata)
#    plt.show
    
    # filter:
    if  filtercondition > 0:
##lowpass filter     
        CutOffFreq_L=cuto/ratio
        SamplingDec=Samplingrate/ratio
        Wn_L=np.float64(CutOffFreq_L)/(SamplingDec/2);
        b, a = signal.butter(2, Wn_L)
        alldata = signal.filtfilt(b, a, np.transpose(alldata))
        alldata=np.transpose(alldata)

### bandstop filter
#    if  filtercondition == 2: #as periodic noise is in this frequency range
#        alldata=cutoff_filt(alldata,4,[390000,550000],Samplingrate,ratio) 
#
### median filter
#    if  filtercondition == 3: #as periodic noise is in this frequency range
#        alldata=signal.medfilt(np.transpose(alldata),25)
#        alldata=np.transpose(alldata)
#        
        
    return [alldata,timevector]


