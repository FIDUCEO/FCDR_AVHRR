import sys
import numpy as np
import matplotlib.pyplot as plt

class read_data(object):

    def set_var(self,noise_type,channel,d):

        if noise_type == 'independent':
            ntype = 1
        elif noise_type == 'structured':
            ntype = 2
        elif noise_type == 'measurement':
            ntype = 3

        if channel == 'Ch1':
            chan = 1
        elif channel == 'Ch2':
            chan = 2
        elif channel == 'Ch3a':
            chan = 3
        elif channel == 'Ch3b':
            chan = 4
        elif channel == 'Ch4':
            chan = 5
        elif channel == 'Ch5':
            chan = 6

        if 1 == ntype:
            if 1 == chan:
                self.independent_ch1 = d
            elif 2 == chan:
                self.independent_ch2 = d
            elif 3 == chan:
                self.independent_ch3a = d
            elif 4 == chan:
                self.independent_ch3b = d
            elif 5 == chan:
                self.independent_ch4 = d
            elif 6 == chan:
                self.independent_ch5 = d
        elif 2 == ntype:
            if 1 == chan:
                self.structured_ch1 = d
            elif 2 == chan:
                self.structured_ch2 = d
            elif 3 == chan:
                self.structured_ch3a = d
            elif 4 == chan:
                self.structured_ch3b = d
            elif 5 == chan:
                self.structured_ch4 = d
            elif 6 == chan:
                self.structured_ch5 = d
        elif 3 == ntype:
            if 1 == chan:
                self.measurement_ch1 = d
            elif 2 == chan:
                self.measurement_ch2 = d
            elif 3 == chan:
                self.measurement_ch3a = d
            elif 4 == chan:
                self.measurement_ch3b = d
            elif 5 == chan:
                self.measurement_ch4 = d
            elif 6 == chan:
                self.measurement_ch5 = d

    def file_length(self,filename):

        with open(filename) as fp:
            for i, l in enumerate(fp):
                pass
        return i + 1

    def read_file(self,avhrr_name,noise_type,channel):

        filename = '{0}_{1}_{2}.dat'.format(avhrr_name,noise_type,channel)
        dtype={'names': ('year', 'month', 'day', 'hour', 'minute', 'channel',\
                             'ndata','minval','maxval','x','x2','mean',\
                             'stdev','var','q_25','median','q_75','nan_fraction'),\
                   'formats': (np.int,np.int,np.int,np.int,np.int,np.int,\
                                   np.float,np.float,np.float,np.float,\
                                   np.float,np.float,np.float,np.float,\
                                   np.float,np.float,np.float,np.float)}
        d = np.loadtxt(filename,dtype=dtype)
        self.set_var(noise_type,channel,d)

    def __init__(self,avhrr_name):

        self.read_file(avhrr_name,'independent','Ch1')
        self.read_file(avhrr_name,'independent','Ch2')
        self.read_file(avhrr_name,'independent','Ch3a')
        self.read_file(avhrr_name,'independent','Ch3b')
        self.read_file(avhrr_name,'independent','Ch4')
        self.read_file(avhrr_name,'independent','Ch5')

        self.read_file(avhrr_name,'structured','Ch1')
        self.read_file(avhrr_name,'structured','Ch2')
        self.read_file(avhrr_name,'structured','Ch3a')
        self.read_file(avhrr_name,'structured','Ch3b')
        self.read_file(avhrr_name,'structured','Ch4')
        self.read_file(avhrr_name,'structured','Ch5')

        self.read_file(avhrr_name,'measurement','Ch1')
        self.read_file(avhrr_name,'measurement','Ch2')
        self.read_file(avhrr_name,'measurement','Ch3a')
        self.read_file(avhrr_name,'measurement','Ch3b')
        self.read_file(avhrr_name,'measurement','Ch4')
        self.read_file(avhrr_name,'measurement','Ch5')

def print_one(instr_name,chan_name,data):

    print '     '+chan_name
    gd = (data['minval'] > 0)
    if np.sum(gd) == 0:
        print 'No Data available'
        return
    iindex = np.arange(len(data['minval']),dtype=np.int32)    
    i_index = iindex[gd]
    ggd = (data['minval'][gd] == data['minval'][gd].min())
    indx = i_index[ggd]
    print ' min    : ',data['minval'][gd].min(),\
        data['year'][indx[0]],\
        data['month'][indx[0]],\
        data['day'][indx[0]],\
        data['hour'][indx[0]],\
        data['minute'][indx[0]]

    minval = data['minval'][gd].min()
    gd = (data['maxval'] > 0)
    iindex = np.arange(len(data['maxval']),dtype=np.int32)    
    i_index = iindex[gd]
    ggd = (data['maxval'][gd] == data['maxval'][gd].max())
    indx = i_index[ggd]
    print ' max    : ',data['maxval'][gd].max(),\
        data['year'][indx[0]],\
        data['month'][indx[0]],\
        data['day'][indx[0]],\
        data['hour'][indx[0]],\
        data['minute'][indx[0]]
    maxval = data['maxval'][gd].max()

# --------------------------------
# plot histogram of minimum values
#    min1 = maxval/(2^8-1)
#    min2 = maxval/(2^16-1)
#    min3 = maxval/(2^24-1)
#    min4 = maxval/(2^32-1)
#    gd = (data['minval'] > 0)
#    min_median = np.median(data['minval'][gd],axis=0)
#    plot_file = instr_name+" "+chan_name+".png"
#    plot_title = instr_name+" "+chan_name+": max = "+str("%02f" %maxval)+" median = "+str("%02f" %min_median)
#
#    fig, ax = plt.subplots()
#    umin = data['minval'][gd].min()
#    umax = data['minval'][gd].max()
#    ubin = np.arange(umin,umax,100)
#    udata = data['minval'][gd]
#    histogram = ax.hist(udata,bins='auto',range=None)
#    ax.tick_params(labelsize=12)
#    ax.set_xlabel("Uncertainty", fontsize=12)
#    ax.set_ylabel("Counts", fontsize=12)
#    ax.axvline(min1, color='k', linestyle='--')
#    ax.axvline(min2, color='r', linestyle='--')
#    ax.axvline(min3, color='b', linestyle='--')
#    ax.axvline(min4, color='m', linestyle='--')
#    plt.title(plot_title)
#    plt.savefig(plot_file)
# --------------------------------

    try:
        ratio = maxval/minval
        nbits = int(np.log(ratio)/np.log(2.))+1
        print ' max/min:', ratio,'   no. bits:',nbits
    except:
        print ' Could not get ratio/no of bits'

    gd = (data['mean'] > -1e20)
    if np.sum(gd) > 0:
        print ' mean   : ',data['mean'][gd].min(),data['mean'][gd].max()

    gd = (data['stdev'] > -1e20)
    if np.sum(gd) > 0:
        print ' stdev  : ',data['stdev'][gd].min(),data['stdev'][gd].max()
     
    gd = (data['median'] > -1e20)
    if np.sum(gd) > 0:
        print ' median : ',data['median'][gd].min(),data['median'][gd].max()

    gd = (data['q_25'] > -1e20) & (data['q_75'] > -1e20)
    if np.sum(gd) > 0:
        rsd = ((data['q_75'][gd] - data['q_25'][gd])/1.349)
        print ' RSD    : ',rsd.min(),rsd.max()
        print ' NaN fraction: ',data['nan_fraction'][indx[0]]

def print_data(avhrr_name):

    d = read_data(avhrr_name)
    
    print_one(avhrr_name,'Channel 1  (independent)',d.independent_ch1)
    print_one(avhrr_name,'Channel 2  (independent)',d.independent_ch2)
    print_one(avhrr_name,'Channel 3a (independent)',d.independent_ch3a)
    print_one(avhrr_name,'Channel 3b (independent)',d.independent_ch3b)
    print_one(avhrr_name,'Channel 4  (independent)',d.independent_ch4)
    print_one(avhrr_name,'Channel 5  (independent)',d.independent_ch5)

    print_one(avhrr_name,'Channel 1  (structured)',d.structured_ch1)
    print_one(avhrr_name,'Channel 2  (structured)',d.structured_ch2)
    print_one(avhrr_name,'Channel 3a (structured)',d.structured_ch3a)
    print_one(avhrr_name,'Channel 3b (structured)',d.structured_ch3b)
    print_one(avhrr_name,'Channel 4  (structured)',d.structured_ch4)
    print_one(avhrr_name,'Channel 5  (structured)',d.structured_ch5)

    print_one(avhrr_name,'Channel 1  (measurement)',d.measurement_ch1)
    print_one(avhrr_name,'Channel 2  (measurement)',d.measurement_ch2)
    print_one(avhrr_name,'Channel 3a (measurement)',d.measurement_ch3a)
    print_one(avhrr_name,'Channel 3b (measurement)',d.measurement_ch3b)
    print_one(avhrr_name,'Channel 4  (measurement)',d.measurement_ch4)
    print_one(avhrr_name,'Channel 5  (measurement)',d.measurement_ch5)

if __name__ == "__main__":
    print_data(sys.argv[1])

