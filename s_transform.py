import time
#import random
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as scisignal
#from scipy import fft, ifft
#from scipy.interpolate import interp1d

def test():
    dt = 0.001
    chirp = signal.chirp_signal(dt)
    #power, freq = dft(chirp, 100, fftmethod='py').dftpy()
    #plt.plot(freq, power)
    #plt.show()
    beginning = time.time()
    specs1 = TimeFrequency(chirp, int(1/dt), method='stft', show=False, savefig=True).plot()
    stftTime = time.time()
    specs2 = TimeFrequency(chirp, sample_rate=int(1/dt), method='dst', show=False, savefig=True).plot()
    dstTime = time.time()
    print('STFT: ' + str(stftTime - beginning))
    print('DST: ' + str(dstTime - stftTime))
        
    return 

class signal: 
    '''Generate sample signals'''

    def chirp_signal(dt):
        '''Generate a chirp signal for testing 
        Input
                        dt              the sampling interval
        Output          
                        x               chirp signal 
        '''
        t = np.arange(0,3,dt)
        f0 = 50
        f1 = 250
        t1 = 2
        x = np.cos(2*np.pi*t*(f0 + (f1 - f0)*np.power(t, 2)/(3*t1**2)))
        fs = 1/dt

        return x

class dftMethods:
    '''Distrete Fourier Transform class
    NOTE: could be removed and combined with class TimeFrequency
    '''

    def __init__(self, data, sample_rate, fftmethod='numpy'):
        self.data = data
        self.sample_rate = sample_rate
        self.fftmethod = fftmethod
        self.length = len(data)

class TimeFrequency:
    '''Time-Frequency Analysis Methods class'''

    def __init__(self, ts, sample_rate=4096, frange=None, frate=1, overlap = None, method = 'dst', show=True, savefig=False):
        '''
        Input
                        ts                  time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
                        L                   the data length included in each DFT
                        frange              the frequency range
                        frate               the frequency sample rate
                        overlap             the length of overlap between DFTs
                        fftmethod           fft/dft implimentation
                        method              spectrogram method - stft/dst/dcst
        '''
        self.ts = ts
        self.sample_rate = sample_rate
        self.frange = frange
        self.frate = frate
        self.overlap = overlap
        self.method = method
        self.show = show
        self.savefig = savefig
        self.length = len(self.ts)

        if self.frange == None:
            # use default frange
            self.frange = [30, 500]
    
    def _window(self, N, nleft=0, nright=0):
        ''' *Planck* Window function
        Imput
                        N               the number of data points included in a window
                        nleft           the number of data points to taper on the left
                        nright          the number of data points to taper on the right
        Return
                        win             the normalized window
        Note: Window scaling factor = 1, normalized
        '''
        from scipy.special import expit

        win = np.ones(N)
        if nleft:
            win[0] *= 0
            zleft = np.array([nleft * (1./k + 1./(k-nleft))
                                for k in range(1, nleft)])
            win[1:nleft] *= expit(-zleft)
        if nright:
            win[N-1] *= 0
            zright = np.array([-nright * (1./(k-nright) + 1./k)
                                for k in range(1, nright)])
            win[N-nright:N-1] *= expit(-zright)

        return win

    def _window_normal(self, freq):
        '''Gaussian Window function
        Input 
                        freq            the number of min frequency bins
        Return
                        win             the normalized gaussian window
        Note: Window scaling factor = 1, normalized
        '''
        gauss = scisignal.gaussian(self.length,std=(freq)/(2*np.pi))
        win = np.hstack((gauss,gauss))[self.length//2:self.length//2+self.length]

        return win

    def stransform(self):
        '''The Stockwell Transform'''
       
        Nfreq = [int(self.frange[0]*self.length/self.sample_rate), int(self.frange[1]*self.length/self.sample_rate)]               # the number of data points for min and max frequencies
        tsVal = np.copy(self.ts)            # copy ts values
        power = np.zeros((int((Nfreq[1]-Nfreq[0])/self.frate)+1,self.length), dtype='c8')            #complex64 C
        tsFFT = np.fft.fft(tsVal)
        vec = np.hstack((tsFFT, tsFFT))
        
        if self.frange[0] == 0:
            power = np.mean(tsVal)*np.ones(self.length)
        else:
            power[0] = np.fft.ifft(vec[Nfreq[0]:Nfreq[0]+self.length]*self._window_normal(Nfreq[0]))
        for i in range(self.frate, (Nfreq[1]-Nfreq[0])+1, self.frate):
            power[int(i/self.frate)] = np.fft.ifft(vec[Nfreq[0]+i:Nfreq[0]+i+self.length]*self._window_normal(Nfreq[0]+i))

        return np.abs(power)

    def get_freq_Hz(self, ks):
        '''Get the frequency label in Hz
        Input
                        ks              the freqeuncy array
        Output
                        ksHz            the frequency array in Hz
        ''' 
        ksHz = ks*self.sample_rate/self.length
        ksHz = [int(i) for i in ksHz]              #to be simplified
          
        return ksHz

    def plot(self):
        '''Generate the spectrogram
        methods: stft / dst / dcst
        '''
        if self.method == 'stft':
            f, t, Sxx = scisignal.spectrogram(
                                self.ts, 
                                self.sample_rate,
                                nperseg=128,
                                noverlap=64, 
                                nfft=5000, 
                                scaling='spectrum'
                            )
           
            plt.figure(figsize=(12,9))
            plt.pcolor(t,f, Sxx, cmap ='jet')
            plt.ylabel('Frequency [Hz]')
            plt.xlabel('Time [sec]')
            #plt.xlim(7.5,15)
            plt.ylim(0,250)
            plt.colorbar()

        elif self.method == 'dst':

            # s table
            sTable = self.stransform()
            # y axis
            self.fscale = int(self.length/self.sample_rate)
            y_sticksN = 10
            ks = np.linspace(self.frange[0], self.frange[1]*self.fscale, y_sticksN)
            ksHz = self.get_freq_Hz(ks)

            # x axis
            x_sticksN = 10
            ts = np.linspace(0, sTable.shape[1], x_sticksN)
            tsSec = ["{:4.2f}".format(i) for i in np.linspace(0, sTable.shape[1]/self.sample_rate, x_sticksN)]

            extent=(0,sTable.shape[1], self.fscale*self.frange[0], self.fscale*self.frange[1])
    
            plt.figure(figsize=(12,9))
            plt.imshow(sTable, origin='lower', extent=extent)
            plt.xticks(ts,tsSec)
            plt.yticks(ks,ksHz) 
            plt.xlabel("Time (sec)")
            plt.ylabel("Freq (Hz)")
            plt.colorbar()
        
        elif self.method == 'dcst':
            pass
        else:
            warnings.warn('Time-Freqeuncy method not supported.')

        if self.show:
            plt.show()
        
        if self.savefig:
            plt.savefig(str(self.method)+str(time.time())+'.png')
        
        return

if __name__ == '__main__':
    test()
