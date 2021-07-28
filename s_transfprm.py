from typing import overload
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
import warnings 
from scipy import signal as sci_signal

def test():
    dt = 0.001
    chirp = signal.chirp_signal(dt)
    #power, freq = dft(chirp, 100, method='py').dftpy()
    #plt.plot(freq, power)
    #plt.show()

    specs = TimeFrequency(chirp, int(1/dt), 100, overlap=50).make_spectrogram()

    return

class signal:   
    '''Generate sample signals'''

    def chirp_signal(dt):
        """Generate a chirp signal for testing 
        Input
                        dt              the sampling interval
        Output          
                        x               chirp signal 
        """
        t = np.arange(0,3,dt)
        f0 = 50
        f1 = 250
        t1 = 2
        x = np.cos(2*np.pi*t*(f0 + (f1 - f0)*np.power(t, 2)/(3*t1**2)))
        fs = 1/dt

        return x

    def constant_signal(dt, freq):
        """Generate a steady signal for testing
        Input
                        freq                frequency
        Output
                        x                   steady signal 
        """
        x = np.cos(np.arange(0, 2*np.pi*freq, dt))

        return x

class dft:
    '''Distrete Fourier Transform class'''

    def __init__(self, data, sample_rate, method='numpy'):
        self.data = data
        self.sample_rate = sample_rate
        self.method = method
        self.length = len(data)

    def get_coeff(self, n):
        '''Get the Fourier coefficients of a DFT
        Input
                        data                the time-domain data
                        n                   the number of segments
        Output
                        coeff               the Fourier coefficient of a DFT
        '''
        
        ks = np.arange(0, self.length, 1)
        coeff = np.sum(self.data*np.exp((1j*2*np.pi*ks*n)/self.length))/self.length       

        return coeff
    
    def get_freq(self):
        '''Get the sampled frequency grid
        variables
                        length              length of the time-domain data          
                        sample_rate         sample_rate, the inverse of sample interval
        return 
                        freq                the frequency array
        '''
      
        ks = np.linspace(0, int(self.length/2), int(self.length/2))     # frequency grid
        freq = ks*self.sample_rate/self.length                          # scaled frequency grid

        return freq

    def dftpy(self):
        '''Get the DFT of input data
        variables
                        data                time-domain signal 
                        sample_rate         sample_rate, the inverse of sample interval
        Output
                        power               frequency-domain signal strength
                        freq                frequency array
        Note:
                        two methods avaliable - 'numpy' and 'py'
                                                'numpy' ==> np.fft.fft (faster)
                                                'py' ==> direct sum DFT (slower)
        '''
        length = int(len(self.data))                    # full frequency length 
        length_w_nyquest = int(length/2)                # Nyquest limit == 1/2 frequency length
        if self.method == 'numpy':                      # Numpy FFT
            power = np.fft.fft(self.data)[:length_w_nyquest]
            freq = np.fft.fftfreq(len(self.data), 1/self.sample_rate)[:length_w_nyquest]

        elif self.method == 'py':                       # direct sum DFT
            power = np.empty(length_w_nyquest)          # empty data frame
            for i in range(length_w_nyquest):
                power[i] = np.abs(self.get_coeff(i)*2, dtype=float)  # loop sum
            freq = self.get_freq()

        return power, freq

class TimeFrequency:
    '''S Transform and Spectrogram class'''

    def __init__(self, data, sample_rate, L, overlap = None, method = 'py'):
        '''
        Input
                        data                time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
                        L                   the data length included in each DFT
                        overlap             the length of overlap between DFTs
        '''
        self.data = data
        self.sample_rate = sample_rate
        self.L = L
        self.overlap = overlap
        self.method = method
        self.power, self.freq = dft(self.data, self.sample_rate, method=self.method).dftpy()

    def s_transform(self):
        '''Generate S Transform using dft/fft
        variables
                        data                time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
                        L                   the data length included in each DFT
                        overlap             the length of overlap between DFTs
        '''
        if self.overlap is None:
            self.overlap = self.L/2
        self.overlap = int(self.overlap)
        sections = np.arange(0, len(self.data), self.L-self.overlap, dtype=int)
        sections  = sections[sections + self.L < len(self.data)]
        xns = np.empty(len(sections))       # fix array
        xns = []
        count = 0
        for begin in sections:
            # short term discrete fourier transform
            ts_window, freq = dft(self.data[begin:begin + self.L], self.sample_rate, method=self.method).dftpy() 
            #xns[count] = ts_window         # fix array
            xns.append(ts_window)
            count+=1
        specX = np.array(xns).T
        spec = 10*np.log10(specX)
        assert spec.shape[1] == len(sections) 
        
        return sections, spec

    def make_spectrogram(self):
        '''Generate the spectrogram
        variables
                        data                time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
        return 
                        plt_spec            image specs
        '''
        data_interval = len(self.data)/self.sample_rate
        
        self.sections, self.spec = self.s_transform()
        plt.figure()
        plt_spec = plt.imshow(self.spec, origin='lower')

        # y axis
        Nyticks = 10
        freq = np.linspace(0, self.spec.shape[0], Nyticks, dtype=int)
    
        downsampled_freq_Hz =  np.empty(Nyticks)
        counter = 0
        for i in range(len(self.freq)):
            if i%int(len(self.freq)/Nyticks) == 0:
                downsampled_freq_Hz[counter] = "{:4.2f}".format(self.freq[i])
                counter+=1
        plt.yticks(freq, downsampled_freq_Hz)
        plt.ylabel("Frequency (Hz)")

        # x axis 
        Nxticks = 10
        ts_spec = np.linspace(0,self.spec.shape[1],Nxticks)
        ts_spec_sec  = ["{:4.2f}".format(i) for i in np.linspace(0,data_interval*self.sections[-1]/len(self.data),Nxticks)]
        plt.xticks(ts_spec,ts_spec_sec)
        plt.xlabel("Time (sec)")

        plt.title("Spectrogram L={} Spectrogram.shape={}".format(self.L,self.spec.shape))
        plt.colorbar(use_gridspec=True)
        plt.show()

        return plt_spec

if __name__ == '__main__':
    test()