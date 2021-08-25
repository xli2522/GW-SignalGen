import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
import warnings 
import time
from scipy import signal as sci_signal

def test():
    dt = 0.001
    chirp = signal.chirp_signal(dt)
    #power, freq = dft(chirp, 100, fftmethod='py').dftpy()
    #plt.plot(freq, power)
    #plt.show()
    beginning = time.time()
    #specs1 = STimeFrequency(chirp, int(1/dt), 100, method='dtft').make_spectrogram()
    #dtftTime = time.time()
    specs2 = STimeFrequency(chirp, int(1/dt), 200, fftmethod='numpy', method='dst').make_spectrogram()
    dstTime = time.time()
    #print('DTFT: ' + str(dtftTime - beginning))
    #print('DST: ' + str(dstTime - dtftTime))
    
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

class dftMethods:
    '''Distrete Fourier Transform class'''

    def __init__(self, data, sample_rate, fftmethod='numpy'):
        self.data = data
        self.sample_rate = sample_rate
        self.fftmethod = fftmethod
        self.length = len(data)

    def get_coeff(self, n):
        '''Get the Fourier coefficients of a DFT
        Input
                        data                the time-domain data
                        n                   the number of segments
        Output
                        dftCoeff               the Fourier coefficient of a DFT
        '''
        ks = np.arange(0, self.length, 1)
        dftCoeff = np.sum(self.data*np.exp((1j*2*np.pi*ks*n)/self.length))/self.length       

        return dftCoeff

    def get_S_coeff(self, n):
        '''Get the Fourier coefficients of a S Transform
        Input
                        data                the time-domain data
                        n                   the number of segments
        Output
                        coeff               the Fourier coefficient of a DFT
        '''
        ks = np.arange(1, self.length+1, 1)             # shifted 1 to avoid 0 division
        ms = np.arange(1, self.length+1, 1)
        js = np.arange(1, self.length+1, 1)

        
        dftsCoeff = 1/self.length*np.sum(self.data*np.exp((1j*2*np.pi*(ks+ms)*n)/self.length))      
        coeff = np.sum(dftsCoeff*np.exp(-2*np.pi**2*ms**2/ks**2)*np.exp(1j*2*np.pi*ms*js/self.length))
       
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

    def dstpy(self):
        '''Get the DST of imput data
        variables
                        data                time-domain signal 
                        sample_rate         sample_rate, the inverse of sample interval
        Output
                        power               frequency-domain signal strength
                        freq                frequency array
                        time                time array
        Associated functions
                        get_freq
                        get_time
        '''
        length = int(len(self.data))                    # full frequency length 
        length_w_nyquest = int(length/2)                # Nyquest limit == 1/2 frequency length
        half_length = int(1/2*length_w_nyquest)         # fold frequency length one more time
        power = np.empty(half_length)              # empty data frame

        if self.fftmethod == 'py':
            for i in range(half_length):
                power[i] = np.abs(self.get_S_coeff(i)*2, dtype=float)  # loop sum
            freq = self.get_freq()
        elif self.fftmethod == 'numpy':
            ks = np.arange(1, self.length+1, 1)             # shifted 1 to avoid 0 division
            ms = np.arange(1, self.length+1, 1)
            js = np.arange(1, self.length+1, 1)
            dfts = np.fft.fft(self.data)
            dsts = np.fft.fft(dfts*np.exp(-2*np.pi**2*ms**2/ks**2)*np.exp(1j*2*np.pi*ms*js/self.length))
            power = dsts
            freq = self.get_freq()
    
        return power, freq

    def dftpy(self):
        '''Get the DFT of input data
        variables
                        data                time-domain signal 
                        sample_rate         sample_rate, the inverse of sample interval
        Output
                        power               frequency-domain signal strength
                        freq                frequency array
        Associated functions
                        get_freq
        Note:
                        two fftmethods avaliable - 'numpy' and 'py'
                                                'numpy' ==> np.fft.fft (faster)
                                                'py' ==> direct sum DFT (slower)
        '''
        length = int(len(self.data))                    # full frequency length 
        length_w_nyquest = int(length/2)                # Nyquest limit == 1/2 frequency length
        if self.fftmethod == 'numpy':                      # Numpy FFT
            power = np.fft.fft(self.data)[:length_w_nyquest]
            freq = np.fft.fftfreq(len(self.data), 1/self.sample_rate)[:length_w_nyquest]

        elif self.fftmethod == 'py':                       # direct sum DFT
            power = np.empty(length_w_nyquest)          # empty data frame
            for i in range(length_w_nyquest):
                power[i] = np.abs(self.get_coeff(i)*2, dtype=float)  # loop sum
            freq = self.get_freq()

        return power, freq

class STimeFrequency:
    '''Transform and Spectrogram class'''

    def __init__(self, data, sample_rate, L, overlap = None, fftmethod = 'py', method = 'dtft'):
        '''
        Input
                        data                time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
                        L                   the data length included in each DFT
                        overlap             the length of overlap between DFTs
                        fftmethod           fft/dft implimentation
                        method              spectrogram method - dtft/dst
        '''
        self.data = data
        self.sample_rate = sample_rate
        self.L = L
        self.overlap = overlap
        self.fftmethod = fftmethod
        self.method = method
        
    def transform(self):
        '''Generate S Transform sections using dft/fft
        variables
                        data                time-domain data
                        sample_rate         sample_rate, the inverse of sample interval
                        L                   the data length included in each DFT
                        overlap             the length of overlap between DFTs
        return 
                        sections
                        spec
        '''
        if self.overlap is None:
            self.overlap = int(self.L/2)
        sections = np.arange(0, len(self.data), self.L-self.overlap, dtype=int)
        sections  = sections[sections + self.L < len(self.data)]
        xns = np.empty(len(sections))       # fix array
        xns = []
        count = 0
        if self.method == 'dtft':
            for begin in sections:
                # short term discrete fourier transform
                ts_window, self.freq = dftMethods(self.data[begin:begin + self.L], self.sample_rate, fftmethod=self.fftmethod).dftpy() 
                #xns[count] = ts_window         # fix array
                xns.append(ts_window)
                count+=1
        elif self.method == 'dst':
            for begin in sections:
                # short term discrete fourier transform
                ts_window, self.freq = dftMethods(self.data[begin:begin + self.L], self.sample_rate, fftmethod=self.fftmethod).dstpy() 
                #xns[count] = ts_window         # fix array
                xns.append(ts_window)
                count+=1
        else:
            warnings.warn('unsupported method')

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
        self.sections, self.spec = self.transform()

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
        #plt.savefig(str(self.method)+'_spectrogram.png')            #+str(time.time())
        #plt.close()

        return plt_spec

if __name__ == '__main__':
    test()