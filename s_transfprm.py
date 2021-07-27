import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
import warnings 

def test():
    chirp = signal.chirp_signal()
    power, freq = dft(chirp, 100, method='py').dftpy()
    plt.plot(freq, power)
    plt.show()

    return

class signal:   

    def chirp_signal():
        """Generate a chirp signal for testing 
        Output          
                        x               chirp signal 
        """
        dt = 0.01
        t = np.arange(0,3,dt)
        f0 = 50
        f1 = 250
        t1 = 2
        x = np.cos(2*np.pi*t*(f0 + (f1 - f0)*np.power(t, 2)/(3*t1**2)))
        fs = 1/dt

        return x

    def constant_signal(freq):
        """Generate a steady signal for testing
        Input
                        freq                frequency
        Output
                        x                   steady signal 
        """
        dt = 0.01
        x = np.cos(np.arange(0, 2*np.pi*freq, dt))

        return x

class dft:
    '''Distrete Fourier Transform class'''

    def __init__(self, data, sample_rate, method='numpy'):
        self.data = data
        self.sample_rate = sample_rate
        self.method = method

    def get_coeff(self, n):
        '''Get the Fourier coefficients of a DFT
        Input
                        data                the time-domain data
                        n                   the number of segments
        Output
                        coeff               the Fourier coefficient of a DFT
        '''
        length = len(self.data)
        ks = np.arange(0, length, 1)
        coeff = np.sum(self.data*np.exp((1j*2*np.pi*ks*n)/length))/length

        return coeff
    
    def get_freq(self):
        length = len(self.data)
        ks = np.linspace(0, int(length/2), int(length/2)) 
        freq = ks*self.sample_rate/length

        return freq

    def dftpy(self):
        '''Get the DFT of input data
        Input
                        data                time-domain signal 
        Output
                        power               frequency-domain signal strength
                        freq                frequency array
        '''
        length = int(len(self.data))
        length_w_nyquest = int(length/2)
        if self.method == 'numpy':
            power = np.fft.fft(self.data)[:length_w_nyquest]
            freq = np.fft.fftfreq(len(self.data), 1/self.sample_rate)[:length_w_nyquest]

        elif self.method == 'py':
            power = np.empty(length_w_nyquest)
            for i in range(length_w_nyquest):
                power[i] = np.abs(self.get_coeff(i)*2)
            freq = self.get_freq()

        return power, freq

if __name__ == '__main__':
    test()