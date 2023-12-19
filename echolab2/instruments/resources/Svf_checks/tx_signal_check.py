
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert


f0 = 92000.0
f1 = 158000.0
tau = 0.002047999994829297
fs = 1500000.0
slope = 0.01061480026692152
txpower = 1000.0
z = 75.0



def ek80_chirp(txpower, fstart, fstop, slope, tau, z, rx_freq):
    '''
    '''

    # Create transmit signal
    sf = 1.0 / rx_freq
    a  = np.sqrt((txpower / 4.0) * (2.0 * z))
    t = np.arange(0, tau, sf)
    nt = t.shape[0]
    nwtx = int(2 * np.floor(slope * nt))
    wtx_tmp = np.hanning(nwtx)
    nwtxh = int(np.round(nwtx / 2))
    wtx = np.concatenate([wtx_tmp[0:nwtxh], np.ones((nt - nwtx)), wtx_tmp[nwtxh:]])
    beta = (fstop - fstart) * (tau ** -1)
    chirp = np.cos(2.0 * np.pi * (beta / 2.0 * (t ** 2) + fstart * t))
    y_tmp = a * chirp * wtx

    # The transmit signal must have a max amplitude of 1
    y = y_tmp / np.max(np.abs(y_tmp))

    return t, y


def chirp(t, f0, t1, f1):
    a = np.pi * (f1 - f0) / t1
    b = 2 * np.pi * f0
    return np.cos(a * t * t + b * t)


def hann(L):
    n = np.arange(0, L, 1)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (L - 1)))


def generateIdealWindowedTransmitSignal(f0, f1, tau, fs, slope):
        nsamples = int(np.floor(tau * fs))
        t = np.linspace(0, nsamples - 1, num=nsamples) * 1 / fs
        y = chirp(t, f0, tau, f1)
        L = int(np.round(tau * fs * slope * 2.0))  # Length of hanning window
        w = hann(L)
        N = len(y)
        w1 = w[0:int(len(w) / 2)]
        w2 = w[int(len(w) / 2):-1]
        i0 = 0
        i1 = len(w1)
        i2 = N - len(w2)
        i3 = N
        y[i0:i1] = y[i0:i1] * w1
        y[i2:i3] = y[i2:i3] * w2
        return y, t



def plotytx(f_0, f_1, tau, f_s, y_tx_n, slope, title):

    # Example of ideal windowed transmit signal with slope 0.5
    y_tx_n05slope, t = generateIdealWindowedTransmitSignal(
        f_0, f_1, tau, f_s, .5)

    plt.figure()
    if title == 'Echolab':
        plt.plot(t * 1000, np.abs(hilbert(y_tx_n[0:-1])), t * 1000, np.abs(hilbert(y_tx_n05slope)))
    else:
        plt.plot(t * 1000, np.abs(hilbert(y_tx_n)), t * 1000, np.abs(hilbert(y_tx_n05slope)))
    plt.title(title)
    #    'Ideal windowed transmit pulse.{:.0f}kHz - {:.0f}kHz, slope {:.3f}'
    #        .format(f_0 / 1000, f_1 / 1000, slope))
    plt.xlabel('Time (ms)')
    plt.ylabel('Envelope ()')
    #plt.savefig('./Paper/Fig_ytx.png',dpi=300)




t, y_tx_n_el = ek80_chirp(txpower, f0, f1, slope, tau, z, fs)
plotytx(f0, f1, tau, fs, y_tx_n_el, slope, 'Echolab')


y_tx_n_no, t =  generateIdealWindowedTransmitSignal(f0, f1, tau, fs, slope)
plotytx(f0, f1, tau, fs, y_tx_n_no, slope, '.no')


plt.figure()
plt.plot(t * 1000, (y_tx_n_el[0:-1])-y_tx_n_no)

plt.show()
