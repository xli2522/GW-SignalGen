# Attempt to create a customized S and Q Transform; Q Transform reference == gwpy 2.0.4
# Reference:
#https://gwpy.github.io/docs/stable/api/gwpy.timeseries.TimeSeries.html#gwpy.timeseries.TimeSeries.q_transform
#https://github.com/gwpy/gwpy/blob/v2.0.4/gwpy/timeseries/timeseries.py#L2120
#https://github.com/gwpy/gwpy/blob/26f63684db17104c5d552c30cdf01248b2ec76c9/gwpy/signal/qtransform.py#L267

import numpy as np
from gwpy.timeseries import TimeSeries
from scipy.interpolate import (interp2d, InterpolatedUnivariateSpline)

def __init__(self, plane, energies, search):
        self.plane = plane
        self.energies = energies
        self.peak = self._find_peak(search)

def spectrogram(self):
        
    frequencies = self.plane.frequencies
    dtype = self.energies[0].dtype
    # build regular Spectrogram from peak-Q data by interpolating each
    # (Q, frequency) `TimeSeries` to have the same time resolution
    
    tres = abs(Segment(outseg)) / 1000.
    xout = np.arange(*outseg, step=tres)
    nx = xout.size
    ny = frequencies.size
    out = Spectrogram(np.empty((nx, ny), dtype=dtype),
                      t0=outseg[0], dt=tres, frequencies=frequencies)
    # record Q in output
    out.q = self.plane.q
    # interpolate rows
    for i, row in enumerate(self.energies):
        xrow = np.arange(row.x0.value, (row.x0 + row.duration).value,
                            row.dx.value)
        interp = InterpolatedUnivariateSpline(xrow, row.value)
        out[:, i] = interp(xout).astype(dtype, casting="same_kind",
                                        copy=False)
    if fres is None:
        return out
    # interpolate the spectrogram to increase its frequency resolution
    # --- this is done because Duncan doesn't like interpolated images
    #     since they don't support log scaling
    interp = interp2d(xout, frequencies, out.value.T, kind='cubic')
    if not logf:
        if fres == "<default>":
            fres = .5
        outfreq = np.arange(
            self.plane.frange[0], self.plane.frange[1], fres,
            dtype=dtype)
    else:
        if fres == "<default>":
            fres = 500
        outfreq = np.geomspace(
            self.plane.frange[0],
            self.plane.frange[1],
            num=int(fres),
        )
    new = type(out)(
        interp(xout, outfreq).T.astype(
            dtype, casting="same_kind", copy=False),
        t0=outseg[0], dt=tres, frequencies=outfreq,
    )
    new.q = self.plane.q
    return new

def q_spectrogram(q_range=[80,100], frange=[30, 500]):

    return 