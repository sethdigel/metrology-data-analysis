import numpy as np
import scipy.optimize

class XyzPlane(object):
    def __init__(self, a, b, c):
        self.pars = a, b, c
    def __call__(self, positions):
        a, b, c = self.pars
        return np.array([a*x + b*y + c for x, y in positions])

def xyz_plane(positions, a, b, c):
    my_plane = XyzPlane(a, b, c)
    return my_plane(positions)

class ContourData(dict):
    def __init__(self, infile):
        super(ContourData, self).__init__()
        self.infile = infile
        key = None
        for line in open(infile):
            if line.startswith('Contour'):
                key = line.strip()
                self[key] = []
            if line.strip() and key is not None:
                try:
                    coords = self._xyz(line)
                    self[key].append(coords)
                except ValueError:
                    pass
            else:
                key = None
        for key in self:
            self[key] = zip(*tuple(self[key]))
        self.npts_sensor = max([len(x[0]) for x in self.values()])
    def _xyz(self, line):
        # Unpack a line and convert z values from mm to microns.
        data = [float(x) for x in line.split()[:3]]
        return data[0], data[1], 1e3*data[2]
    def all_points(self):
        xx, yy, zz = [], [], []
        for arr in self.values():
            xx.extend(arr[0])
            yy.extend(arr[1])
            zz.extend(arr[2])
        return xx, yy, zz
    def ref_data(self):
        xx, yy, zz = [], [], []
        for arr in self.values():
            if len(arr[0]) < self.npts_sensor:
                xx.extend(arr[0])
                yy.extend(arr[1])
                zz.extend(arr[2])
        return np.array(zip(xx, yy)), np.array(zz)
    def sensor_data(self):
        for arr in self.values():
            if len(arr[0]) == self.npts_sensor:
                return np.array(zip(arr[0], arr[1])), np.array(arr[2])
    def sensorPlane_fit(self, nsigma=3, p0=(0, 0, 0)):
        return self.xyzPlane_fit(*self.sensor_data(), nsigma=nsigma, p0=p0)
    def refPlane_fit(self, nsigma=3, p0=(0, 0, 0)):
        return self.xyzPlane_fit(*self.ref_data(), nsigma=nsigma, p0=p0)
    def xyzPlane_fit(self, positions, zvals, nsigma=3, p0=(0, 0, 0)):
        """
        Fit xyz data, clipping the initial fit at the nsigma level to remove
        outlier points.
        """
        # Initial fit
        pars, _ = scipy.optimize.curve_fit(xyz_plane, positions, zvals, p0=p0)
        dz = xyz_plane(positions, *pars) - zvals
        mean, stdev = np.mean(dz), np.std(dz)

        # Refit the reference data within nsigma*stdev of the mean.
        index = np.where((dz > mean-nsigma*stdev) & (dz < mean+nsigma*stdev))
        pars, _ = scipy.optimize.curve_fit(xyz_plane, positions[index],
                                           zvals[index], p0=pars)

        # Return a XyzPlane functor initialized with the fitted parameters.
        return XyzPlane(*pars)
    def plot_statistics(self, xyzPlane, plane_data=None, nsigma=4, title=None):
        if plane_data is None:
            pos, z = data.sensor_data()
        else:
            pos, z = plane_data
        dz = z - xyzPlane(pos)
        mean, stdev = np.mean(dz), np.std(dz)

        # Trim outliers at nsigma and recompute mean and stdev.
        index = np.where((dz > mean-nsigma*stdev) & (dz < mean+nsigma*stdev))
        mean, stdev = np.mean(dz[index]), np.std(dz[index])

        win0 = plot.histogram(dz[index],
                              xname=r'$z - z_{\rm model}$ $(\mu{\rm m})$',
                              yname='entries/bin')
        plot.pylab.annotate('mean=%.3f\nstdev=%.3f\n%i-sigma clip' 
                            % (mean, stdev, nsigma),
                            (0.05, 0.8), xycoords='axes fraction')
        if title is None:
            title = self.infile
        win0.set_title(title)

if __name__ == '__main__':
    import pylab_plotter as plot
    plot.pylab.ion()

    infile = '112-05_150811a.DAT'
#    infile = '113-11_Abs_Hgt_01.DAT'
    data = ContourData(infile)

    xx, yy, zz = data.all_points()
    win = plot.xyplot(xx, yy, xrange=(-40, 80), yrange=(-40, 80),
                      xname='x (mm)', yname='y (mm)')
    win.set_title(infile)
    plot.save('xy_points_%s.png' % infile.split('.DAT')[0])

    #
    # Reference plane fit statistics
    #
    refPlane = data.refPlane_fit()
    data.plot_statistics(refPlane, plane_data=data.ref_data(),
                         title='Reference plane residuals, %s' % infile)
    plot.save('refplane_stats_%s.png' % infile.split('.DAT')[0])

    #
    # Absolute height statistics
    #
    data.plot_statistics(refPlane, title='Sensor Absolute Height, %s' % infile)
    plot.save('abs_height_stats_%s.png' % infile.split('.DAT')[0])

    #
    # Flatness Statistics
    #
    sensorPlane = data.sensorPlane_fit()
    data.plot_statistics(sensorPlane, title='Sensor Flatness, %s' % infile)
    plot.save('flatness_stats_%s.png' % infile.split('.DAT')[0])
