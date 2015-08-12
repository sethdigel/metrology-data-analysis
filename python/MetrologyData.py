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

class PointCloud(object):
    def __init__(self, x, y, z):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
    def data(self):
        """
        Return the xyz data repackaged in a format appropriate for
        fitting using scipy.optimize.curve_fit and the xyz_plane function.
        """
        return np.array(zip(self.x, self.y)), self.z
    def __add__(self, other):
        result = PointCloud([], [], [])
        result.x = np.concatenate((self.x, other.x))
        result.y = np.concatenate((self.y, other.y))
        result.z = np.concatenate((self.z, other.z))
        return result
    def xyzPlane_fit(self, nsigma=3, p0=(0, 0, 0)):
        """
        Fit a plane to the xyz data, clipping the initial fit at the 
        nsigma level to remove outlier points.
        """
        positions = np.array(zip(self.x, self.y))
        
        # Initial fit
        pars, _ = scipy.optimize.curve_fit(xyz_plane, positions, self.z, p0=p0)
        dz = xyz_plane(positions, *pars) - self.z
        mean, stdev = np.mean(dz), np.std(dz)

        # Refit the reference data within nsigma*stdev of the mean.
        index = np.where((dz > mean-nsigma*stdev) & (dz < mean+nsigma*stdev))
        pars, _ = scipy.optimize.curve_fit(xyz_plane, positions[index],
                                           self.z[index], p0=pars)

        # Return a XyzPlane functor initialized with the fitted parameters.
        return XyzPlane(*pars)

class MetrologyData(object):
    def __init__(self, infile):
        self.infile = infile
    def sensor_data(self):
        raise NotImplementedError('Subclass has not implemented this function.')
    def plot_statistics(self, plane_functor, plane_data=None, 
                        nsigma=4, title=None, zoffset=0):
        if plane_data is None:
            pos, z = self.sensor_data()
        else:
            pos, z = plane_data
        dz = z - plane_functor(pos) + zoffset
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
        return win0, dz

class OgpData(MetrologyData):
    def __init__(self, infile):
        super(OgpData, self).__init__(infile)
        self.infile = infile
        self._read_data()
    def _read_data(self):
        # Read the "Contour" data blocks into a local dict, convert
        # each to a PointCloud object, sort by mean y-value, and
        # finally, set the sensor and reference datasets.
        data = dict()
        key = None
        for line in open(self.infile):
            if line.startswith('Contour'):
                key = line.strip()
                data[key] = []
            if line.strip() and key is not None:
                try:
                    coords = self._xyz(line)
                    data[key].append(coords)
                except ValueError:
                    pass
            else:
                key = None
        # Convert to PointCloud objects.
        for key in data:
            data[key] = PointCloud(*zip(*tuple(data[key])))
        if len(data) != 3:
            raise RuntimeError("Expected exactly 3 Contour data blocks in the OGP data. %i found." % len(data))
        # Identify sensor and reference point clouds by mean y-values.
        # The sensor dataset is in the middle.
        yavgs = sorted([np.mean(cloud.y) for cloud in data.values()])
        ref_clouds = []
        for cloud in data.values():
            if np.mean(cloud.y) == yavgs[1]:
                self.sensor = cloud
            else:
                ref_clouds.append(cloud)
        self.reference = ref_clouds[0] + ref_clouds[1]
    def _xyz(self, line):
        # Unpack a line and convert z values from mm to microns.
        data = [float(x) for x in line.split()[:3]]
        return data[0], data[1], 1e3*data[2]
    def sensor_data(self):
        return self.sensor.data()
    def sensorPlane_fit(self, nsigma=3, p0=(0, 0, 0)):
        return self.sensor.xyzPlane_fit(nsigma=nsigma, p0=p0)
    def refPlane_fit(self, nsigma=3, p0=(0, 0, 0)):
        return self.reference.xyzPlane_fit(nsigma=nsigma, p0=p0)

if __name__ == '__main__':
    import pylab_plotter as plot
    plot.pylab.ion()

    infile = '112-05_150811a.DAT'
    ogpData = OgpData(infile)

    #
    # Reference plane fit statistics
    #
    refPlane = ogpData.refPlane_fit()
    ogpData.plot_statistics(refPlane, plane_data=ogpData.reference.data(),
                         title='Reference plane residuals, %s' % infile)
    plot.save('refplane_stats_%s.png' % infile.split('.DAT')[0])

    #
    # Absolute height statistics
    #
    ogpData.plot_statistics(refPlane,
                            title='Sensor Absolute Height, %s' % infile,
                            zoffset=-1.2)
    plot.save('abs_height_stats_%s.png' % infile.split('.DAT')[0])

    #
    # Flatness Statistics
    #
    sensorPlane = ogpData.sensorPlane_fit()
    ogpData.plot_statistics(sensorPlane, title='Sensor Flatness, %s' % infile)
    plot.save('flatness_stats_%s.png' % infile.split('.DAT')[0])
