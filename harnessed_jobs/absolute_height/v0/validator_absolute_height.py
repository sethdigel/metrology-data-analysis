#!/usr/bin/env python
import os
import glob
import numpy as np
import lcatr.schema
import siteUtils
from MetrologyData import md_factory

sensor_id = siteUtils.getUnitId()
ccd_vendor = siteUtils.getCcdVendor()

png_files = glob.glob('%(sensor_id)s_abs_height*.png'% locals())
results = [lcatr.schema.fileref.make(x) for x in png_files]

txt_files = glob.glob('%(sensor_id)s_abs_height*.txt'% locals())
results.extend([lcatr.schema.fileref.make(x) for x in txt_files])

#
# Extract numerical results from pickled MetrologyData object, if it exists.
#
pickle_file = 'abs_height.pickle'
if os.path.isfile(pickle_file):
    ZNOM = dict(ITL=12998., e2v=13000.)
    sensorData = md_factory.load(pickle_file)
    dzdx, dzdy, z0 = sensorData.plane_functor.pars
    zmean = np.mean(sensorData.sensor.z)
    znom = ZNOM[ccd_vendor]
    znom_residual_025 = sensorData.quantiles['0.025']
    znom_residual_975 = sensorData.quantiles['0.975']
    results.append(lcatr.schema.valid(lcatr.schema.get('absolute_height'),
                                      dzdx=dzdx, dzdy=dzdy, z0=z0,
                                      zmean=zmean, znom=znom, 
                                      znom_residual_025=znom_residual_025,
                                      znom_residual_975=znom_residual_975))

results.append(siteUtils.packageVersions())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
