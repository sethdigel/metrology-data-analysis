#!/usr/bin/env python
import glob
import lcatr.schema
import siteUtils

sensor_id = siteUtils.getUnitId()

png_files = glob.glob('%(sensor_id)s_abs_height*.png'% locals())
results = [lcatr.schema.fileref.make(x) for x in png_files]

txt_files = glob.glob('%(sensor_id)s_abs_height*.txt'% locals())
results.extend([lcatr.schema.fileref.make(x) for x in txt_files])

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
