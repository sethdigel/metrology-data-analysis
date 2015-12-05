import os
import fnmatch
from DataCatalog import DataCatalog
import siteUtils

def _folder(sensor_id, root_folder='LSST/vendorData'):
    ccd_manu = sensor_id.split('-')[0]
    my_folder = os.path.join(root_folder, ccd_manu, sensor_id)
    return my_folder

def get_met_scan_data(sensor_id, pattern, root_folder='LSST/vendorData',
                      site='slac.lca.archive', sort=False,
                      description='Metrology Scan Files:'):
    dc = DataCatalog(folder=_folder(sensor_id), site=site)

    query = '&&'.join(('DATA_PRODUCT=="MET_SCAN"',
                       'TEST_CATEGORY=="MET"',
                       'LSST_NUM=="%(sensor_id)s"'))
    query = query % locals()

    datasets = dc.find_datasets(query)
    file_list = []
    for item in datasets.full_paths():
        if fnmatch.fnmatch(os.path.basename(item), pattern):
            file_list.append(item)
    if sort:
        file_list = sorted(file_list)
    siteUtils.print_file_list(description, file_list)
    return file_list

if __name__ == '__main__':
#    get_met_scan_data('e2v-CCD250-11093-10-04', '*.csv')
    get_met_scan_data('ITL-3800C-033', '*.txt')
