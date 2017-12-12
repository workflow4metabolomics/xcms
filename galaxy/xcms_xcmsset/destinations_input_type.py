from galaxy.jobs import JobDestination
import logging
import os

'''
This file must be placed in lib/galaxy/jobs/rules/
'''

def input_type(job):
    '''
    This function checks the input file format/extension and decide which destination in the job_conf.xml using
     - If it's a zip file, we will launch the job in multi-thread mode (-pe thread 8)
     - If it's an individual file (mzxml, mzml, mzdata or netcdf), the job will use for instance (-pe thread 1)
    '''
    log = logging.getLogger(__name__)
    inp_data = dict( [ ( da.name, da.dataset ) for da in job.input_datasets ] )
    inp_data.update( [ ( da.name, da.dataset ) for da in job.input_library_datasets ] )
    # for the backward compatibility < 2.1.0 
    if 'input' in inp_data:
        input_extension = inp_data[ "input" ].extension
        log.debug("The input extension is %s" % input_extension)
        if input_extension in ["mzxml","mzml","mzdata","netcdf"]:
            return 'thread1-mem_free8'
    # zip file
    return 'thread9-mem_free8'
