# -*- coding: utf-8 -*-
"""
xcms datatypes

Author: lecorguille@sb-roscoff.fr
"""

import logging
import os,os.path,re
from galaxy.datatypes.data import *
from galaxy.datatypes.sniff import *
from galaxy.datatypes.binary import * 

log = logging.getLogger(__name__)
    

class RdataXcmsGeneric (RData):
    """Class describing a sequence"""
    edam_data = "data_0943"
    
    """Add metadata elements"""
    MetadataElement( name="xset", default=0, desc="XCMS object", readonly=True, visible=False, optional=True, no_value=0 )

    def set_meta( self, dataset, **kwd ):
        """
        Set XCMS object information in dataset.
        
        Example:
            An "xcmsSet" object with 4 samples

            Time range: 2506-4484 seconds (41.8-74.7 minutes)
            Mass range: 200.1-600 m/z
            Peaks: 32720 (about 8180 per sample)
            Peak Groups: 8157
            Sample classes: KO, WT

            Profile settings: method = bin
                              step = 0.1

            Memory usage: 4.25 MB

        """
        try:
            cmd="R --slave -e \"load('%s'); print(xset); \" 2> /dev/null"%(file_name)
            xset=subprocess.check_output(cmd, shell=True)
        except:
            xset="Either the RData don't contain a xset object or R is not installed"
        dataset.metadata.xset=xset

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            
            if dataset.metadata.xset:
                dataset.peek = str( dataset.metadata.xset )
            else:
                dataset.peek = ""
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

class RdataXcmsRaw( RdataXcmsGeneric ):
    file_ext = "rdata.xcms.raw"
    dataset.blurb = "Output of xcms.xcmsSet\nInput of xcms.group and xcms.retcor"
    
    #TODO: sniff


class RdataXcmsGroup( RdataXcmsGeneric ):
    file_ext = "rdata.xcms.group"
    dataset.blurb = "Output of xcms.group\nInput of xcms.group, xcms.retcor and xcms.fillpeaks"
    
    #TODO: sniff


class RdataXcmsRetcor( RdataXcmsGeneric ):
    file_ext = "rdata.xcms.retcor"
    dataset.blurb = "Output of xcms.retcor\nInput of xcms.group"
    
    #TODO: sniff


class RdataXcmsFillpeaks( RdataXcmsGeneric ):
    file_ext = "rdata.xcms.fillpeaks"
    dataset.blurb = "Output of xcms.fillpeaks\nInput of camera.annotate"
    
    #TODO: sniff


