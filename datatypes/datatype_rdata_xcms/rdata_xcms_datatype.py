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
    

class RdataXcmsRaw( Binary ):
    file_ext = "rdata.xcms.raw"
    blurb = "Output of xcms.xcmsSet\nInput of xcms.group and xcms.retcor"
    
    #TODO: sniff


class RdataXcmsGroup( Binary ):
    file_ext = "rdata.xcms.group"
    blurb = "Output of xcms.group\nInput of xcms.group, xcms.retcor and xcms.fillpeaks"
    
    #TODO: sniff


class RdataXcmsRetcor( Binary ):
    file_ext = "rdata.xcms.retcor"
    blurb = "Output of xcms.retcor\nInput of xcms.group"
    
    #TODO: sniff


class RdataXcmsFillpeaks( Binary ):
    file_ext = "rdata.xcms.fillpeaks"
    blurb = "Output of xcms.fillpeaks\nInput of camera.annotate"
    
    #TODO: sniff


