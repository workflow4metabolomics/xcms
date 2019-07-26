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
    

class RdataCameraPositive( Binary ):
    file_ext = "rdata.camera.positive"
    blurb = "Output of CAMERA.annotate using positif acquisition mode\nInput of ProbMetab or CAMERA.combinexcAnnos"
    
    #TODO: sniff


class RdataCameraNegative( Binary ):
    file_ext = "rdata.camera.negative"
    blurb = "Output of CAMERA.annotate using negatif acquisition mode\nInput of ProbMetab or CAMERA.combinexcAnnos"
    
    #TODO: sniff


class RdataCameraQuick( Binary ):
    file_ext = "rdata.camera.quick"
    blurb = "Output of CAMERA.annotate using the quick mode"
    
    #TODO: sniff


