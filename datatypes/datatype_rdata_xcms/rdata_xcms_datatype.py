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


class RdataMsnbaseRaw( RData ):
    file_ext = "rdata.msnbase.raw"
    blurb = "Output of MSnbase.readMSData\nInput of xcms.findChromPeaks"

    #TODO: sniff

class RdataXcmsFindChromPeaks( RData ):
    file_ext = "rdata.xcms.findchrompeaks"
    blurb = "Output of xcms.findChromPeaks\nInput of xcms.groupChromPeaks and xcms.adjustRtime"

    #TODO: sniff

class RdataXcmsGroup( RData ):
    file_ext = "rdata.xcms.group"
    blurb = "Output of xcms.groupChromPeaks\nInput of xcms.groupChromPeaks, xcms.adjustRtime, xcms.fillChromPeaks and camera.annotate"

    #TODO: sniff


class RdataXcmsRetcor( RData ):
    file_ext = "rdata.xcms.retcor"
    blurb = "Output of xcms.adjustRtime\nInput of xcms.groupChromPeaks"

    #TODO: sniff


class RdataXcmsFillpeaks( RData ):
    file_ext = "rdata.xcms.fillpeaks"
    blurb = "Output of xcms.fillChromPeaks\nInput of camera.annotate"

    #TODO: sniff
