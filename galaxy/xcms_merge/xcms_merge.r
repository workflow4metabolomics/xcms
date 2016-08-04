rm(xset_merged);for(image in images) { load(image); if (!exists("xset_merged")) xset_merged=xset else xset_merged=c(xset_merged,xset); print(xset) }
