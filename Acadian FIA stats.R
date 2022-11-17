#**********************************************
##      Introduction                              ----
# This script will assess the proportion of FIA plots harvested in a region
# and the covariates that can be applied as filters to impact that proportion.
#      Last edited: EPB, 8-9-21

#**********************************************
rm(list = ls())
library(rFIA)
library(tidyverse)
library(dplyr)
library(tibble)
library(geosphere)
detach("package:plyr", unload=TRUE)

path.FIA.storage = "C:\\Users\\ethan.belair.TNC\\Documents\\R\\FIA Data"
path.FIA.storage.modified = "C:\\Users\\ethan.belair.TNC\\Documents\\R\\FIA Data\\Modified FIA Tables"
path.input = "C:\\Users\\ethan.belair.TNC\\Box\\Count\\Data\\"


#**********************************************
##      Import Data                               ----
#**********************************************

# FIPS codes and region codes
fips = read.csv(paste0(path.input, "Misc\\StateFIPS.csv"))
fips$State = as.factor(fips$State)
fips$FRA_Region = as.factor(fips$FRA_Region)
fips$FRA_Region_abbr = as.factor(fips$FRA_Region_abbr)

# Reference table containing
# FIA REF_SPECIES table
# Density Reduction Factors from Harmon et al. 2011
# Storage Factors from ???
ref.table = read.csv(paste0(path.input, "REF Table.csv"), na.strings = c("", "#N/A"))

# FIA PLOT and COND tables, for linking covariates
# Download FIA data 
data = readFIA(states = c('ME', 'NH', 'VT'), 
               dir = path.FIA.storage,
               tables = c("PLOT", "COND", "TREE", "POP_STRATUM", "POP_PLOT_STRATUM_ASSGN")) # 

##Data reset/backup
#data.temp = data
#data = data.temp


#**********************************************
##      Process Data                              ----
#**********************************************

# Append Specific Gravity from REF_SPECIES to tree table
data$TREE$SG = ref.table[match(data$TREE$SPCD, ref.table$SPCD), "WOOD_SPGR_GREENVOL_DRYWT"]

## Stem Biomass
data$TREE$stem.bio.ac = data$TREE$VOLCFGRS * data$TREE$TPA_UNADJ * data$TREE$SG

## Boardfoot volume per acre
data$TREE$volbfnet.ac = data$TREE$VOLBFNET * data$TREE$TPA_UNADJ
data$TREE$volbfnet.ac.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "volbfnet.ac"]

## Gross cubic volume per acre
data$TREE$volcfgrs.ac = data$TREE$VOLCFGRS * data$TREE$TPA_UNADJ

## Sound merchantable cubic volume per acre
data$TREE$volcfsnd.ac = data$TREE$VOLCFSND * data$TREE$TPA_UNADJ

## C per acre - above - and below-ground
data$TREE$c.per.ac = (data$TREE$CARBON_AG * data$TREE$TPA_UNADJ) #+ (data$TREE$CARBON_BG * data$TREE$TPA_UNADJ)

data$TREE$t.co2e.per.ha = data$TREE$c.per.ac * 2.47 *(44/12) * (1/2200)


##### Add previous plot CN to allow later removal of older plot measurements  
# Add the previous plot sequence to each tree from the plot table
# this will be used to link individual trees and plot level summaries to their previous measurements
data$TREE$PREV_PLT_CN = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "PREV_PLT_CN"]
data$TREE$prev.prev.plt.cn = data$PLOT[match(data$TREE$PREV_PLT_CN, data$PLOT$CN), "PREV_PLT_CN"]





#**********************************************
##      Subset Data                               ----
#**********************************************

## Remove periodic inventory points
# Append KINDCD to tree and cond tables
data$TREE$KINDCD = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "KINDCD"]
data$COND$KINDCD = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "KINDCD"]
# Subset by KINDCD
data$PLOT = data$PLOT[which(data$PLOT$KINDCD %in% c(1,2,3)), ]
data$TREE = data$TREE[which(data$TREE$KINDCD %in% c(1,2,3)), ]
data$COND = data$COND[which(data$COND$KINDCD %in% c(1,2,3)), ]


## Remove split condition plots
# This will inevitably cause issues with plots that are "unsplit" at one meausrement and "split" at another
data$COND = data$COND[which(data$COND$CONDPROP_UNADJ == 1), ]
data$PLOT = data$PLOT[which(data$PLOT$CN %in% unique(data$COND$PLT_CN)),]
data$TREE = data$TREE[which(data$TREE$PLT_CN %in% unique(data$COND$PLT_CN)),]



## Subset to appropriate region for analysis (Acadian Forest region)
# Subset PLOT table

# Create temp layers from the PLOT table for each state, containing only the selected counties FFCP operates in
temp.me = subset(data$PLOT, STATECD == 23)
temp.me = subset(temp.me, COUNTYCD %in% c(3,7,9,17,19,21,25,29))

temp.nh = subset(data$PLOT, STATECD == 33)
temp.nh = subset(temp.nh, COUNTYCD %in% c(3,5,7,9,19))

temp.vt = subset(data$PLOT, STATECD == 50)
temp.vt = subset(temp.vt, COUNTYCD %in% c(3,5,9,13,15,17,19,21,23,25))


# combine the separate temp layers
plot = rbind( temp.me, temp.nh, temp.vt)

# Replace the PLOT table with this output from rbind
data$PLOT = plot

# Remove temp layers
rm(temp.me, temp.nh, temp.vt, plot)

# subset COND and TREE tables
data$COND = data$COND[which(data$COND$PLT_CN %in% unique(data$PLOT$CN)),]
data$TREE = data$TREE[which(data$TREE$PLT_CN %in% unique(data$PLOT$CN)),]


## Remove all plots which are not the most recent measurement
# The variables prev.plt.cn and prev.plt.cn were previously added to each record in the tree table.
# This variable indicates the previous plot sequence number for all plots which were measured in previous years.
# Thus, if a plot is listed as another plot's prev.plt.cn, it is ineligible for inclusion

# Create a vector including all prev.plt.cn values
prev.plt.cn = unique(data$TREE$PREV_PLT_CN)

# add variable to tree, plot, and cond tables indicating whether a plot's CN not listed in prev.plt.cn
data$TREE$prev.meas = ifelse(data$TREE[,'PLT_CN'] %in% prev.plt.cn, 1, 0)
data$PLOT$prev.meas = ifelse(data$PLOT[,'CN'] %in% prev.plt.cn, 1, 0)
data$COND$prev.meas = ifelse(data$COND[,'PLT_CN'] %in% prev.plt.cn, 1, 0)


data$TREE = data$TREE[which(data$TREE$prev.meas == 0), ]
data$PLOT = data$PLOT[which(data$PLOT$prev.meas == 0), ]
data$COND = data$COND[which(data$COND$prev.meas == 0), ]


# Repeat process for all prev.prev.plt.cn values
prev.prev.plt.cn = unique(data$TREE$prev.prev.plt.cn)

# add variable to tree, plot, and cond tables indicating whether a plot's CN not listed in prev.prev.plt.cn
data$TREE$prev.prev.meas = ifelse(data$TREE[,'PLT_CN'] %in% prev.prev.plt.cn, 1, 0)
data$PLOT$prev.prev.meas = ifelse(data$PLOT[,'CN'] %in% prev.prev.plt.cn, 1, 0)
data$COND$prev.prev.meas = ifelse(data$COND[,'PLT_CN'] %in% prev.prev.plt.cn, 1, 0)


data$TREE = data$TREE[which(data$TREE$prev.prev.meas == 0), ]
data$PLOT = data$PLOT[which(data$PLOT$prev.prev.meas == 0), ]
data$COND = data$COND[which(data$COND$prev.prev.meas == 0), ]

# remove vectors
rm(prev.plt.cn, prev.prev.plt.cn) 





#**********************************************
##      Analyse Data                               ----
#**********************************************

## Percent of Maine privately owned

# create maine dataset
maine = data
maine$PLOT = maine$PLOT[which(maine$PLOT$STATECD == 23), ]
maine$TREE = maine$TREE[which(maine$TREE$STATECD == 23), ]
maine$COND = maine$COND[which(maine$COND$STATECD == 23), ]

table(maine$COND$OWNGRPCD)
2103/(2103+20+23+150)

table(data$COND$OWNGRPCD)
2710/(2710+492+34+226)



## Average Stocking on Acadian Forest region
tree = data$TREE[which(data$TREE$DIA >= 5), ]
summ = tree %>%
  group_by(PLT_CN, STATECD, COUNTYCD, PLOT, INVYR) %>%  
  summarize(vol.cf.grs = sum(volcfgrs.ac, na.rm = TRUE),
            vol.cf.snd = sum(volcfsnd.ac, na.rm = TRUE))

cu.ft.ac = mean(summ$vol.cf.snd) #1880 cubic feet per acre
cords = cu.ft.ac/85 #22.1 cords/acre
cu.m.ha = (cu.ft.ac / 35.15) * 2.47 # 132 cubic meters per ha



## Estimate 2 of C per ha for high stocking stands
tree = data$TREE[which(data$TREE$DIA >= 5), ]
summ = tree %>%
  group_by(PLT_CN, STATECD, COUNTYCD, PLOT, INVYR) %>%  
  summarize(vol.cf.grs = sum(volcfgrs.ac, na.rm = TRUE),
            vol.cm.snd = sum(volcfsnd.ac, na.rm = TRUE)/35.15,
            c.per.ha = sum(t.co2e.per.ha, na.rm = TRUE))

summ.est2 = summ[which(summ$vol.cm.snd >= 143 & summ$vol.cm.snd < 155), ]
mean(summ.est2$c.per.ha)
 