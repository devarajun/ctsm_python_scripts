#!/usr/bin/python
#***********************************************************************
#Create new lat, lon grids and overlay the the original data from METNorDIC
#***********************************************************************

#import required modules
import datetime, netCDF4, numpy as np, os, numpy.ma as ma
import glob
import xarray as xr
import pandas as pd
#from scipy.interpolate import griddata
#from scipy.ndimage.interpolation import shift

#list site names in respective forcing files
sitename = ['fns_v1','isk_v1','aas_v1','hisl_v1']
metsitename = ['Finse','Iskoras','Aas','Hisaasen_lower']

gswp_varname = ['FSDS','PRECTmms','PSRF','TBOT','WIND','QBOT','FLDS']
metnord_varname = ['DIR_SWdown','Rainf','Snowf','PSurf','Tair','Wind','Qair','LWdown']

#directories
gswp_dir  = '/cluster/projects/nn2806k/dev/sites_forcing/'
metnord_dir='/tos-project3/NS9066K/ET_MIP/forc/METNordicpy/'

##*************************************************************************************
#output files creation
# Define global attributes to be written in the output netcdf-files
outdir = '/cluster/projects/nn2806k/dev/sites_forcing/metnordic_clm'
atts = dict(web = "--",date = datetime.datetime.now().strftime("%d/%m/%Y %H:%M"),
            resolution = "Regular lat/lon single point")


#loop through the site and extract variables, write variables
i=0
for site in sitename:
        i = i+1
        print('GSWP3 site',site,i)
        listfiles_solr = glob.glob(gswp_dir+site+"/clmforc.GSWP3.c2011.0.5x0.5.Solr*.1972-01.nc")
        #listfiles_pr   = glob.glob(gswp_dir+site+"/clmforc.GSWP3.c2011.0.5x0.5.Prec*.1972-01.nc")
        #listfiles_tpq  = glob.glob(gswp_dir+site+"/clmforc.GSWP3.c2011.0.5x0.5.TPQWL*.1972-01.nc")
        solrf  = xr.open_dataset(listfiles_solr[0])
        #prf   = xr.open_dataset(listfiles_pr[0])
        #tpqf  = xr.open_dataset(listfiles_tpq[0])
        metncin = xr.open_dataset(metnord_dir+'FORCING_'+metsitename[i-1]+'_201410_to_202110_update.nc')
        print('METNordic site',metsitename[i-1])
        gswp_fsds = solrf['FSDS']
        met_fsds  = metncin['DIR_SWdown']
        met_prect = metncin['total_precipitation']
        met_psrf  = metncin['PSurf']
        met_tbot  = metncin['Tair']
        met_wind  = metncin['Wind']
        met_qbot  = metncin['Qair']
        met_flds  = metncin['LWdown']
        print(gswp_fsds.shape)
        print(met_fsds.shape)
        met_fsds_nparray = np.zeros((len(met_fsds),1,1))
        met_fsds_nparray[:,0,0] = np.array(met_fsds.data,dtype='float32')
        met_prect_nparray = np.zeros((len(met_prect),1,1))
        met_prect_nparray[:,0,0]=np.array(met_prect.data,dtype='float32')
        met_psrf_nparray = np.zeros((len(met_psrf),1,1))
        met_psrf_nparray[:,0,0]=np.array(met_psrf.data,dtype='float32')
        met_tbot_nparray = np.zeros((len(met_tbot),1,1))
        met_tbot_nparray[:,0,0]=np.array(met_tbot.data,dtype='float32')
        met_wind_nparray = np.zeros((len(met_wind),1,1))
        met_wind_nparray[:,0,0]=np.array(met_wind.data,dtype='float32')
        met_qbot_nparray = np.zeros((len(met_qbot),1,1))
        met_qbot_nparray[:,0,0]=np.array(met_qbot.data,dtype='float32')
        met_flds_nparray = np.zeros((len(met_flds),1,1))
        met_flds_nparray[:,0,0]=np.array(met_flds.data,dtype='float32')
        #time_values = np.array(met_fsds.time.data,dtype='float32')
        #initialize dataset with multiple dimensions
        ds = xr.Dataset(
             data_vars=dict(
                 FSDS=(["time","lat","lon"], met_fsds_nparray,{'units':'W/m**2','long_name':'total incident solar radiation'}),
                 PRECTmms=(["time","lat","lon"], met_prect_nparray,{'units':'W/m2','long_name':'total precipitation'}),
                 PSRF=(["time","lat","lon"], met_psrf_nparray,{'units':'Pa','long_name':'surface pressure at the lowest atm level'}),
                 TBOT=(["time","lat","lon"], met_tbot_nparray,{'units':'kg/kg','long_name':'temperature at the lowest atm level'}),
                 WIND=(["time","lat","lon"], met_wind_nparray,{'units':'m/s','long_name':'wind at the lowest atm level'}),
                 QBOT=(["time","lat","lon"], met_qbot_nparray,{'units':'kg/kg','long_name':'specific humidity at the lowest atm level'}),
                 FLDS=(["time","lat","lon"], met_flds_nparray,{'units':'W/m**2','long_name':'incident longwave radiation'})),
                 coords= dict(time=(["time"],met_fsds.time.data),lat=solrf['lat'],lon=solrf.lon),
                 attrs=atts,
        )
        ds_1=xr.merge([ds,solrf.EDGEE,solrf.EDGEW,solrf.EDGES,solrf.EDGEN,solrf.LONGXY,solrf.LATIXY])
        #write the dataset to netcdf file
        ds_1.to_netcdf('clmforc.METNordic_'+metsitename[i-1]+'_201410_to_202110.nc',mode='w')

