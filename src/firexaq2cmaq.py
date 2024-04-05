# # FIREX-AQ 2 CMAQ
# ### Grid FIREX-AQ observations to CMAQ grid

# Inputs:
# * FIREX-AQ data
# * Met files
# * GRIDDESC file

# James East
# 2022-08-09

import PseudoNetCDF as pnc
import pandas as pd
import numpy as np
import xarray as xr
import sys

startdate = sys.argv[1] # YYYYMMDD
enddate   = sys.argv[2] # YYYYMMDD
mygrid    = sys.argv[3] # 108NHEMI2
#mygrid    = '108NHEMI2'
plumet    = sys.argv[4] # 10
#plumet    = 10 # near plume minutes to exclude

gf = pnc.pncopen('/home/jeast/GRIDDESC',format='griddesc',GDNAM=mygrid)

# path, can read this from file/command line
pathtofxaq = '../data/MERGES/60_SECOND.DC8_MRG/firexaq-mrg60-dc8_merge_20190722_R1_thru20190905.ict'
fxaq = pnc.pncopen(pathtofxaq,format='ffi1001')

# list all variables and pick ones I want
keepv = []
keepvu = []
for v in list(fxaq.variables.keys()):
    #if ('NO' in v) or ('O3' in v):
    if 'RYERSON' in v:
        print(v)
        keepv.append(v)
        keepvu.append(f'{v}_{fxaq.variables[v].units}')
    if 'CO_' in v and 'DISKIN' in v:
        print(v)
        keepv.append(v)
        keepvu.append(f'{v}_{fxaq.variables[v].units}')
altitudevar = 'MSL_GPS_Altitude_YANG'
print('Keeping variables: ',keepv,flush=True)

firstdate = pd.to_datetime(startdate, format='%Y%m%d')
lastdate = pd.to_datetime(enddate, format='%Y%m%d')
drange = pd.date_range(start=firstdate, end=lastdate, freq='D')
#sdate = pd.to_datetime(f'{m3d.SDATE}',format='%Y%j') # assumes STIME = 0

def grid1day(sdate):
    '''Grid 1 day of FIREX-AQ to CMAQ'''
    # create a dataframe of FRIEXAQ data
    dfx = pd.DataFrame()
    dfx['time'] = pd.to_datetime('20190101')+pd.to_timedelta(fxaq.variables['Fractional_Day']-1,unit='D')
    dfx['time'] = dfx['time'].dt.round('1min')
    dfx['lat'] = fxaq.variables['Latitude_YANG'].array()
    dfx['lon'] = fxaq.variables['Longitude_YANG'].array()
    dfx['smokeflag'] = fxaq.variables['Smoke_flag_SCHWARZ'].array()
    dfx[f'altitude_{fxaq.variables[altitudevar].units}'] = fxaq.variables[altitudevar].array()
    for v,vu in zip(keepv,keepvu):
        dfx[vu] = fxaq.variables[v].array()

    # lat/lon to CMAQ grid J/I
    dfx['I'],dfx['J'] = gf.ll2ij(dfx['lon'],dfx['lat'])

    # time to CMAQ TSTEP -- floor divide
    dfx['T'] = np.array([(t - sdate).total_seconds() // 3600 for t in dfx['time']], dtype='i')

    # subset only the day we want
    dfxd = dfx.query('T >= 0 and T <= 23').copy(deep=True)
    if dfxd.shape[0] == 0:
        print(f'No observations on {sdate.strftime("%Y%m%d")}')
        return

    # assign to CMAQ layer based on altitude and ZH variable
    m3df = sdate.strftime('/work/MOD3EVAL/jeast/NO2ASSIM/CMAQ/input/2019_hemi/mcip/METCRO3D_%Y%m%d.nc4')
    m3d = xr.open_dataset(m3df)
    dfxd['K'] = np.ma.masked_less(
        m3d['ZF'].values[dfxd['T'],:,dfxd['J'],dfxd['I']] - dfxd['altitude_m'].to_numpy()[:,None],0
    ).argmin(1)

    # create flags for in/out of plume
    dfxd['out_plume'] = np.where(dfxd['smokeflag'].isna(),True,False)
    dfxd['in_plume'] = ~dfxd['out_plume'].copy(deep=True)

    # find a better way to do this
    # identify all minutes within "plumet" minutes of being in a plume
    tinvalid_range = np.array([
        pd.date_range(
            start = dfxd.where(dfxd['in_plume']).dropna()['time'][i]-pd.to_timedelta(f'{plumet}min'),
            end = dfxd.where(dfxd['in_plume']).dropna()['time'][i]+pd.to_timedelta(f'{plumet}min'),
            freq='1min'
        )
        for i in dfxd.where(dfxd['in_plume']).dropna()['time'].index
    ])

    # keep only the unique minutes
    tinvalid = np.unique(tinvalid_range)

    # all points outside plume but within X minutes of plume
    dfxd['near_plume'] = dfxd['out_plume'] & dfxd['time'].isin(tinvalid) 

    # only keep out of plume and not near plume points
    #dtmp = dfxd.query('out_plume & ~near_plume').groupby(by=['T','K','J','I']).mean()[keepv]

    print(dfxd)

    return dfxd


def save1day(dtmp, sdate, fnameappend=None):
    '''Save 1 day of FIREX AQ to IOAPI'''
    # open an IOAPI for a template
    m3df = sdate.strftime('/work/MOD3EVAL/jeast/NO2ASSIM/CMAQ/input/2019_hemi/mcip/METCRO3D_%Y%m%d.nc4')
    ioapif = pnc.pncopen(m3df, format='ioapi')

    # Make data the shape of an ioapi file
    blankij=pd.DataFrame(
        np.zeros((24*ioapif.dimensions['LAY'].size*ioapif.dimensions['ROW'].size*ioapif.dimensions['COL'].size)),
        columns=['blank'],
        index=pd.MultiIndex.from_product(
            [
                np.arange(0,24,dtype=int),
                np.arange(0,ioapif.dimensions['LAY'].size,dtype=int),
                np.arange(0,ioapif.dimensions['ROW'].size,dtype=int),
                np.arange(0,ioapif.dimensions['COL'].size,dtype=int)
            ],
            names=['T','K','J','I']
        ),
    )
    
    dfout = (
        blankij.drop(['blank'],axis='columns')
    ).merge(
        dtmp,
        how='outer',
        left_index=True,
        right_on=['T','K','J','I']
    )
    dout = dfout.to_xarray()
    
    # create template
    tmpf = ioapif.subsetVariables(['ZF']).copy()
    tmpf.TSTEP = 10000
    setattr(tmpf,'SDATE',int(sdate.strftime('%Y%j')))
    tmpf.updatetflag(overwrite=True)
    outf = tmpf.sliceDimensions(TSTEP=slice(0,24))
    
    # make variable in new file
    for v in dout.data_vars: 
        print(v)
        evar = outf.createVariable(v, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))
        evar.setncatts(dict(
          units=v.split('_')[-1], long_name=v, var_desc=f'FIREX-AQ variable {v} from {pathtofxaq.split("/")[-1]}'
        ))
        evar[:] = dout[v].data
    
    # Get rid of FAKE file variables
    del(outf.variables['ZF'])
    
    # Update TFLAG to be consistent with variables
    outf.updatetflag(tstep=10000, overwrite=True)
    
    # Remove VAR-LIST so that it can be inferred
    delattr(tmpf, 'VAR-LIST')
    setattr(outf,'VAR-LIST','    '.join([v for v in dout.data_vars]))
    outf.updatemeta()
    setattr(outf,'FILEDESC',f'FIREX-AQ gridded to CMAQ {mygrid}')
    
    # save the file to disk
    if fnameappend == None:
        outf.save(f'./gridded/FIREX-AQ.{mygrid}.{sdate.strftime("%Y%m%d")}.nc', format='NETCDF4_CLASSIC')
    else:
        outf.save(f'./gridded/FIREX-AQ.{mygrid}.{fnameappend}.{sdate.strftime("%Y%m%d")}.nc', format='NETCDF4_CLASSIC')
    

if __name__ == '__main__':
    for sdate in drange:
        print(sdate.strftime('%Y%m%d'))
        dfxd = grid1day(sdate) # grid 1 day
        if dfxd is None:
            continue

        # only keep out of plume and not near plume points
        doutplume = dfxd.query('out_plume & ~near_plume').groupby(by=['T','K','J','I']).mean()[keepvu]
        save1day(doutplume, sdate, fnameappend='outofplume') # save file
        
        # only near plume
        dnearplume = dfxd.query('near_plume').groupby(by=['T','K','J','I']).mean()[keepvu]
        save1day(dnearplume, sdate, fnameappend='nearplume') # save file

        # only in plume
        dinplume = dfxd.query('in_plume').groupby(by=['T','K','J','I']).mean()[keepvu]
        save1day(dinplume, sdate, fnameappend='inplume') # save file
    
