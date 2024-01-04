#!/usr/bin/env python
# coding: utf-8

# In[31]:

 
import csv, time
from itertools import groupby
import os,sys
import astropy.io.ascii as ascii
from astropy.table import Table
from isochrones import StarModel, get_ichrone, SingleStarModel
import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
from astropy import coordinates 
import astropy.units as u
#from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astroquery.gaia import Gaia
Simbad.add_votable_fields('sptype','flux(V)','flux(G)','flux(B)','distance','plx','plx_error','ids')

def GetGAIAData(GaiaDR3SourceIDs):
    # gets the GAIA data for the provided GaiaDR3SourceIDs's
 
    qry = "SELECT * FROM gaiadr3.gaia_source gs WHERE gs.source_id in (" + GaiaDR3SourceIDs + ");"
    
    job = Gaia.launch_job_async( qry )
    tblGaia = job.get_results()       #Astropy table
    return tblGaia
spec_st_param = np.genfromtxt("parameters.csv", dtype=None, delimiter=',', names=True ) 
#spec_st_param.dtype.names
#Simbad.list_votable_fields()


# In[2]:


### HARPS_RVBank ver02
st_name = spec_st_param['Name']

teff    = spec_st_param['teff']
e_teff  = spec_st_param['eteff']

logg    = spec_st_param['logg']
e_logg  = spec_st_param['elogg']
 
feh     = spec_st_param['F']
e_feh   = spec_st_param['eF']


# In[3]:


st_name.shape


# In[4]:


print(type(spec_st_param))


# In[5]:


spec_st_param_df = pd.DataFrame(spec_st_param)


# In[6]:


spec_st_param_df


# In[7]:


spec_st_param_df.shape


# In[8]:


spec_st_param_df.columns.to_list()


# In[9]:


bp   = []
rp   = []
ag   = []
ID   = []
px   = []
e_px = []

for d in range(len(st_name[0:3])):

    print(st_name[d].decode("utf-8"))
    try:
        #### Simbad 
        result_table = Simbad.query_object(st_name[d].decode("utf-8"), wildcard=True)
        
        Ids = Simbad.query_objectids(st_name[d].decode("utf-8"), cache=False)
        Gaia_DR3_id = [i for i in Ids["ID"] if i.startswith('Gaia DR3')][0] #Gets the DR3 ID
        
        ##### Gaia  
        Gaia_table = GetGAIAData(Gaia_DR3_id[8:])
 
        bp.append(Gaia_table['phot_bp_mean_mag'][0])
        rp.append(Gaia_table['phot_rp_mean_mag'][0])
        ID.append(Gaia_DR3_id[8:])
        ag.append(Gaia_table['ag_gspphot'][0]) #https://arxiv.org/pdf/2206.07937.pdf
        px.append(Gaia_table['parallax'][0])
        e_px.append(Gaia_table['parallax_error'][0])
    except:    
        bp.append(np.nan)
        rp.append(np.nan)
        ID.append(np.nan)
        ag.append(np.nan)
        px.append(np.nan)
        e_px.append(np.nan)
        


# In[10]:


print(type(Gaia_table))


# In[11]:


gaia_table_pd = Gaia_table.to_pandas()


# In[12]:


gaia_table_pd


# In[14]:


Gaia_table.dtype.names


# In[15]:


#Ids = Simbad.query_objectids('BD+072474', cache=False)
#Ids


# In[16]:


#Gaia_table.dtype.names



# In[18]:


spec_st_param_df['eep_calc'] = np.nan
spec_st_param_df['eep_e_calc'] = np.nan
spec_st_param_df['age_calc']= np.nan
spec_st_param_df['age_e_calc'] = np.nan
spec_st_param_df['feh_calc'] = np.nan
spec_st_param_df['feh_e_calc'] = np.nan
spec_st_param_df['mass_calc'] = np.nan
spec_st_param_df['mass_e_calc'] = np.nan
spec_st_param_df['initial_mass_calc']= np.nan
spec_st_param_df['initial_mass_e_calc'] = np.nan
spec_st_param_df['radius_calc'] = np.nan
spec_st_param_df['radius_e_calc'] = np.nan
spec_st_param_df['density_calc'] = np.nan
spec_st_param_df['density_e_calc'] = np.nan
spec_st_param_df['logTeff_calc']= np.nan
spec_st_param_df['logTeff_e_calc'] = np.nan
spec_st_param_df['Teff_calc'] = np.nan
spec_st_param_df['Teff_e_calc'] = np.nan
spec_st_param_df['logg_calc'] = np.nan
spec_st_param_df['logg_e_calc'] = np.nan
spec_st_param_df['logL_calc']= np.nan
spec_st_param_df['logL_e_calc'] = np.nan
spec_st_param_df['Mbol_calc'] = np.nan
spec_st_param_df['Mbol_e_calc'] = np.nan
spec_st_param_df['delta_nu_calc']= np.nan
spec_st_param_df['delta_nu_e_calc'] = np.nan
spec_st_param_df['nu_max_calc'] = np.nan
spec_st_param_df['nu_max_e_calc'] = np.nan
spec_st_param_df['phase_calc'] = np.nan
spec_st_param_df['phase_e_calc'] = np.nan
spec_st_param_df['dm_deep_calc']= np.nan
spec_st_param_df['dm_deep_e_calc'] = np.nan
spec_st_param_df['BP_mag_calc'] = np.nan
spec_st_param_df['BP_mag_e_calc'] = np.nan
spec_st_param_df['RP_mag_calc']= np.nan
spec_st_param_df['RP_mag_e_calc'] = np.nan
spec_st_param_df['parallax_calc'] = np.nan
spec_st_param_df['parallax_e_calc'] = np.nan
spec_st_param_df['distance_calc'] = np.nan
spec_st_param_df['distance_e_calc'] = np.nan
spec_st_param_df['AV_calc']= np.nan
spec_st_param_df['AV_e_calc'] = np.nan

 

check_params = {'Teff': (teff[0],e_teff[0]), 'logg': (logg[0], e_logg[0]), 'feh': (feh[0], e_feh[0]),'BP': (bp[0], 0.001), 
              'RP': (rp[0], 0.001),'parallax': (px[0], e_px[0])}#, 'AV': (0.101,0.009)} 


 

# In[19]:


mist = get_ichrone('mist', bands=['BP', 'RP'])
for i in range(len(st_name)):

    print("isoChrones running on:", st_name[i])
    all_params = {'Teff': (teff[i],e_teff[i]), 'logg': (logg[i], e_logg[i]), 'feh': (feh[i], e_feh[i]),'BP': (bp[i], 0.001), 
              'RP': (rp[i], 0.001),'parallax': (px[i], e_px[i])}#, 'AV': (0.101,0.009)} 

    params = {}
    for key,value in all_params.items():
        if all(el is not None for el in value):
            params[key] = value
            print("Adding {} to table".format(key))
        else:
            print("Param NA")
        
    mod = SingleStarModel(mist,**params)        

    mod.fit_multinest(n_live_points=1000,
    basename=None,
    verbose=False,
    refit=True,
    overwrite=True)

    derived_samples = mod.derived_samples.describe()

    
    spec_st_param_df['eep_calc'].iloc[i] = derived_samples['eep'][1]
    spec_st_param_df['eep_e_calc'].iloc[i] = derived_samples['eep'][2]
    spec_st_param_df['age_calc'].iloc[i] = derived_samples['age'][1]
    spec_st_param_df['age_e_calc'].iloc[i] = derived_samples['age'][2]
    spec_st_param_df['feh_calc'].iloc[i] = derived_samples['feh'][1]
    spec_st_param_df['feh_e_calc'].iloc[i] = derived_samples['feh'][2]
    spec_st_param_df['mass_calc'].iloc[i] = derived_samples['mass'][1]
    spec_st_param_df['mass_e_calc'].iloc[i] = derived_samples['mass'][2]
    spec_st_param_df['initial_mass_calc'].iloc[i] = derived_samples['initial_mass'][1]
    spec_st_param_df['initial_mass_e_calc'].iloc[i] = derived_samples['initial_mass'][2]
    spec_st_param_df['radius_calc'].iloc[i] = derived_samples['radius'][1]
    spec_st_param_df['radius_e_calc'].iloc[i] = derived_samples['radius'][2]
    spec_st_param_df['density_calc'].iloc[i] = derived_samples['density'][1]
    spec_st_param_df['density_e_calc'].iloc[i] = derived_samples['density'][2]
    spec_st_param_df['logTeff_calc'].iloc[i] = derived_samples['logTeff'][1]
    spec_st_param_df['logTeff_e_calc'].iloc[i] = derived_samples['logTeff'][2]
    spec_st_param_df['Teff_calc'].iloc[i] = derived_samples['Teff'][1]
    spec_st_param_df['Teff_e_calc'].iloc[i] = derived_samples['Teff'][2]
    spec_st_param_df['logg_calc'].iloc[i] = derived_samples['logg'][1]
    spec_st_param_df['logg_e_calc'].iloc[i] = derived_samples['logg'][2]
    spec_st_param_df['logL_calc'].iloc[i] = derived_samples['logL'][1]
    spec_st_param_df['logL_e_calc'].iloc[i] = derived_samples['logL'][2]
    spec_st_param_df['Mbol_calc'].iloc[i] = derived_samples['Mbol'][1]
    spec_st_param_df['Mbol_e_calc'].iloc[i] = derived_samples['Mbol'][2]
    spec_st_param_df['delta_nu_calc'].iloc[i] = derived_samples['delta_nu'][1]
    spec_st_param_df['delta_nu_e_calc'].iloc[i] = derived_samples['delta_nu'][2]
    spec_st_param_df['nu_max_calc'].iloc[i] = derived_samples['nu_max'][1]
    spec_st_param_df['nu_max_e_calc'].iloc[i] = derived_samples['nu_max'][2]
    spec_st_param_df['phase_calc'].iloc[i] = derived_samples['phase'][1]
    spec_st_param_df['phase_e_calc'].iloc[i] = derived_samples['phase'][2]
    spec_st_param_df['dm_deep_calc'].iloc[i] = derived_samples['dm_deep'][1]
    spec_st_param_df['dm_deep_e_calc'].iloc[i] = derived_samples['dm_deep'][2]
    spec_st_param_df['BP_mag_calc'].iloc[i] = derived_samples['BP_mag'][1]
    spec_st_param_df['BP_mag_e_calc'].iloc[i] = derived_samples['BP_mag'][2]
    spec_st_param_df['RP_mag_calc'].iloc[i] = derived_samples['RP_mag'][1]
    spec_st_param_df['RP_mag_e_calc'].iloc[i] = derived_samples['RP_mag'][2]
    spec_st_param_df['parallax_calc'].iloc[i] = derived_samples['parallax'][1]
    spec_st_param_df['parallax_e_calc'].iloc[i] = derived_samples['parallax'][2]
    spec_st_param_df['distance_calc'].iloc[i] = derived_samples['distance'][1]
    spec_st_param_df['distance_e_calc'].iloc[i] = derived_samples['distance'][2]
    spec_st_param_df['AV_calc'].iloc[i] = derived_samples['AV'][1]
    spec_st_param_df['AV_e_calc'].iloc[i] = derived_samples['AV'][2]


# In[22]:

 


# In[74]:


spec_st_param_df.head(10)


# In[78]:


spec_st_param_df.Teff_calc[3]


# In[ ]:




