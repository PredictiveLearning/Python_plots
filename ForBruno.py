
# coding: utf-8

# In[1]:

import SciServer.CasJobs as cj


# In[2]:

query = '''select GalaxyID,
        lastProgenitorId
into MyDB.Bruno2
   from MRscPlanck1
  where redshift >3 and stellarmass/0.673>10.
    and Sfr/(stellarmass*1.e10/0.673)< 1.e-10'''

DB = 'Henriques2015a'
jobId = cj.submitJob(query,DB)

cj.waitForJob(jobId)



# In[3]:

df=cj.executeQuery("select * from Bruno1","MyDB")
df


# In[ ]:



