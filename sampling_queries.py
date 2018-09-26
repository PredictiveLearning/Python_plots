
# coding: utf-8

# In[1]:

import SciServer.CasJobs as cj

query = '''select GalaxyID,
        lastProgenitorId
into MyDB.Bruno1
   from MRscPlanck1
  where redshift >3 and stellarmass/0.673>10.
    and Sfr/(stellarmass*1.e10/0.673)< 1.e-10'''

DB = 'Henriques2015a'
jobId = cj.submitJob(query,DB)

cj.waitForJob(jobId)

df=cj.executeQuery("select * from Bruno1","MyDB")


# In[ ]:



