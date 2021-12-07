import os
import stochpy
import random
import numpy as numpy
from scipy import stats

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
n_samples=end_time*60
nruns=1000
cases = numpy.zeros([nruns, 5])

#con_traj=numpy.empty([end_time,nruns])
acq_traj=numpy.empty([n_samples,nruns])

# Run is a single run of the model that returns the number of incident cases
def Baserun(iterations):
    for k in range(0, iterations):
        #pdict = {'epsilon1':random.uniform(0.01,0.39),
             #'epsilon2':random.uniform(0.1,48)
            #} 
        model = stochpy.SSA()
        model.Model(model_file='MRSA_model.psc', dir=workingdir)
        model.Endtime(end_time)
        #model.ChangeParameter('epsilon1',pdict['epsilon1'])
        #model.ChangeParameter('epsilon2',1/pdict['epsilon2'])
        model.DoStochSim()
        model.GetRegularGrid(n_samples)
        outcomes = model.data_stochsim_grid.species
        cases[k,0] = outcomes[16][0][-1]
        cases[k,1] = outcomes[17][0][-1]
        cases[k,2] = outcomes[28][0][-1]
        cases[k,3] = outcomes[36][0][-1]
        cases[k,4] = outcomes[42][0][-1]
        for t in range(0,n_samples):
            acq_traj[t,k]=outcomes[16][0][t]
           #con_traj[t,k]=outcomes[16][0][t]+outcomes[19][0][t]+outcomes[21][0][t]+outcomes[23][0][t]+outcomes[25][0][t]+outcomes[27][0][t]+outcomes[35][0][t]+outcomes[37][0][t]+outcomes[38][0][t]+outcomes[39][0][t]+outcomes[40][0][t]+outcomes[41][0][t]
    return cases
  
#Fit delta
base_sweep = Baserun(nruns)
print("Complete")

print("Saving Files")
numpy.savetxt('acq_trajNoColAdm12.csv',acq_traj,delimiter=',')