import sys
import dadi
import numpy
import scipy
import pyOpt
import dadiFunctions

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
except:
    raise ImportError('mpi4py is required for parallelization')

def IM_sym_fixed_ancGrowth(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,p_misid)

    Isolation-with-migration model with exponential pop growth before and after split.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m: symmetric migration rate for each pop (2*Na*m)
    Tg: Time of growth in ancestral population minus T
    nua: Final size of ancestral population
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,m,Tg,nua = params
    p_misid=0.001
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    nua_func = lambda t: numpy.exp(numpy.log(nua) * t/(Tg+T))
    phi = dadi.Integration.one_pop(phi, xx, Tg, nua_func)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def readLFromFSFile(fsFileName):
    with open(fsFileName) as fsFile:
        lines = []
        for line in fsFile:
            if not line.strip().startswith("#"):
                lines.append(line)
        ns = [int(x)-1 for x in lines[0].split()[:2]]
        L = sum([int(x) for x in lines[1].split()])
    return L, ns

inFileName = sys.argv[1]
swarmSize = int(sys.argv[2])
gensPerYear = float(sys.argv[3])

L, ns = readLFromFSFile(inFileName)
data = dadi.Spectrum.from_file(inFileName)

firstSize=120
pts_l = [firstSize,firstSize+10,firstSize+20]

func = IM_sym_fixed_ancGrowth

upper_bound = [1e2, 1e2, 1e2, 1e2, 2, 20, 2, 1e2]
lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 0, 0, 0, 1e-2]

p1=dadiFunctions.makeRandomParams(lower_bound,upper_bound)

func_ex = dadi.Numerics.make_extrap_func(func)

# Instantiate Optimization Problem 

def objfunc(x):
    f = dadi.Inference._object_func(x, data, func_ex, pts_l, 
                                      lower_bound=lower_bound,
                                          upper_bound=upper_bound)
    g=[]
    fail = 0
    return f,g,fail
    
opt_prob = pyOpt.Optimization('dadi optimization',objfunc)
opt_prob.addVar('nu1_0','c',lower=lower_bound[0],upper=upper_bound[0],value=p1[0])
opt_prob.addVar('nu2_0','c',lower=lower_bound[1],upper=upper_bound[1],value=p1[1])
opt_prob.addVar('nu1','c',lower=lower_bound[2],upper=upper_bound[2],value=p1[2])
opt_prob.addVar('nu2','c',lower=lower_bound[3],upper=upper_bound[3],value=p1[3])
opt_prob.addVar('T','c',lower=lower_bound[4],upper=upper_bound[4],value=p1[4])
opt_prob.addVar('m','c',lower=lower_bound[5],upper=upper_bound[5],value=p1[5])
opt_prob.addVar('Tg','c',lower=lower_bound[6],upper=upper_bound[6],value=p1[6])
opt_prob.addVar('nua','c',lower=lower_bound[7],upper=upper_bound[7],value=p1[7])
opt_prob.addObj('f')

if myrank == 0:
    print opt_prob

#optimize
psqp = pyOpt.ALPSO(pll_type='DPM')
psqp.setOption('printOuterIters',1)
#psqp.setOption('maxOuterIter',1)
#psqp.setOption('stopCriteria',0)
psqp.setOption('SwarmSize',swarmSize)
psqp(opt_prob)
print opt_prob.solution(0)

popt = numpy.zeros(len(p1))
for i in opt_prob._solutions[0]._variables:
    popt[i]= opt_prob._solutions[0]._variables[i].__dict__['value']

model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, data)
mu=3.5e-9
if myrank == 0:
    print 'Optimized log-likelihood:', ll_opt
    print 'AIC:', -2*ll_opt + 2*len(popt)
    #scaled estimates
    theta0 = dadi.Inference.optimal_sfs_scaling(model, data)
    Nref= theta0 / mu / L / 4

    print 'Nref:',Nref
    paramsTxt =['nu1_0','nu2_0','nu1','nu2','T', '2Nref_m', "Tg", "nua"]
    scaledParams = [Nref*popt[0],Nref*popt[1],Nref*popt[2],Nref*popt[3],2*Nref/gensPerYear*popt[4], popt[5], 2*Nref/gensPerYear*popt[4]+2*Nref/gensPerYear*popt[6], Nref*popt[7]]
    for i in range(len(paramsTxt)):
        print paramsTxt[i],':',str(scaledParams[i])
    print ""
    print repr(popt)

############### 
# Now refine the optimization using Local Optimizer
# Instantiate Optimizer (SLSQP) 
# Instantiate Optimizer (SLSQP)
slsqp = pyOpt.SLSQP()
# Solve Problem (With Parallel Gradient)
if myrank == 0:
    print 'going for second optimization'

slsqp(opt_prob.solution(0),sens_type='FD',sens_mode='pgc')
print opt_prob.solution(0).solution(0)
opt = numpy.zeros(len(p1))
for i in opt_prob._solutions[0]._solutions[0]._variables:
    popt[i]= opt_prob._solutions[0]._solutions[0]._variables[i].__dict__['value']

model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, data)
if myrank == 0:      
    print 'After Second Optimization'
    print 'Optimized log-likelihood:', ll_opt
    print 'AIC:', -2*ll_opt + 2*len(popt)

    #scaled estimates
    theta0 = dadi.Inference.optimal_sfs_scaling(model, data)
    print 'with u = %e' %(mu)
    Nref= theta0 / mu / L / 4

    print 'Nref:',Nref
    paramsTxt =['nu1_0','nu2_0','nu1','nu2','T', '2Nref_m', "Tg", "nua"]
    scaledParams = [Nref*popt[0],Nref*popt[1],Nref*popt[2],Nref*popt[3],2*Nref/gensPerYear*popt[4], popt[5], 2*Nref/gensPerYear*popt[4]+2*Nref/gensPerYear*popt[6], Nref*popt[7]]
    for i in range(len(paramsTxt)):
        print paramsTxt[i],':',str(scaledParams[i])
    print ""
    print repr(popt)
