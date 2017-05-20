#A collcetion of functions for dealing with Dadi models
# A. Kern


import dadi
import numpy
import scipy
import pylab
#import nlopt



######### Demographic stuff

def three_epoch_fixed_duration(params, ns, pts):
    """
    params = (nuB,nuF,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TF: Time since bottleneck recovery (in units of 2*Na generations) 

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TF = params
    TB=500/(2.0*10000) #TB is the duration of the bottleneck

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def OutOfAfricaGrowB((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
                  TAf, TB, TEuNA, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA))
    
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
                               m12=0, m21=0)
    nuEu0 = nuEu_func(TB)
   
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu_func, nu3=nuNA_func, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)

    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)


def OutOfAfrica_admix((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
                  TAf, TB, TEuNA,T_ad, p_ad, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA+T_ad))

    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
                               m12=0, m21=0)
    nuEu0 = nuEu_func(TB)

    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA+T_ad))
    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA+T_ad)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu_func, nu3=nuNA_func, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)
    nuEu0 = nuEu_func(TEuNA)
    nuNA0 = nuNA_func(TEuNA)
    phi = dadi.PhiManip.phi_3D_admix_1_and_2_into_3(phi, p_ad,0, xx,xx,xx)
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/T_ad)
    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/T_ad)
    phi = dadi.Integration.three_pops(phi, xx, T_ad, nu1=nuAf, 
                                 nu2=nuEu_func, nu3=nuNA_func, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)
    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)


def OutOfAfrica_mig_Af_NA((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
                  TAf, TB, TEuNA,mNA_Af,mAf_NA, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)

	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA))

	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
	                               m12=0, m21=0)
	    nuEu0 = nuEu_func(TB)

	    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=0, m13=mAf_NA, m21=0, m23=0,
	                                 m31=mNA_Af, m32=0)

	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)
	
def OutOfAfrica_mig((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
	                  TAf, TB, TEuNA,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)

	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA))

	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
	                               m12=mAf_Eu, m21=mEu_Af)
	    nuEu0 = nuEu_func(TB)

	    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)

	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica_mig_noAncient((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
	                  TAf, TB, TEuNA,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)

	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA))

	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
	                               m12=0, m21=0)
	    nuEu0 = nuEu_func(TB)

	    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)

	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)


def OutOfAfrica_mig_admix((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
	                  TAf, TB, TEuNA,T_ad,p_ad,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)

	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA+T_ad))

	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
	                               m12=mAf_Eu, m21=mEu_Af)
	    nuEu0 = nuEu_func(TB)

	    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA+T_ad))
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA+T_ad)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)
	    nuEu0 = nuEu_func(TEuNA)
	    nuNA0 = nuNA_func(TEuNA)
	    phi = dadi.PhiManip.phi_3D_admix_1_and_2_into_3(phi, p_ad,0, xx,xx,xx)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/T_ad)
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/T_ad)
	    phi = dadi.Integration.three_pops(phi, xx, T_ad, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)
	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica2((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
                  TAf, TB, TEuNA, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)


    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu0, 
                               m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu_func, nu3=nuNA_func, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)
    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica2_mig((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
                  TAf, TB, TEuNA,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)


    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu0, 
                               m12=mAf_Eu, m21=mEu_Af)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA))
    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu_func, nu3=nuNA_func, 
                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
                                 m31=mNA_Af, m32=mNA_Eu)
    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica2_mig_admix((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
	                  TAf, TB, TEuNA,T_ad,p_ad,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)
	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)
	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)


	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu0, 
	                               m12=mAf_Eu, m21=mEu_Af)
	    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TEuNA+T_ad))
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA+T_ad)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)
	    nuEu0 = nuEu_func(TEuNA)
	    nuNA0 = nuNA_func(TEuNA)
	    phi = dadi.PhiManip.phi_3D_admix_1_and_2_into_3(phi, p_ad,0, xx,xx,xx)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/T_ad)
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/T_ad)
	    phi = dadi.Integration.three_pops(phi, xx, T_ad, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)   
	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica3((nuAf, nuEu, nuNA, 
                  TAf, TB, TEuNA, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)


    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu, 
                               m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu, nu3=nuNA, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)
    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)


def OutOfAfrica4((nuAf, nuEu, nuNA, 
                  TAf, TB, TEuNA, p_misid), (n1,n2,n3), pts):
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)


    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu, 
                               m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
                                 nu2=nuEu, nu3=nuNA, 
                                 m12=0, m13=0, m21=0, m23=0,
                                 m31=0, m32=0)
    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def OutOfAfrica_mig_admix2((nuAf, nuEu0, nuEu, nuNA0, nuNA, 
	                  TAf, TB, TEuNA,T_ad,mAf_Eu,mAf_NA,mEu_Af,mEu_NA,mNA_Af,mNA_Eu, p_misid), (n1,n2,n3), pts):
	    xx = dadi.Numerics.default_grid(pts)

	    phi = dadi.PhiManip.phi_1D(xx)
	    phi = dadi.Integration.one_pop(phi, xx, TAf, nu=nuAf)

	    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/(TB+TEuNA+T_ad))

	    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2=nuEu_func, 
	                               m12=mAf_Eu, m21=mEu_Af)
	    nuEu0 = nuEu_func(TB)

	    phi = dadi.PhiManip.phi_2D_to_3D_admix(phi,p_ad,xx,xx,xx)
	    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/TEuNA)
	    nuNA_func = lambda t: nuNA0*(nuNA/nuNA0)**(t/TEuNA)
	    phi = dadi.Integration.three_pops(phi, xx, TEuNA, nu1=nuAf, 
	                                 nu2=nuEu_func, nu3=nuNA_func, 
	                                 m12=mAf_Eu, m13=mAf_NA, m21=mEu_Af, m23=mEu_NA,
	                                 m31=mNA_Af, m32=mNA_Eu)

	    fs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
	    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)
	
###################################################
################## Two Populations 

##two population model with misorientation
def IM_misorient_5epoch(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu1_1,nu1_2,nu1_3,nu1_4,nu2_0,nu2_1,nu2_2,nu2_3,nu2_4,t0,t1,t2,t3,t4,m12,m21,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu1_1,nu1_2,nu1_3,nu1_4,nu2_0,nu2_1,nu2_2,nu2_3,nu2_4,t0,t1,t2,t3,t4,m12,m21,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, t0, nu1_0, nu2_0,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, t1, nu1_1, nu2_1,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, t2, nu1_2, nu2_2,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, t3, nu1_3, nu2_3,
                               m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, t4, nu1_4, nu2_4,
                               m12=m12, m21=m21) 
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)	
	
##two population model with misorientation
def IM_misorient(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,m12,m21,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

def IM_const(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,T,m12,m21)

    Isolation-with-migration model.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu1_0,nu2_0,T,m12,m21 = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_0, nu2_0, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM_noMig_const(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,T,m12,m21)

    Isolation model.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu1_0,nu2_0,T = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_0, nu2_0, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split.
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu1_0,nu2_0,nu1,nu2,T,m12,m21 = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM_ancGrowth(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T)

    Isolation-with-migration model with exponential pop growth before and after split

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: migration rate for pop 2 to pop 1 (2*Na*m)
    m21: migration rate for pop 1 to pop 2 (2*Na*m)
    Tg: Time of growth in ancestral population minus T
    nua: Final size of ancestral population
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu1_0,nu2_0,nu1,nu2,T,m12,m21,Tg,nua = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    nua_func = lambda t: numpy.exp(numpy.log(nua) * t/(Tg+T))
    phi = dadi.Integration.one_pop(phi, xx, Tg, nua_func)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM_ancGrowth_const1(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T)

    Isolation-with-migration model with exponential pop growth before split and after in pop 2

    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: migration rate for pop 2 to pop 1 (2*Na*m)
    m21: migration rate for pop 1 to pop 2 (2*Na*m)
    Tg: Time of growth in ancestral population minus T
    nua: Final size of ancestral population
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu2_0,nu1,nu2,T,m12,m21,Tg,nua = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    nua_func = lambda t: numpy.exp(numpy.log(nua) * t/(Tg+T))
    phi = dadi.Integration.one_pop(phi, xx, Tg, nua_func)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM_noMig(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split.
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    try:
        nu1_0,nu2_0,nu1,nu2,T = params
    except ValueError:
        print "params:", params, len(params)
        raise ValueError

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def contraction_and_growth(params, ns, pts):
    """
    params = (nuB,nuG,TB,TG)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuG: Ratio of contemporary to bottleneck pop size
    TB: Length of bottleneck (in units of 2*Na generations) 
    TG: Time of growth onset (in units of 2*Na generations) 

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuG,TB,TG = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuG/nuB) * t/TG)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TG, nu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

##two population model with misorientation
def IM_misorient_noMig(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

##two population model with misorientation
def IM_misorient_admix(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad,p_ad,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad,p_ad,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T+t_ad))
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/(T+t_ad))
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad, xx,xx)
    nu1_0 = nu1_func(t_ad)
    nu2_0 = nu2_func(t_ad)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad)
    phi = dadi.Integration.two_pops(phi, xx, t_ad, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)	

##two population model with misorientation
def IM_misorient_doubleAdmix(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad1,p_ad1,t_ad2,p_ad2,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad1,p_ad1,t_ad2,p_ad2,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T+t_ad1+t_ad2))
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/(T+t_ad1+t_ad2))
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad1, xx,xx)
    nu1_0 = nu1_func(t_ad1)
    nu2_0 = nu2_func(t_ad1)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad1+t_ad2)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad1+t_ad2)
    phi = dadi.Integration.two_pops(phi, xx, t_ad1, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad2, xx,xx)
    nu1_0 = nu1_func(t_ad2)
    nu2_0 = nu2_func(t_ad2)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad2)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad2)
    phi = dadi.Integration.two_pops(phi, xx, t_ad2, nu1_func, nu2_func,
                               m12=m12, m21=m21)


    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

##two population model with misorientation
def IM_misorient_doubleAdmix_noMig(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad1,p_ad1,t_ad2,p_ad2,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 

    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,t_ad1,p_ad1,t_ad2,p_ad2,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T+t_ad1+t_ad2))
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/(T+t_ad1+t_ad2))
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=0, m21=0)

    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad1, xx,xx)
    nu1_0 = nu1_func(t_ad1)
    nu2_0 = nu2_func(t_ad1)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad1+t_ad2)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad1+t_ad2)
    phi = dadi.Integration.two_pops(phi, xx, t_ad1, nu1_func, nu2_func,
                               m12=0, m21=0)
    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad2, xx,xx)
    nu1_0 = nu1_func(t_ad2)
    nu2_0 = nu2_func(t_ad2)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad2)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad2)
    phi = dadi.Integration.two_pops(phi, xx, t_ad2, nu1_func, nu2_func,
                               m12=0, m21=0)


    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)
	
##two population model with misorientation
def IM_misorient_noMig_admix(params, ns, pts):
    """
    ns = (n1,n2)
    params = (nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad,p_ad,p_misid)

    Isolation-with-migration model with exponential pop growth.

    nu1_0: Size of pop 1 after split.
    nu2_0: Size of pop 2 after split. 
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 

    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1_0,nu2_0,nu1,nu2,T,t_ad,p_ad,p_misid = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T+t_ad))
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/(T+t_ad))
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=0, m21=0)

    phi = dadi.PhiManip.phi_2D_admix_1_into_2(phi, p_ad, xx,xx)
    nu1_0 = nu1_func(t_ad)
    nu2_0 = nu2_func(t_ad)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/t_ad)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/t_ad)
    phi = dadi.Integration.two_pops(phi, xx, t_ad, nu1_func, nu2_func,
                               m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return (1-p_misid)*fs + p_misid * dadi.Numerics.reverse_array(fs)

##########################
#######
####  Helper functions

def makeRandomParams(lower,upper):
    pNew=numpy.zeros(len(lower))
    for i in range(len(lower)):
        pNew[i]= numpy.random.uniform(lower[i],upper[i])
    return pNew


def plot2file_3d_comp_multinom(model, data, filename,vmin=None, vmax=None,
                          resid_range=None, fig_num=None,
                          pop_ids=None, residual='Anscombe', adjust=True):
    """
    Multinomial comparison between 3d model and data.


    model: 3-dimensional model SFS
    data: 3-dimensional data SFS
    vmin, vmax: Minimum and maximum values plotted for sfs are vmin and
                vmax respectively.
    resid_range: Residual plot saturates at +- resid_range.
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    adjust: Should method use automatic 'subplots_adjust'? For advanced
            manipulation of plots, it may be useful to make this False.

    This comparison is multinomial in that it rescales the model to optimally
    fit the data.
    """
    model = dadi.Inference.optimally_scaled_sfs(model, data)

    plot2file_3d_comp_Poisson(model, data,filename, vmin=vmin, vmax=vmax,
                         resid_range=resid_range, fig_num=fig_num,
                         pop_ids=pop_ids, residual=residual,
                         adjust=adjust)

def plot2file_3d_comp_Poisson(model, data,filename, vmin=None, vmax=None,
                         resid_range=None, fig_num=None, pop_ids=None, 
                         residual='Anscombe', adjust=True):
    """
    Poisson comparison between 3d model and data.


    model: 3-dimensional model SFS
    data: 3-dimensional data SFS
    vmin, vmax: Minimum and maximum values plotted for sfs are vmin and
                vmax respectively.
    resid_range: Residual plot saturates at +- resid_range.
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    adjust: Should method use automatic 'subplots_adjust'? For advanced
            manipulation of plots, it may be useful to make this False.
    """
    if data.folded and not model.folded:
        model = model.fold()

    masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)

    if fig_num is None:
        f = pylab.gcf()
    else:
        f = pylab.figure(fig_num, figsize=(8,10))

    pylab.clf()
    if adjust:
        pylab.subplots_adjust(bottom=0.07, left=0.07, top=0.95, right=0.95)

    modelmax = max(masked_model.sum(axis=sax).max() for sax in range(3))
    datamax = max(masked_data.sum(axis=sax).max() for sax in range(3))
    modelmin = min(masked_model.sum(axis=sax).min() for sax in range(3))
    datamin = min(masked_data.sum(axis=sax).min() for sax in range(3))
    max_toplot = max(modelmax, datamax)
    min_toplot = min(modelmin, datamin)

    if vmax is None:
        vmax = max_toplot
    if vmin is None:
        vmin = min_toplot
    extend = dadi.Plotting._extend_mapping[vmin <= min_toplot, vmax >= max_toplot]

    # Calculate the residuals
    if residual == 'Anscombe':
        resids = [dadi.Inference.\
                  Anscombe_Poisson_residual(masked_model.sum(axis=2-sax), 
                                            masked_data.sum(axis=2-sax), 
                                            mask=vmin) for sax in range(3)]
    elif residual == 'linear':
        resids =[dadi.Inference.\
                 linear_Poisson_residual(masked_model.sum(axis=2-sax), 
                                         masked_data.sum(axis=2-sax), 
                                         mask=vmin) for sax in range(3)]
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)


    min_resid = min([r.min() for r in resids])
    max_resid = max([r.max() for r in resids])
    if resid_range is None:
        resid_range = max((abs(max_resid), abs(min_resid)))
    resid_extend = dadi.Plotting._extend_mapping[-resid_range <= min_resid, 
                                   resid_range >= max_resid]

    if pop_ids is not None:
        if len(pop_ids) != 3:
            raise ValueError('pop_ids must be of length 3.')
        data_ids = model_ids = resid_ids = pop_ids
    else:
        data_ids = masked_data.pop_ids
        model_ids = masked_model.pop_ids

        if model_ids is None:
            model_ids = data_ids

        if model_ids == data_ids:
           resid_ids = model_ids
        else:
            resid_ids = None

    for sax in range(3):
        marg_data = masked_data.sum(axis=2-sax)
        marg_model = masked_model.sum(axis=2-sax)

        curr_ids = []
        for ids in [data_ids, model_ids, resid_ids]:
            if ids is None:
                ids = ['pop0', 'pop1', 'pop2']

            if ids is not None:
                ids = list(ids)
                del ids[2-sax]

            curr_ids.append(ids)

        ax = pylab.subplot(4,3,sax+1)
        plot_colorbar = (sax == 2)
        dadi.Plotting.plot_single_2d_sfs(marg_data, vmin=vmin, vmax=vmax, pop_ids=curr_ids[0],
                           extend=extend, colorbar=plot_colorbar)

        pylab.subplot(4,3,sax+4, sharex=ax, sharey=ax)
        dadi.Plotting.plot_single_2d_sfs(marg_model, vmin=vmin, vmax=vmax, 
                           pop_ids=curr_ids[1], extend=extend, colorbar=False)

        resid = resids[sax]
        pylab.subplot(4,3,sax+7, sharex=ax, sharey=ax)
        dadi.Plotting.plot_2d_resid(resid, resid_range, pop_ids=curr_ids[2],
                      extend=resid_extend, colorbar=plot_colorbar)

        ax = pylab.subplot(4,3,sax+10)
        flatresid = numpy.compress(numpy.logical_not(resid.mask.ravel()), 
                                   resid.ravel())
        ax.hist(flatresid, bins=20, normed=True)
        ax.set_yticks([])
    pylab.savefig(filename, bbox_inches='tight')

################################################
## MS stuff
## and discoal... and msAdmix....

##########



def IM_misorient_admix_core(params):
    """
    msAdmix core command for IM_misorient_admix.
    """
    nu1_0,nu2_0,nu1,nu2,T,m12,m21,t_ad,p_ad,p_misid = params

    alpha1 = numpy.log(nu1/nu1_0)/T
    alpha2 = numpy.log(nu2/nu2_0)/T
    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-eg 0 1 %(alpha1)f -eg 0 2 %(alpha2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-eA %(t_ad)f 2 1 %(p_ad)f "\
            "-ej %(T)f 2 1 -en %(T)f 1 1"

    sub_dict = {'nu1':nu1, 'nu2':nu2, 'alpha1':2*alpha1, 'alpha2':2*alpha2,
                'm12':2*m12, 'm21':2*m21, 'T': T/2, 't_ad':t_ad/2, 'p_ad':p_ad}

    return command % sub_dict

def msAdmix_command(theta, ns, core, iter, recomb=0, rsites=None):
    """
    Generate ms command for simulation from core.

    theta: Assumed theta
    ns: Sample sizes
    core: Core of ms command that specifies demography.
    iter: Iterations to run ms
    recomb: Assumed recombination rate
    rsites: Sites for recombination. If None, default is 10*theta.
    """
    if len(ns) > 1:
        ms_command = "msAdmix %(total_chrom)i %(iter)i -t %(theta)f -I %(numpops)i "\
                "%(sample_sizes)s %(core)s"
    else:
        ms_command = "msAdmix %(total_chrom)i %(iter)i -t %(theta)f  %(core)s"

    if recomb:
        ms_command = ms_command + " -r %(recomb)f %(rsites)i"
        if not rsites:
            rsites = theta*10
    sub_dict = {'total_chrom': numpy.sum(ns), 'iter': iter, 'theta': theta,
                'numpops': len(ns), 'sample_sizes': ' '.join(map(str, ns)),
                'core': core, 'recomb': recomb, 'rsites': rsites}

    return ms_command % sub_dict

