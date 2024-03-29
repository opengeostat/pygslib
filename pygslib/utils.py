import pygslib
import numpy as np

def gcos(table, vardhole, varmodel, cutoff, model, cbb, decwt = None,
        ltail=1, utail=4, ltpar=1, utpar=1.5, K=30,
        ymin=None, ymax=None, ndisc = 1000):
    """
    gcos(table, mask, vardhole, varmodel, cutoff, model, cbb, decwt = None)
    
    Produce global change of support validation (GCOS)
    
    Inputs
    ------
    table: pandas dataframe
        drillhole database
    vardhole: str
        variable on drillhole
    varmodel:
        variable in bmodel
    cutoff:
        cutoff list
    model: pandas dataframe
        block model 
    cbb:
        block variance
    decwt: str default None
        declustering weight variable name
    ltail, utail, ltpar, utpar: floats
        lower and upper tail model. The default parameters are ltail=1, utail=4, ltpar=1, utpar=1.5
    K: integer
        number of polynomials. Default is 30
    ymin, ymax: 
        gaussian values validity interval (default None)
    ndisc: int 
        discretization intervals in cdf. Default is 1000
    """
    
    # firt the discrete gausian model 
    zin = table[vardhole]
    
    if decwt is None:
        w = np.ones(table[vardhole].shape[0])
    else:
        w = table[decwt]
        
    zmin = table[vardhole].min(),
    zmax = table[vardhole].max(),  
    
    PCI, H, raw, zana, gauss, z, P, raw_var, PCI_var, fig1, \
    _, _, _, _, _, _, _, _ = pygslib.nonlinear.anamor(
       zin, w, ltail=ltail, utail=utail, ltpar=ltpar, utpar=utpar, K=K,
       ymin=ymin, ymax=ymax, ndisc = ndisc)


    r = pygslib.nonlinear.get_r(Var_Zv = cbb, PCI = PCI)
    
    ZV, PV, fig2, \
    _, _, _, _, _, _, _, _  = pygslib.nonlinear.anamor_blk( PCI, H, r = r, gauss = gauss, Z = z,
                  ltail=ltail, utail=utail, ltpar=ltpar, utpar=utpar,
                  raw=raw, zana=zana)
    
    # calculate GTC
    tt = []
    gg = []
    label = []

    # calculate GTC from gaussian in block support
    
    cutoff[cutoff>ZV.max()] = ZV.max() - 0.000001 # there is an error of not like this
    cutoff[cutoff<ZV.min()] = ZV.min() + 0.000001
    t,ga,gb = pygslib.nonlinear.gtcurve (cutoff = cutoff, z=ZV, p=PV)
    
    
    
    tt.append(t)
    gg.append(ga)
    label.append(f'DGM {vardhole}')

    fig = pygslib.nonlinear.plotgt(cutoff = cutoff, t = tt, g = gg, label = label)
    
    
    
    # to compare global resources with the one estimated we calculate the CDF of the blocks

    # cdf of kriging estimate
    for v in varmodel:
        parameters_probplt = {
                'iwt'  : 0,                            
                'va'   : model[v][model[v].notnull()].values,    
                'wt'   : np.ones(model[v][model[v].notnull()].shape[0])} 

    
        binval_ok,cl_ok,xpt025,xlqt,xmed,xuqt,xpt975,xmin,xmax, \
                xcvr,xmen,xvar,error = pygslib.gslib.__plot.probplt(**parameters_probplt)
    

        # calculate GTC ok
        cutoff[cutoff>cl_ok.max()] = cl_ok.max() - 0.000001 # there is an error of not like this
        cutoff[cutoff<cl_ok.min()] = cl_ok.min() + 0.000001
        t,ga,gb = pygslib.nonlinear.gtcurve (cutoff = cutoff, z=cl_ok,
                       p=binval_ok, varred = 1, ivtyp = 2, zmin = 0, zmax = None,
                      ltail = ltail, ltpar = ltpar, middle = 1, mpar = 1, utail = utail,
                      utpar =  utpar, maxdis = ndisc)
        
        
        tt.append(t)
        gg.append(ga)
        label.append(f'Estimate {v}')
        
    return pygslib.nonlinear.plotgt(cutoff = cutoff, t = tt, g = gg, label = label) 
