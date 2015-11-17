import numpy as np


# dict CONFIG will be imported by pylitho.py

CONFIG = {
    'MPI' : True,
    'MPI_VARIATION_TYPE' : 20,
    'MPI_VARIATION_PARAMS' : np.array((100.,300.,5e3,50e3,50,300)),
    'MPI_VIRT_PROCS' : 10*10*10,


    'OUTPUT_FILE' : True,
    'OUTPUT_FILE_EVERY_TSTEP' : 1000,
    'OUTDIR' : "output",  # no trailing slash!
    'MODELNAME' : "001-ini",
    'OUTPUT_OVERWRITE' : True,

    'RESTART' : False,
    'RESTART_INDIR' : "output",
    'RESTART_MODELNAME' : "test1",
    'RESTART_TSTEP' : 9450,
    'RESTART_POST_MOD' : 0,   # post-restart modifications

    'NX' : 200,
    'L_KM' : (0, 100), # in km
    'MAXTIME_MA' : 500, # in Ma
    'MAXRUNTIME_MA' : 500, # in Ma
    'TSTEP_MULTI' : 0.5,

    'DIFFSCHANGE_ACCURACY' : 1e-10,  # max diffus. change, iteration criteria

    'MOHO_DEPTH_KM' : 35,

    # ****
    # define type of initial T field
    # if restart==True, use type -1 to prevent overwriting
    # ****
    #'TINI_TYPE' : -1,
    'TINI_TYPE' : 0,
    #'TINI_TYPE' : 1,
    #'TINI_TYPE' : 10,

    # ****
    # define type of k=k(T) relation
    # ****

    # comment/uncomment as needed

    'KT_RELATION_TYPE' : 0,
    'KT_RELATION_PARAMS' : np.array(()),

    #'KT_RELATION_TYPE' : 1,
    #'KT_RELATION_PARAMS' : np.array((1e-3),ndmin=1),

    #'KT_RELATION_TYPE' : 2,
    #'KT_RELATION_PARAMS' : np.array((1e-3, 5e-10)),


    # ****
    # define type of cp=cp(T) relation
    # ****

    'CT_RELATION_TYPE' : 0,
    'CT_RELATION_PARAMS' : np.array(()),

    #'CT_RELATION_TYPE' : 1,
    #'CT_RELATION_PARAMS' : np.array((8.95e-10, -2.13e-6, 0.00172, 0.716, 750.0)),


    # ****
    # define type of k0 field
    # ****
    'K0_TYPE' : 10,

    # ****
    # define type of cp0 field
    # ****
    #'C0_TYPE' : 0,
    'C0_TYPE' : 10,

    # ****
    # define type of rho field
    # ****
    #'RHO0_TYPE' : 0,
    'RHO0_TYPE' : 10,


    # ****
    # define type of H field
    # ****
    #'H0_TYPE' : 0,
    'H0_TYPE' : 10,
    #'H0_TYPE' : 2,

    'BND_BOT_TYPE' : 0,
    'BND_BOT_HFLOW' : 0.0,
    'BND_BOT_TEMP' : 1350.0,
    'BND_TOP_TEMP' : 0,
    'EROSION_SPEED_M_MA' : 0,
    'EROSION_SPEED_TYPE' : 1,
    'MAX_TIME_TO_ERODE_MA' : 0,
}
