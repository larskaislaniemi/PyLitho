import numpy as np


# dict CONFIG will be imported by pylitho.py

CONFIG = {
    'OUTPUT_FILE' : True,
    'OUTPUT_FILE_EVERY_TSTEP' : 1000,
    'OUTDIR' : "output",  # no trailing slash!
    'MODELNAME' : "test1-re",
    'OUTPUT_OVERWRITE' : True,

    'RESTART' : True,
    'RESTART_INDIR' : "output",
    'RESTART_MODELNAME' : "test1",
    'RESTART_TSTEP' : 252000,
    'RESTART_POST_MOD' : 1,   # post-restart modifications

    'NX' : 130,
    'L_KM' : (15, 115), # in km
    'MAXTIME_MA' : 1001, # in Ma
    'MAXRUNTIME_MA' : 1, # in Ma
    'TSTEP_MULTI' : 0.5,

    'DIFFSCHANGE_ACCURACY' : 1e-10,  # max diffus. change, iteration criteria

    'MOHO_DEPTH_KM' : 35,

    # ****
    # define type of initial T field
    # if restart==True, use type -1 to prevent overwriting
    # ****
    'TINI_TYPE' : -1,
    #'TINI_TYPE' : 0,
    #'TINI_TYPE' : 1,
    #'TINI_TYPE' : 10,

    # ****
    # define type of k=k(T) relation
    # ****

    # comment/uncomment as needed

    #'KT_RELATION_TYPE' : 0,
    #'KT_RELATION_PARAMS' : np.array(()),

    #'KT_RELATION_TYPE' : 1,
    #'KT_RELATION_PARAMS' : np.array((1e-3),ndmin=1),

    'KT_RELATION_TYPE' : 2,
    'KT_RELATION_PARAMS' : np.array((1e-3, 5e-10)),


    # ****
    # define type of cp=cp(T) relation
    # ****

    #'CT_RELATION_TYPE' : 0,
    #'CT_RELATION_PARAMS' : np.array(()),

    'CT_RELATION_TYPE' : 1,
    'CT_RELATION_PARAMS' : np.array((8.95e-10, -2.13e-6, 0.00172, 0.716, 750.0)),


    # ****
    # define type of k0 field
    # ****
    'K0_TYPE' : 0,

    # ****
    # define type of cp0 field
    # ****
    #'C0_TYPE' : 0,
    'C0_TYPE' : 1,

    # ****
    # define type of rho field
    # ****
    #'RHO0_TYPE' : 0,
    'RHO0_TYPE' : 1,


    # ****
    # define type of H field
    # ****
    #'H0_TYPE' : 0,
    'H0_TYPE' : 1,
    #'H0_TYPE' : 2,


    'BND_BOT_HFLOW' : 35e-3, #35e-3,
    'BND_TOP_TEMP' : 0,
    'EROSION_SPEED_M_MA' : 0,
}
