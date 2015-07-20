import numpy as np


# dict CONFIG will be imported by pylitho.py

CONFIG = {
    'MPI' : True,
    'MPI_VARIATION_TYPE' : 1,
    'MPI_VARIATION_PARAMS' : np.array((5e-4,15e-4)),
    'MPI_VIRT_PROCS' : 20,


    'OUTPUT_FILE' : True,
    'OUTPUT_FILE_EVERY_TSTEP' : 100,
    'OUTDIR' : "output",  # no trailing slash!
    'MODELNAME' : "mpitest1",
    'OUTPUT_OVERWRITE' : True,

    'RESTART' : False,
    'RESTART_INDIR' : "output",
    'RESTART_MODELNAME' : "test1",
    'RESTART_TSTEP' : 9450,
    'RESTART_POST_MOD' : 0,   # post-restart modifications

    'NX' : 100,
    'L_KM' : (30, 130), # in km
    'MAXTIME_MA' : 30, # in Ma
    'MAXRUNTIME_MA' : 30, # in Ma
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

    #'KT_RELATION_TYPE' : 0,
    #'KT_RELATION_PARAMS' : np.array(()),

    'KT_RELATION_TYPE' : 1,
    'KT_RELATION_PARAMS' : np.array((1e-3),ndmin=1),

    #'KT_RELATION_TYPE' : 2,
    #'KT_RELATION_PARAMS' : np.array((1e-3, 5e-10)),


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

    'BND_BOT_TYPE' : 1,
    'BND_BOT_HFLOW' : 35e-3,
    'BND_TOP_TEMP' : 0,
    'EROSION_SPEED_M_MA' : 0,
    'EROSION_SPEED_TYPE' : 1,
}
