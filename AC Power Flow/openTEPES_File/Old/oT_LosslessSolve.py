# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.0.0 - April 26, 2020

# solve lossless first
if pIndNetLosses:

    # deactivate loss constraints
    mTEPES.eLineLosses1.deactivate()
    mTEPES.eLineLosses2.deactivate()

    #%% solving the problem
    mTEPES.write('openTEPES_LL_'+CaseName+'.lp', io_options={'symbolic_solver_labels':True})  # create lp-format file
    Solver = SolverFactory('gurobi')  # select solver
    Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
    # Solver.options['ResultFile'  ] = 'openTEPES_'+CaseName+'.ilp'                           # introduced to show results of IIS
    Solver.options['Method'        ] = 2  # barrier method
    Solver.options['MIPGap'        ] = 0.03
    Solver.options['Threads'       ] = round((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)
    # Solver.options['TimeLimit'     ] = 3600
    # Solver.options['IterationLimit'] = 3600000

    SolvingTime = time.time() - StartTime
    StartTime = time.time()
    print('Solving lossless                      status ... ', round(SolvingTime), 's')

    # activate loss constraints
    mTEPES.eLineLosses1.activate()
    mTEPES.eLineLosses2.activate()
