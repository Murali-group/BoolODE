def genDropouts(opts):
    """
    :param outPrefix: default='', Prefix for output files.
    :type outPrefix:  str
    :param expr: Path to expression data file
    :type expr: str
    :param pseudo: Path to pseudotime file
    :type pseudo: str
    :param refNet: Path to reference network file
    :type refNet: str
    :param nCells: default=100, Number of cells to sample.    
    :type nCells: int
    :param dropout: default=False, Carry out dropout analysis?
    :type dropout: bool
    :param drop-cutoff: default=0.5, Specify quantile cutoff on gene expression    
    :type drop-cutoff: float
    :param drop-prob: default=0.5, Specify the probability of dropping a gene below quantile q     
    :type drop-prob: float
    :param samplenum:, default=None, Sample Number
    :type samplenum: int
    """
    
    if opts.samplenum is None:
        print('Please specify sample number')
        sys.exit()
    if opts.dropout:
        dropoutCutoffs = opts.drop_cutoff
    else:
        dropoutCutoffs = 0

    ## Read the ExpressionData.csv file
    expDF = pd.read_csv(opts.expr, index_col=0)
    ## Read the PseudoTime.csv file
    PTDF = pd.read_csv(opts.pseudo, index_col=0)
    ## Read the refNetwork.csv file
    refDF = pd.read_csv(opts.refNet, index_col=0)

    ## Generate output path for the dropout datasets
    path = opts.outPrefix + '-' +str(opts.nCells) +'-' + str(opts.samplenum)
    if dropoutCutoffs != 0:
        path += '-' + str(int(100*dropoutCutoffs))+  '-' + str( str(opts.drop_prob))
    if not os.path.exists(path):
        os.makedirs(path)
    # Dropout here
    # copy over PT and refNetwork files
    refDF.to_csv(path + '/refNetwork.csv')
    PTDF.to_csv(path+'/PseudoTime.csv')        
    # Drop-out genes if they are less than the 
    # percentile value @ "dc" with 50% chance
    DropOutDF = expDF.copy()
    if dropoutCutoffs != 0:
        quantileExp = expDF.quantile(q = dropoutCutoffs, axis = 'columns')
        for idx, row in tqdm(expDF.iterrows()):
            for col in expDF.columns:
                if row[col] < quantileExp.loc[idx]:
                    cointoss = np.random.random()
                    if cointoss < opts.drop_prob:
                        DropOutDF.loc[idx,col] = 0.0

    DropOutDF.to_csv(path + '/ExpressionData.csv')
