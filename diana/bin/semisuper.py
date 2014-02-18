


class Experiment(object):

    @profile
    def __init__(self, df):
        self.df = df

    def logSummary(self):
        logging.info("summary input file:")
        logging.info("   %d lines" % len(self.df))
        logging.info("   %d transition groups" %  len(self.df.tg_id.unique()))
        logging.info("   %d scores including main score" %  (len(self.df.columns.values) - 6))

    def __getitem__(self, *args):
        return self.df.__getitem__(*args)

    def __setitem__(self, *args):
        return self.df.__setitem__(*args)

    def __setattr__(self, name, value):
        if name not in ["df", ]:
            raise Exception("for setting table columns use '[...]' syntax")
        object.__setattr__(self, name, value)


    def setAndRerank(self, colName, scores):
        self.df[colName] = scores
        self.rankBy(colName)

    def rankBy(self, score_col_name):
        flags = find_top_ranked(self.df.tg_num_id.values, self.df[score_col_name].values)
        self.df.is_top_peak = flags


    def decoyPeaks(self):
        return Experiment(self.df[self.df.is_decoy == True])

    def targetPeaks(self):
        return Experiment(self.df[self.df.is_decoy == False])


    def topTestPeaks(self):
        df = self.df
        return Experiment(df[(df.is_train == False) & (df.is_top_peak == True)])

    def topDecoyPeaks(self):
        ix_top = self.df.is_top_peak == True
        return Experiment(self.df[(self.df.is_decoy == True) & ix_top])

    def topTargetPeaks(self):
        ix_top = self.df.is_top_peak == True
        return Experiment(self.df[(self.df.is_decoy == False) & ix_top])


    def featureMatrix(self, use_main_score):
        min_col = 5 if use_main_score else 6
        return self.df.iloc[:, min_col:-1].values

    def filter_(self, idx):
        return Experiment(self.df[idx])

    @profile
    def add_peak_group_rank(self):
        ids = self.df.tg_num_id.values
        scores = self.df.d_score.values
        peak_group_ranks = rank(ids, scores)
        self.df["peak_group_rank"] = peak_group_ranks

    @profile
    def denoteSplitForXval(self, fraction, isTest):
        df = self.df
        decoyIDs = df[df.decoy == True].tgID.unique()
        targetIDs = df[df.decoy == False].tgID.unique()
        
        if not isTest:
            random.shuffle(decoyIDs)
            random.shuffle(targetIDs)
        else:
            decoyIDs = sorted(decoyIDs)
            targetIDs = sorted(targetIDs)
        
        decoyIDs 	= decoyIDs[:int(len(decoyIDs) * fraction) + 1]
        targetIDs 	= targetIDs[:int(len(targetIDs) * fraction) + 1]
        learnIDs 	= np.concatenate((decoyIDs, targetIDs))
        ix_learn 	= df.tgID.isin(learnIDs)
        df.trainPeak[ix_learn] = True
        df.trainPeak[~ix_learn] = False

    def trainPeaks(self):
        df = self.df[self.df.trainPeak == True]
        return Experiment(df)



class AbstractSemiSupervisedLearner(object):
    
    def startSemiSupervisedLearning(self, train):
        raise NotImplementedError()
    
    def iterSemiSupervisedLearning(self, train):
        raise NotImplementedError()
    
    def averaged_learner(self, params):
        raise NotImplementedError()
    
    def score(self, df, params):
        raise NotImplementedError()
    
    def learnRandomized(self, experiment):
        assert isinstance(experiment, Experiment)
        
        num_iter 	= CONFIG.get("semi_supervised_learner.num_iter")
        logging.info("start learnRandomized")
        
        fraction 	= CONFIG.get("xeval.fraction")
        isTest 		= CONFIG.get("is_test", False)
        experiment.denoteSplitForXval(fraction, isTest)
        train 		= experiment.trainPeaks()
        
        train.rankBy("main_score")
        
        params, clf_scores = self.startSemiSupervisedLearning(train)
        
        train.setAndRerank("classifier_score", clf_scores)
        
        # semi supervised iteration:
        for inner in range(numIter):
            params, clf_scores = self.iterSemiSupervisedLearning(train)
            train.setAndRerank("classifier_score", clf_scores)
        
        # nach semi supervsised iter: classfiy full dataset
        clf_scores = self.score(experiment, params)
        mu, nu = mean_and_std_dev(clf_scores)
        experiment.setAndRerank("classifier_score", clf_scores)
        
        td_scores = experiment.topDecoyPeaks()["classifier_score"]
        
        mu, nu = mean_and_std_dev(td_scores)
        experiment["classifier_score"] = (experiment["classifier_score"] - mu) / nu
        experiment.rank_by("classifier_score")
        
        top_test_peaks = experiment.get_top_test_peaks()
        
        top_test_target_scores = top_test_peaks.get_target_peaks()["classifier_score"]
        top_test_decoy_scores = top_test_peaks.get_decoy_peaks()["classifier_score"]
        
        logging.info("end learnRandomized")
        
        return top_test_target_scores, top_test_decoy_scores, params


class StandardSemiSupervisedLearner(AbstractSemiSupervisedLearner):
    
    def __init__(self, inner_learner):
        assert isinstance(inner_learner, AbstractLearner)
        self.inner_learner = inner_learner
    
    def select_train_peaks(self, train, sel_column, fdr, lambda_):
        assert isinstance(train, Experiment)
        assert isinstance(sel_column, basestring)
        assert isinstance(fdr, float)
        
        tt_peaks = train.get_top_target_peaks()
        tt_scores = tt_peaks[sel_column]
        td_peaks = train.get_top_decoy_peaks()
        td_scores = td_peaks[sel_column]
        
        # find cutoff fdr from scores and only use best target peaks:
        cutoff = find_cutoff(tt_scores, td_scores, lambda_, fdr)
        best_target_peaks = tt_peaks.filter_(tt_scores >= cutoff)
        return td_peaks, best_target_peaks
    
    def startSemiSupervisedLearning(self, train):
        fdr = CONFIG.get("semi_supervised_learner.initial_fdr")
        lambda_ = CONFIG.get("semi_supervised_learner.initial_lambda")
        td_peaks, bt_peaks = self.select_train_peaks(
            train, "main_score", fdr, lambda_)
        model = self.inner_learner.learn(td_peaks, bt_peaks, False)
        w = model.get_parameters()
        clf_scores = model.score(train, False)
        clf_scores -= np.mean(clf_scores)
        return w, clf_scores
    
    def iterSemiSupervisedLearning(self, train):
        fdr = CONFIG.get("semi_supervised_learner.iteration_fdr")
        lambda_ = CONFIG.get("semi_supervised_learner.iteration_lambda")
        td_peaks, bt_peaks = self.select_train_peaks( train, "classifier_score", fdr, lambda_)
        
        model = self.inner_learner.learn(td_peaks, bt_peaks, True)
        w = model.get_parameters()
        clf_scores = model.score(train, True)
        return w, clf_scores
    
    def averaged_learner(self, params):
        return self.inner_learner.averaged_learner(params)

    def score(self, df, params):
        self.inner_learner.set_parameters(params)
        return self.inner_learner.score(df, True)
