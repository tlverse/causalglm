# Helpers

family_list <- list(OR = "binomial", RR = "poisson", CATE = "gaussian")

superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(list(  Lrnr_glmnet$new(),  Lrnr_glm$new(), Lrnr_gam$new(), Lrnr_earth$new() ,
                                                                  Lrnr_ranger$new() ,  Lrnr_xgboost$new(verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,max_depth =5) ),  full_fit = TRUE) , Lrnr_cv_selector$new(loss_function_least_squares))
superlearner_RR <- make_learner(Pipeline, Lrnr_cv$new(list(  Lrnr_glmnet$new(family = "poisson"),  Lrnr_glm$new(family = poisson()), Lrnr_gam$new(family = poisson()),  
                                                             Lrnr_xgboost$new(verbose=0,max_depth =3 , objective = "count:poisson"), Lrnr_xgboost$new(verbose=0,max_depth =4, objective = "count:poisson" ), Lrnr_xgboost$new(verbose=0,max_depth =5,  objective = "count:poisson") ) , full_fit = TRUE), Lrnr_cv_selector$new(loss_function_least_squares))

learner_list_A = list(HAL = Lrnr_hal9001$new(max_degree = 2, smoothness_orders  = 1, num_knots = c(10,3)), SuperLearner = superlearner_default, glmnet =  Lrnr_glmnet$new(), glm = Lrnr_glm$new(), gam = Lrnr_gam$new(), mars = Lrnr_earth$new() ,
                      ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,max_depth =3, eval_metric = "logloss" ), Lrnr_xgboost$new(verbose=0,max_depth =4, eval_metric = "logloss" ), Lrnr_xgboost$new(verbose=0,max_depth =5, eval_metric = "logloss") ), full_fit = TRUE), Lrnr_cv_selector$new(loss_loglik_binomial)))

learner_list_Y0W = list(SuperLearner = superlearner_default, glmnet =  Lrnr_glmnet$new(), glm = Lrnr_glm$new(), gam = Lrnr_gam$new(), mars = Lrnr_earth$new() ,
                        ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,max_depth =5) ), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error)))

learner_list_Y0W_RR = list(SuperLearner = superlearner_RR, glmnet =  Lrnr_glmnet$new(family = "poisson"), glm = Lrnr_glm$new(family = poisson()), gam = Lrnr_gam$new(family = poisson()), 
                           xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,max_depth =3,  objective = "count:poisson", eval_metric = "error"), Lrnr_xgboost$new(verbose=0,max_depth =4, objective = "count:poisson", eval_metric = "error"), Lrnr_xgboost$new(verbose=0,max_depth =5, objective = "count:poisson", eval_metric = "error") ), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error)))


