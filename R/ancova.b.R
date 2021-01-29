
ancovaClass <- R6::R6Class(
    "ancovaClass",
    inherit = ancovaBase,
    public = list(
        initialize=function(...) {
            super$initialize(...)
            require('nortest')
        }
    ),
    private = list(
        .init = function() {

            if (self$parent$options$norm){
 
                table <- self$parent$results$assump$norm

                table$setTitle('Normality tests')
                table$getColumn('t[sw]')$setVisible(TRUE)

                table$addColumn(
                    name='t[ks]',
                    title='',
                    type='text')
                table$addColumn(
                    name='s[ks]',
                    title='statistic')
                table$addColumn(
                    name='p[ks]',
                    title='p',
                    format='pvalue')

                table$addColumn(
                    name='t[ad]',
                    title='',
                    type='text')
                table$addColumn(
                    name='s[ad]',
                    title='statistic')
                table$addColumn(
                    name='p[ad]',
                    title='p',
                    format='pvalue')

                table$setRow(rowNo=1, values=list(
                    `t[ks]`='Kolmogorov-Smirnov',
                    `t[ad]`='Anderson-Darling'))

                table$setNote('moretests', 'Additional results provided by <em>moretests</em>')
            }

            # redo Levene's and add Bartlett

            if ( self$parent$options$homo) {

                thomo <- jmvcore::Table$new(
                    options=self$options,
                    name='homo2',
                    title='Homogeneity of Variances Tests',
                    columns=list(
                        list(
                            name='t[lv]',
                            title='',
                            type='text',
                            content="Levene's"),
                        list(
                            name='f[lv]',
                            title='Statistic'),
                        list(
                            name='df[lv]',
                            title='df',
                            type='integer'),
                        list(
                            name='df2[lv]',
                            title='df2',
                            type='integer'),
                        list(
                            name='p[lv]',
                            title='p',
                            format='zto,pvalue'),
                        list(
                            name='t[bt]',
                            title='',
                            type='text',
                            content="Bartlett's"),
                        list(
                            name='f[bt]',
                            title='Statistic'),
                        list(
                            name='df[bt]',
                            title='df',
                            type='integer'),
                        list(
                            name='df2[bt]',
                            title='df2'),
                        list(
                            name='p[bt]',
                            title='p',
                            format='zto,pvalue')
                    )
                )

                self$parent$results$assump$homo$setVisible(FALSE)
                self$parent$results$assump$insert(1, thomo)

                thomo$addRow(1, values=list(
                    `t[lv]`="Levene's",
                    `t[bt]`="Bartlett's"))

                thomo$setNote('moretests', 'Additional results provided by <em>moretests</em>')
            }
        },
        .run = function() {

            if (self$parent$options$norm) {

                residuals <- self$parent$residuals
                if (! is.null(residuals)) {

                    table <- self$parent$results$assump$norm
                    residuals <- scale(residuals)
                    values <- list()

                    res <- try(ks.test(residuals, 'pnorm', 0, 1), silent=TRUE)

                    if ( ! jmvcore::isError(res)) {
                        values[['s[ks]']] <- res$statistic
                        values[['p[ks]']] <- res$p.value
                    }
                    else {
                        values[['s[ks]']] <- NaN
                        values[['p[ks]']] <- ''
                    }

                    res <- try(nortest::ad.test(residuals), silent=TRUE)

                    if ( ! jmvcore::isError(res)) {
                        values[['s[ad]']] <- res$statistic
                        values[['p[ad]']] <- res$p.value
                    }
                    else {
                        values[['s[ad]']] <- NaN
                        values[['p[ad]']] <- ''
                    }

                    table$setRow(rowNo=1, values)
                }
            }

            # Homogeneity

            if (self$parent$options$homo) {

                thomo <- self$parent$results$assump$get('homo2')

                data <- self$data

                dep <- self$parent$options$dep
                factors <- self$parent$options$factors

                # Levene's
                modelTerms <- private$.modelTerms()
                formula <- jmvcore::constructFormula(dep, modelTerms)
                formula <- stats::as.formula(formula)
                model <- stats::aov(formula, data)
                residuals <- model$residuals

                data[[".RES"]] <- abs(residuals)

                rhs <- paste0('`', factors, '`', collapse=':')
                formula <- as.formula(paste0('.RES ~', rhs))

                model <- stats::lm(formula, data)
                summary <- stats::anova(model)

                values <- list()
                footnote <- NULL

                if (is.na(summary[1,"F value"])) {
                    values[['f[lv]']] <- NaN
                    values[['df[lv]']] <- ""
                    values[['df2[lv]']] <- ""
                    values[['p[lv]']] <- ""
                } else {
                    values[['f[lv]']] <- summary[1,"F value"]
                    values[['df[lv]']] <- summary[1,"Df"]
                    values[['df2[lv]']] <- summary[2,"Df"]
                    values[['p[lv]']] <- summary[1,"Pr(>F)"]
                }
                
                # Bartlett's
                formula <- as.formula(paste0(dep, '~', rhs))
                res <- try(private$.bartlett(formula, data))

                if ( ! jmvcore::isError(res) && ! is.na(res$statistic)) {
                    values[['f[bt]']] <- res$statistic
                    values[['df[bt]']] <- res$parameter
                    values[['df2[bt]']] <- ""
                    values[['p[bt]']] <- res$p.value
                }
                else {
                    values[['f[bt]']] <- NaN
                    values[['df[bt]']] <- ""
                    values[['df2[bt]']] <- ""
                    values[['p[bt]']] <- ""
                }

                thomo$setRow(rowNo=1, values)

            }
        },
        .bartlett = function (formula, data) {
            # modified from stats:::bartlett.test.formula
            m <- match.call()
            if (is.matrix(eval(m$data, parent.frame())))
                m$data <- as.data.frame(data)

            m[[1L]] <- quote(stats::model.frame)
            mf <- eval(m, parent.frame())
            DNAME <- as.character(formula)[3]
            if (length(mf) != 2L){
                fac<-do.call(paste, c(mf[,-1], sep='-'))
                mf<-data.frame(mf[,1],fac)
            }
            names(mf) <- NULL
            y <- do.call("bartlett.test", as.list(mf))
            y$data.name <- DNAME
            y
        },
        .modelTerms=function() {
            modelTerms <- self$parent$options$modelTerms
            if (length(modelTerms) == 0)
                modelTerms <- private$.ff()

            lengths <- vapply(modelTerms, length, 1)
            modelTerms <- modelTerms[order(lengths)]

            modelTerms
        },
        .ff=function() {
            factors <- self$parent$options$factors
            if (length(factors) > 1) {
                formula <- as.formula(paste('~', paste(paste0('`', factors, '`'), collapse='*')))
                terms   <- attr(stats::terms(formula), 'term.labels')
                modelTerms <- sapply(terms, function(x) as.list(strsplit(x, ':')), USE.NAMES=FALSE)
            } else {
                modelTerms <- as.list(factors)
            }
            for (i in seq_along(modelTerms)) {
                term <- modelTerms[[i]]
                quoted <- grepl('^`.*`$', term)
                term[quoted] <- substring(term[quoted], 2, nchar(term[quoted])-1)
                modelTerms[[i]] <- term
            }
            covs <- NULL
            if ('covs' %in% names(self$parent$options))
                covs <- self$parent$options$covs

            for (covariate in covs)
                modelTerms[[ length(modelTerms) + 1 ]] <- covariate

            modelTerms
        }
    )
)