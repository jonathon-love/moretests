
linRegClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "linRegClass",
    inherit = linRegBase,
    public = list(
        initialize=function(...) {
            super$initialize(...)
            require('nortest')
            require('lmtest')
        }
    ),
    private = list(
        .init = function() {

            if (! self$parent$options$norm)
                return()
 
            termsAll <- private$.getModelTerms()
            for (i in seq_along(termsAll)) {

                thomo <- jmvcore::Table$new(
                    options=self$options,
                    name='homo',
                    title='Heteroskedasticity Tests',
                    rows=1,
                    clearWith=list(
                        "dep",
                        "blocks"),
                    columns=list(
                        list(
                            name='t[bp]',
                            title='',
                            type='text',
                            content='Breusch-Pagan'),
                        list(
                            name='s[bp]',
                            title='Statistic'),
                        list(
                            name='p[bp]',
                            title='p',
                            format='zto,pvalue'),
                        list(
                            name='t[gq]',
                            title='',
                            type='text',
                            content='Goldfeld-Quandt'),
                        list(
                            name='s[gq]',
                            title='Statistic'),
                        list(
                            name='p[gq]',
                            title='p',
                            format='zto,pvalue'),
                        list(
                            name='t[hm]',
                            title='',
                            type='text',
                            content='Harrison-McCabe'),
                        list(
                            name='s[hm]',
                            title='Statistic'),
                        list(
                            name='p[hm]',
                            title='p',
                            format='zto,pvalue')
                    )
                )

                self$parent$results$models$get(key=i)$assump$insert(1, thomo)
                thomo$addRow(1, values=list(
                    `t[sw]`="Breusch-Pagan",
                    `t[ks]`="Goldfeld-Quandt",
                    `t[ad]`="Harrison-McCabe"))

                thomo$setNote('moretests', 'Additional results provided by <em>moretests</em>')
 
                table <- jmvcore::Table$new(
                    options=self$options,
                    name='norm2',
                    title='Normality Tests',
                    rows=1,
                    clearWith=list(
                        "dep",
                        "blocks"),
                    columns=list(
                        list(
                            name='t[sw]',
                            title='',
                            type='text',
                            content='Shapiro-Wilk'),
                        list(
                            name='s[sw]',
                            title='Statistic'),
                        list(
                            name='p[sw]',
                            title='p',
                            format='zto,pvalue'),
                        list(
                            name='t[ks]',
                            title='',
                            type='text',
                            content='Kolmogorov-Smirnov'),
                        list(
                            name='s[ks]',
                            title='Statistic'),
                        list(
                            name='p[ks]',
                            title='p',
                            format='zto,pvalue'),
                        list(
                            name='t[ad]',
                            title='',
                            type='text',
                            content='Anderson-Darling'),
                        list(
                            name='s[ad]',
                            title='Statistic'),
                        list(
                            name='p[ad]',
                            title='p',
                            format='zto,pvalue')
                    )
                )

                self$parent$results$models$get(key=i)$assump$norm$setVisible(FALSE)
                self$parent$results$models$get(key=i)$assump$insert(1, table)              
                table$addRow(1, values=list(
                    `t[sw]`="Shapiro-Wilk",
                    `t[ks]`="Kolmogorov-Smirnov",
                    `t[ad]`="Anderson-Darling"))

                table$setNote('moretests', 'Additional results provided by <em>moretests</em>')
          }
        },
        .run = function() {

            if (! self$parent$options$norm)
                return()

            groups <- self$parent$results$models
            termsAll <- private$.getModelTerms()
            .residuals <- private$.computeResiduals()

            for (i in seq_along(termsAll)) {
                residuals <- .residuals[[i]]
                residuals <- scale(residuals)
                table <- groups$get(key=i)$assump$get('norm2')

                res <- try(shapiro.test(residuals), silent=TRUE)

                if (jmvcore::isError(res)) {
                    values <- list(`s[sw]`=NaN, `p[sw]`='')
                } else {
                    values <- list(`s[sw]`=res$statistic, `p[sw]`=res$p.value)
                }

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


                # Heteroskedasticity Tests
                thomo <- groups$get(key=i)$assump$get('homo')
                model <- self$parent$models[[i]]

                # Breusch-Pagan   
                res <- try(lmtest::bptest(model), silent=TRUE) 
                if (jmvcore::isError(res)) {
                    values <- list(`s[bp]`=NaN, `p[bp]`='')
                } else {
                    values <- list(`s[bp]`=res$statistic, `p[bp]`=res$p.value)
                }
                # Goldfeld-Quandt  
                res <- try(lmtest::gqtest(model), silent=TRUE)
                if ( ! jmvcore::isError(res)) {
                    values[['s[gq]']] <- res$statistic
                    values[['p[gq]']] <- res$p.value
                }
                else {
                    values[['s[gq]']] <- NaN
                    values[['p[gq]']] <- ''
                }

                # Harrison-McCabe   
                res <- try(lmtest::hmctest(model), silent=TRUE)
                if ( ! jmvcore::isError(res)) {
                    values[['s[hm]']] <- res$statistic
                    values[['p[hm]']] <- res$p.value
                }
                else {
                    values[['s[hm]']] <- NaN
                    values[['p[hm]']] <- ''
                }

                thomo$setRow(rowNo=1, values)
            }
        },
        .modelTerms = NULL,
        .getModelTerms = function() {
            if (is.null(private$.modelTerms)) {
                blocks <- self$parent$options$blocks

                terms <- list()
                if (is.null(blocks)) {
                    terms[[1]] <- c(self$parent$options$covs, self$parent$options$factors)
                } else {
                    for (i in seq_along(blocks)) {
                        terms[[i]] <- unlist(blocks[1:i], recursive = FALSE)
                    }
                }
                private$.modelTerms <- terms
            }

            return(private$.modelTerms)
        },
        .computeResiduals = function() {
            res <- list()
            for (i in seq_along(self$parent$models))
                res[[i]] <- self$parent$models[[i]]$residuals

            return(res)
        }
    )
)