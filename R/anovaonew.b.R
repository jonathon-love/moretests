
anovaOneWClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "anovaOneWClass",
    inherit = anovaOneWBase,
    public = list(
        initialize=function(...) {
            super$initialize(...)
            require('nortest')
        }
    ),
    private = list(
        .init = function() {

            if (self$parent$options$norm){
 
                 table <- jmvcore::Table$new(
                    options=self$options,
                    name='norm2',
                    title='Normality Tests',
                    columns=list(
                        list(
                            name='dep',
                            title='',
                            type='text'),
                        list(
                            name='t[sw]',
                            title='',
                            type='text',
                            content='Shapiro-Wilk'),
                        list(
                            name='s[sw]',
                            title='statistic'),
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
                            title='statistic'),
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
                            title='statistic'),
                        list(
                            name='p[ad]',
                            title='p',
                            format='zto,pvalue')
                    )
                )

                self$parent$results$assump$norm$setVisible(FALSE)
                self$parent$results$assump$insert(1, table)

                for (dep in self$parent$options$deps)
                    table$addRow(rowKey=dep, values=list(dep=dep))

                table$setNote('moretests', 'Additional results provided by <em>moretests</em>')
           }

            if (self$parent$options$eqv) {

                tableEqv <- jmvcore::Table$new(
                    options=self$options,
                    name='eqv2',
                    title='Homogeneity of Variances Tests',
                    columns=list(
                        list(
                            name='dep',
                            title='',
                            type='text'),
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
                            name='t[vr]',
                            title='',
                            type='text',
                            content="Bartlett's"),
                        list(
                            name='f[vr]',
                            title='Statistic'),
                        list(
                            name='df[vr]',
                            title='df',
                            type='integer'),
                        list(
                            name='df2[vr]',
                            title='df2'),
                        list(
                            name='p[vr]',
                            title='p',
                            format='zto,pvalue')
                    )
                )

                self$parent$results$assump$eqv$setVisible(FALSE)
                self$parent$results$assump$insert(1, tableEqv)

                for (dep in self$parent$options$deps)
                    tableEqv$addRow(rowKey=dep, values=list(dep=dep))

                tableEqv$setNote('moretests', 'Additional results provided by <em>moretests</em>')
            }
        },
        .run = function() {

            groupVarName <- self$parent$options$group
            depVarNames <- self$parent$options$deps
            data <- self$data

            if (is.null(groupVarName) || length(depVarNames) == 0)
                return()

            # Normality

            if (self$parent$options$norm) {

                table <- self$parent$results$assump$get('norm2')

                for (depName in depVarNames) {

                    dataOWTest <- data.frame(
                        dep=jmvcore::toNumeric(data[[depName]]),
                        group=data[[groupVarName]])

                    if (self$parent$options$miss == "perAnalysis")
                        dataOWTest <- jmvcore::naOmit(dataOWTest)

                    residuals <- rstandard(lm(dep ~ group, data=dataOWTest))
                    if (! is.null(residuals)) {

                        residuals <- scale(residuals)
                        values <- list()
                        footnote <- NULL

                        if (length(dataOWTest$dep) < 3) {
                            values[['s[sw]']] <- NaN
                            values[['p[sw]']] <- ''
                            footnote <- 'Too few samples to compute statistic (N < 3)'
                        }

                        if (length(dataOWTest$dep) > 5000) {
                            values[['s[sw]']] <- NaN
                            values[['p[sw]']] <- ''
                            footnote <- 'Too many samples to compute statistic (N > 5000)'
                        } else {
                            res <- try(shapiro.test(residuals), silent=TRUE)

                            if ( ! jmvcore::isError(res)) {
                                values[['s[sw]']] <- res$statistic
                                values[['p[sw]']] <- res$p.value
                            }
                            else {
                                values[['s[sw]']] <- NaN
                                values[['p[sw]']] <- ''
                            }
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

                        table$setRow(rowKey=depName, values)

                        if ( ! is.null(footnote))
                            table$addFootnote(rowKey=dep, 's', footnote)
                    }
                }
            }

            # Homogeneity

            if (self$parent$options$eqv) {

                tableEqv <- self$parent$results$assump$get('eqv2')

                for (depName in depVarNames) {

                    dataOWTest <- data.frame(
                        dep=jmvcore::toNumeric(data[[depName]]),
                        group=data[[groupVarName]])

                    if (self$parent$options$miss == "perAnalysis")
                        dataOWTest <- jmvcore::naOmit(dataOWTest)

                    values <- list()    

                    levene <- try(car::leveneTest(dep ~ group, data=dataOWTest, "mean"), silent=TRUE)

                    if (jmvcore::isError(levene) || is.na(levene[1,"F value"])) {
                        values[['f[lv]']] <- NaN
                        values[['df[lv]']] <- ""
                        values[['df2[lv]']] <- ""
                        values[['p[lv]']] <- ""
                    } else {
                        values[['f[lv]']] <- levene[1,"F value"]
                        values[['df[lv]']] <- levene[1,"Df"]
                        values[['df2[lv]']] <- levene[2,"Df"]
                        values[['p[lv]']] <- levene[1,"Pr(>F)"]
                    }

                    res <- try(bartlett.test(dep ~ group, data = dataOWTest), silent=TRUE)

                    if ( ! jmvcore::isError(res) && ! is.na(res$statistic)) {
                        values[['f[vr]']] <- res$statistic
                        values[['df[vr]']] <- res$parameter
                        values[['df2[vr]']] <- ""
                        values[['p[vr]']] <- res$p.value
                    }
                    else {
                        values[['f[vr]']] <- NaN
                        values[['df[vr]']] <- ""
                        values[['df2[vr]']] <- ""
                        values[['p[vr]']] <- ""
                    }

                    tableEqv$setRow(rowKey=depName, values)

                   if (jmvcore::isError(levene) || is.na(levene[1,"F value"]))
                        tableEqv$addFootnote(rowKey=depName, "f[lv]", "F-statistic could not be calculated")
                    if (jmvcore::isError(res) || is.na(res$statistic))
                        tableEqv$addFootnote(rowKey=depName, "f[vr]", "F-statistic could not be calculated")
                }
            }
        }
    )
)