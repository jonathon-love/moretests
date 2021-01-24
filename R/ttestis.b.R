
ttestISClass <- R6::R6Class(
    "ttestISClass",
    inherit = ttestISBase,
    public = list(
        initialize=function(...) {
            super$initialize(...)
            require('nortest')
        }
    ),
    private = list(
        .init = function() {

            if ( self$parent$options$norm) {

            table <- jmvcore::Table$new(
                options=self$options,
                name='norm2',
                title='Tests of Normality',
                columns=list(
                    list(
                        name='var',
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

            self$parent$results$assum$norm$setVisible(FALSE)
            self$parent$results$assum$insert(1, table)

            for (var in self$parent$options$vars)
                table$addRow(rowKey=var, values=list(var=var))

            table$setNote('moretests', 'Additional results provided by <em>moretests</em>')
            } 

            # redo Levene's and add ratio of variances

            if ( self$parent$options$eqv) {

            tableEqv <- jmvcore::Table$new(
                options=self$options,
                name='eqv2',
                title='Homogeneity of Variances Tests',
                columns=list(
                    list(
                        name='var',
                        title='',
                        type='text'),
                    list(
                        name='t[lv]',
                        title='',
                        type='text',
                        content="Levene's"),
                    list(
                        name='f[lv]',
                        title='F'),
                    list(
                        name='df[lv]',
                        title='df',
                        format='integer'),
                    list(
                        name='df2[lv]',
                        title='df2',
                        format='integer'),
                    list(
                        name='p[lv]',
                        title='p',
                        format='zto,pvalue'),
                    list(
                        name='t[vr]',
                        title='',
                        type='text',
                        content="Variance ratio"),
                    list(
                        name='f[vr]',
                        title='F'),
                    list(
                        name='df[vr]',
                        title='df',
                        format='integer'),
                    list(
                        name='df2[vr]',
                        title='df2',
                        format='integer'),
                    list(
                        name='p[vr]',
                        title='p',
                        format='zto,pvalue')
                )
            )

            self$parent$results$assum$eqv$setVisible(FALSE)
            self$parent$results$assum$insert(1, tableEqv)

            for (var in self$parent$options$vars)
                tableEqv$addRow(rowKey=var, values=list(var=var))

            tableEqv$setNote('moretests', 'Additional results provided by <em>moretests</em>')
            }
        },
        .run = function() {

            if (self$parent$options$norm) {

            groupVarName <- self$parent$options$group
            depVarNames <- self$parent$options$vars

            if (is.null(groupVarName) || length(depVarNames) == 0)
                return()

            data <- self$data
            table <- self$parent$results$assum$get('norm2')

            for (depName in depVarNames) {

                dataTTest <- data.frame(
                    dep=jmvcore::toNumeric(data[[depName]]),
                    group=data[[groupVarName]])

                if (self$parent$options$miss == "perAnalysis")
                    dataTTest <- jmvcore::naOmit(dataTTest)

                values <- list()
                footnote <- NULL
                # values[["name"]] <- depName

                residuals <- tapply(dataTTest$dep, dataTTest$group, function(x) x - mean(x))
                residuals <- unlist(residuals, use.names=FALSE)
                residuals <- scale(residuals)

                if (length(dataTTest$dep) < 3) {
                    next()
                }

                if (length(dataTTest$dep) > 5000) {
                    values[['s[sw]']] <- NaN
                    values[['p[sw]']] <- ''
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

                res <- try(ks.test(residuals, 'pnorm', 0, 1), silent=FALSE)

                if ( ! jmvcore::isError(res)) {
                    values[['s[ks]']] <- res$statistic
                    values[['p[ks]']] <- res$p.value
                }
                else {
                    values[['s[ks]']] <- NaN
                    values[['p[ks]']] <- ''
                }

                res <- try(nortest::ad.test(residuals), silent=FALSE)

                if ( ! jmvcore::isError(res)) {
                    values[['s[ad]']] <- res$statistic
                    values[['p[ad]']] <- res$p.value
                }
                else {
                    values[['s[ad]']] <- NaN
                    values[['p[ad]']] <- ''
                }

                table$setRow(rowKey=depName, values)
            }
            }

            if (self$parent$options$eqv) {

            groupVarName <- self$parent$options$group
            depVarNames <- self$parent$options$vars

            if (is.null(groupVarName) || length(depVarNames) == 0)
                return()

            data <- self$data
            tableEqv <- self$parent$results$assum$get('eqv2')

            for (depName in depVarNames) {

                dataTTest <- data.frame(
                    dep=jmvcore::toNumeric(data[[depName]]),
                    group=data[[groupVarName]])

                if (self$parent$options$miss == "perAnalysis")
                    dataTTest <- jmvcore::naOmit(dataTTest)

                values <- list()
                footnote <- NULL

                levene <- try(car::leveneTest(dep ~ group, data=dataTTest, "mean"), silent=TRUE)

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

                dep <- dataTTest$dep
                group <- as.numeric(factor(dataTTest$group))

                res <- try(var.test(dep[group==1], dep[group==2]))

                if ( ! jmvcore::isError(res) && ! is.na(res$statistic)) {
                    values[['f[vr]']] <- res$statistic
                    values[['df[vr]']] <- res$parameter[1]
                    values[['df2[vr]']] <- res$parameter[2]
                    values[['p[vr]']] <- res$p.value
                }
                else {
                    values[['f[vr]']] <- NaN
                    values[['df[vr]']] <- ""
                    values[['df2[vr]']] <- ""
                    values[['p[vr]']] <- ""
                }

                tableEqv$setRow(rowKey=depName, values)

                if (jmvcore::isError(levene) || is.na(levene[1,"F value"]) || jmvcore::isError(res) || is.na(res$statistic))
                    tableEqv$addFootnote(rowKey=depName, "f", "F-statistic could not be calculated")
            }
            }
        }
    )
)
