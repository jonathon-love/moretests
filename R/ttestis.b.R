

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

            if ( ! self$parent$options$norm)
                return()

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
        },
        .run = function() {

            if ( ! self$parent$options$norm)
                return()

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
    )
)
