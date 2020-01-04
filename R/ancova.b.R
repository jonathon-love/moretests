
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

            if ( ! self$parent$options$norm)
                return()

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
        },
        .run = function() {

            if ( ! self$parent$options$norm)
                return()

            residuals <- self$parent$residuals
            if (is.null(residuals))
                return()

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
        })
)
