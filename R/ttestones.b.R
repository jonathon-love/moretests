
ttestOneSClass <- R6::R6Class(
    "ttestOneSClass",
    inherit = ttestOneSBase,
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

            self$parent$results$normality$setVisible(FALSE)
            self$parent$results$insert(2, table)

            for (var in self$parent$options$vars)
                table$addRow(rowKey=var, values=list(var=var))

            table$setNote('moretests', 'Additional results provided by <em>moretests</em>')
        },
        .run = function() {

            if ( ! self$parent$options$norm)
                return()

            data <- self$data
            if (self$parent$options$miss == 'listwise')
                data <- jmvcore::naOmit(data)

            table <- self$parent$results$get('norm2')

            for (name in self$parent$options$vars) {

                column <- jmvcore::toNumeric(data[[name]])
                column <- na.omit(column)
                residuals <- scale(column)
                n <- length(residuals)
                values <- list()

                swf <- NULL

                if (n < 3) {
                    swf <- list(rowKey=name, "s[sw]", "Too few observations (N < 3) to compute statistic")
                    values[['s[sw]']] <- NaN
                    values[['p[sw]']] <- ''
                } else if (n > 5000) {
                    swf <- list(rowKey=name, "s[sw]", "Too many observations (N > 5000) to compute statistic")
                    values[['s[sw]']] <- NaN
                    values[['p[sw]']] <- ''
                } else {
                    res <- try(shapiro.test(residuals), silent=TRUE)
                    if (jmvcore::isError(res)) {
                        values[['s[sw]']] <- NaN
                        values[['p[sw]']] <- ''
                    } else {
                        values[['s[sw]']] <- res$statistic
                        values[['p[sw]']] <- res$p.value
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

                table$setRow(rowKey=name, values)
                if ( ! is.null(swf))
                    do.call(table$addFootnote, swf)
            }
        })
)
