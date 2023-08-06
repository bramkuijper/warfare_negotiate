#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) {
    if (!is.null(x))
    {
        sym(x)
    }
}

xvar <- "time_step"

# structure of the graph in JSON
plots_json <- paste0('[
    {"xvar" : "',xvar,'",
    "yvar" : ["mean_a_bel","mean_a_brav"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["mean_b_bel","mean_b_brav"]
    },
    {"xvar" : "',xvar,'",
    "yvar" : ["mean_belligerence","mean_bravery"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_a_bel","var_a_brav"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_b_bel","var_b_brav"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["var_belligerence","var_bravery"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : "n_escalated_contests"
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : "global_fecundity_total"
    }
    ]'
)

#file.name <- "sim_asr_20230509_101417900895_0"

# if the file name not provided raise error
if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

if (nrow(data.tibble) > 50000)
{
    data.tibble <- data.tibble[data.tibble$time %% 10 == 0,]
}

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line
        ,col_names=c("name","value")
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

plot.structure <- fromJSON(plots_json, simplifyVector = FALSE)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l))

plot.list.idx <- 1

# first plot stuff for the tibble as a whole
for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        # get rid of aes_string like this: 
        # https://stackoverflow.com/questions/74414272/how-to-replace-the-deprecated-ggplot2-function-aes-string-accepting-an-arbitrar/74414389#74414389 

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values),
                FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args))
    } else {
        
        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar
                        ),
                FUN=transform.sym)

        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

plot.title <- ""


wrap_plots(plot.list,ncol=1) + plot_annotation(
        title=plot.title)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l)
