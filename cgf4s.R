# code by Rodrigo Miranda (rodrigo.qmiranda@gmail.com)
# and Josicleda Galvincio (josicleda@gmail.com)

for (pkg in c("future", "RANN", "dplyr", "caret", "Cubist", 'qmap', 'NMOF')) {
    suppressMessages(suppressWarnings(library(pkg, character.only=TRUE)))
}


dataset.load <- function(files, header=TRUE) {
    # this loads all data from txt files
    # into a list of dataframe objects

    # RECENT CHANGES #
    # 1. changed read.table(file=x, header=FALSE, col.names=c( sub('\\..*$','', basename(x)) )) to
    # read.delim(file=x, header=FALSE, quote = "\"", sep=',', dec = '.', col.names=c( sub('\\..*$','', basename(x)) ))
    #
    # 2. this next function was completely rewritten
    files_df <- parallel::mcmapply(function(x) {y <- read.delim(file=x, header=FALSE, quote = "\"", sep=',', dec = '.')
                                                y <- y[, 1] # 1 for tmax and 2 for tmin
                                                y <- data.frame(y)
                                                colnames(y) <- c( sub('\\..*$','', basename(x)) )
                                                return(y)},
                                   files,
                                   USE.NAMES = FALSE,
                                   mc.cores = parallel::detectCores(),
                                   mc.cleanup = TRUE)

    # this merges all dataframes
    # into a single one
    combined_df <- do.call("cbind", files_df)
    if (header == TRUE) {     # removes the first line which is just the starting date (YYYYMMMDD)
        combined_df <- combined_df[2:nrow(combined_df),]
    }
    combined_df[combined_df == -99] <- NA    # sets -99 as NA

    # this is to release memory
    rm(list=setdiff(ls(), "combined_df")); gc()

    return(combined_df)
}


files.order <- function(files) {
    # this loads all data from txt files
    # into a list of dataframe objects

    # RECENT CHANGES #
    # 1. changed y <- read.table(file=x, header=FALSE) to
    # y <- read.delim(file=f, header=FALSE, quote = "\"", sep=',', dec = '.')
    #
    # 2. changed y <- y[2:nrow(y),]
    # to y <- y[2:nrow(y),1]
    files.na_pct <- parallel::mcmapply(function(x) {y <- read.delim(file=x, header=FALSE, quote = "\"", sep=',', dec = '.')
                                                    y <- y[2:nrow(y), 1] # 1 for tmax and 2 for tmin
                                                    y[y == -99] <- NA
                                                    na_pct <- sum(is.na(y))/length(y)
                                                    return(na_pct)},
                                       files,
                                       USE.NAMES = FALSE,
                                       mc.cores = parallel::detectCores(),
                                       mc.cleanup = TRUE)

    files <- files[order(files.na_pct)]

    # this is to release memory
    rm(list=setdiff(ls(), "files")); gc()

    return(files)
}


fork.load <- function(fork) {
    f <- read.table(file=fork, header=TRUE, sep=',')

    return(f)
}


model.nn <- function(n, m, x, y, train=FALSE) {
    n1 <- n[1]
    n2 <- n[2]

    y_mod.n1 = predict(m, x, neighbors=n1)
    y_mod.n2 = predict(m, x, neighbors=n2)

    if (max(y_mod.n1) > 0) {
        for (i in 1:2) {
            rl <- rle(y_mod.n1)
            rl.lst <- lapply(split(seq_along(y_mod.n1), rep(seq_along(rl$values), rl$lengths)), range)
            rl.df <- data.frame(do.call(rbind, rl.lst))
            colnames(rl.df) <- c("begin", "end")
            rl.df$diff <- rl.df$end - rl.df$begin + 1
            rl.df$values <- rl$values
            rl.df <- rl.df[rl.df$values>0 & rl.df$diff>1,]

            if (nrow(rl.df) > 0) {
                for(j in 1:nrow(rl.df)) {
                    r <- rl.df[j,]
                    if (i == 1) {
                        y_mod.n1[r$begin:r$end] <- y_mod.n2[r$begin:r$end]
                    } else {
                        y_mod.n1[r$begin:r$end] <- 0
                    }
                }
            }
        }
    }

    if (train == TRUE) {
        return(cor(y, y_mod.n1, use="complete.obs", method="pearson"))
    } else {
        return(y_mod.n1)
    }
}


model.run <- function(f) {
    # this temporarily changes
    # the current directory
    # from a specific dir related
    # to that iteration
    wd <- getwd()
    setwd(file.path(wd, iter))

    # this gets the basename of the file
    # without extension
    var <- sub('\\..*$','', basename(f))

    if (file.exists(basename(f))) {
        print(c("Importing: ", var))
        y <- read.table(file=basename(f), header=FALSE)
        y[y == -99] <- NA    # sets -99 as NA
        combined_df[,var] <<- y     # updates the main dataset
        names(y) <- var
        setwd(wd)
        rm(list=setdiff(ls(), "y")); gc()   # this is to release memory
        return(y)
    } else {
        # this subsets the fork
        # dataset based on the
        # name of the station
        sub_forks <- subset(forks, Name==var)
        id <- sub_forks[1,1]  # this transforms the real id into an index (there was a +1 here before)

        # this subssets the closest
        # dataset based on the id in
        # the column X1
        sub_closest <- closest[id,,drop=F]  # this ',,' is important to keep
        ids <- sub_closest[1,1:k]  # this transforms the indices into real ids (as a vector) (there was a -1 here before)

        # this creates the vector
        # with the names of all
        # closest station to the target
        sub_forks <- subset(forks, ID %in% ids)
        nam <- as.vector(sub_forks$Name)

        # this pre-processes the
        # combined_df into a new dataset
        # with only the closest stations
        n_combined_df <- dplyr::select(combined_df, all_of(nam))

        # this creates the x matrix
        # used by the model
        # x <- dplyr::select(n_combined_df, -one_of(var))
        if (iter == 1) {     # this conditional keeps the modelled data across iterations (2nd experiment)
            x <- dplyr::select(n_combined_df, -one_of(var))
        } else {
            x <- n_combined_df
        }

        # this creates the y vector
        # used by the model

        # changed y <- read.table(file=f, header=FALSE)
        # to y <- read.delim(file=x, header=FALSE, quote = "\"", sep=',', dec = '.')
        y <- read.delim(file=f, header=FALSE, quote = "\"", sep=',', dec = '.')
        # y <- y[2:nrow(y),]    # removes the first line which is just the starting date (YYYYMMMDD)
        y <- y[2:nrow(y), 1] # 1 for tmax and 2 for tmin
        y[y == -99] <- NA    # sets -99 as NA

        # this inserts the months as
        # as a column in the dataset
        # I guess this is the 1st exp
        dts <- seq(idt, idt+nrow(x)-1, "days")  # this could be weekly instead of monthly
        mon <- as.numeric(format(dts, "%m"))
        x$mon <- mon

        # # if yearly precitation was less than
        # # 50 mm, this puts NA to whole year.
        # # this idea is based on the thought that
        # # some stations has more zeros than should
        # # creating some very dry years that are
        # # actually wrong data (50 mm is desert climate)
        # yr <- as.numeric(format(dts, "%Y"))
        # agr <- aggregate(y, by=list(yr), sum)
        # agr <- subset(agr, agr['x'] < 50)[,1]
        # y <- replace(y, yr %in% agr, NA)

        try ({
            # this counts the NA
            # in the target vector and
            # cancels processing if
            # higher than a threshold
            na_pct <- sum(is.na(y))/length(y)
            if (na_pct > 0.99) {
                stop("NA count is too high")
            }

            # this removes NA values based on
            # ths observed dataset for training
            x_train <- x[!is.na(y),]
            y_train <- y[!is.na(y)]

            # this part uses only the
            # cubist module for training
            ctrl = cubistControl(unbiased=FALSE,
                                rules=10,
                                extrapolation=5,
                                label="outcome")
            m <- cubist(x=x_train, y=y_train, control=ctrl)     # , committees=100)

            # this is to export various
            # types of results
            cat(m$output, file=sub('.txt', '_out.txt', basename(f)))

            # for (i in 1:2) {
            #     if (exists("y_mod")) {
            #         if (round(cor(y, y_mod, use="complete.obs", method="pearson"), 2) < 0.99) {
            #             print("Running deeper calibration")
            #             n1 <- seq(1, 9, 1)
            #         } else {
            #             break
            #         }
            #     } else if (i == 1) {
            #         n1 <- c(1)
            #     } else {
            #         break
            #     }

            n1 <- c(1)
            n <- gridSearch(fun=model.nn,
                            levels=list(n1=n1, n2=seq(1, 9, 1)),
                            m=m, x=x, y=y, train=TRUE,
                            printDetail=FALSE,
                            method='multicore',
                            mc.control=list(mc.cores = parallel::detectCores(), mc.cleanup = TRUE))
        
            n.best <- match(max(n$values), n$values)
            n <- n$levels[n.best]
            n <- unlist(n)

            y_mod <- model.nn(n, m, x, y)

            print(c(var, n, round(cor(y, y_mod, use="complete.obs", method="pearson"),2)))
            # }

            if (max(y_mod) > 0) {
                qm_fit <- qmap::fitQmap(y, y_mod,
                                        method="PTF",
                                        transfun="scale",
                                        cost="RSS",
                                        qstep=0.01,    # (<- this could be calibrated)
                                        wett.day=FALSE)     # <- TRUE for precipitation
                y_mod <- qmap::doQmap(y_mod, qm_fit)

                # # this calculates and
                # # exports stats
                # # regarding the simulation
                # # of dry and wet days
                # y.bin <- ifelse(y > 0, 1, 0)[!is.na(y)]
                # y_mod.bin <- ifelse(y_mod > 0, 1, 0)[!is.na(y)]
                # u <- union(y_mod.bin, y.bin)
                # if (length(unique(u)) > 1) {
                #     result <- confusionMatrix(table(factor(y_mod.bin, u), factor(y.bin, u)))$byClass
                #     result.names <- paste(names(result), collapse=",")
                #     result.data <- paste(result, collapse=",")
                #     cat(paste(c(result.names, result.data), collapse = "\n"), file=sub('.txt', '_bin.txt', basename(f)))
                # }
            }
        })



            # print(c('var'=basename(f),
            #         'na_pct'=round(sum(is.na(y))/length(y), 2),
            #         'n'=n,
            #         'cor'=round(cor(y, y_mod, use="complete.obs", method="pearson"), 2),
            #         'rmse'=round(sqrt(mean((y-y_mod)^2, na.rm=TRUE)), 2),
            #         'pct_rmse'=round(sqrt(mean((y-y_mod)^2, na.rm=TRUE))/mean(y, na.rm=TRUE), 2)))


            # y_mod <- ifelse(is.na(y), y_mod, y)


        if (exists("y_mod")) {
            cat(y_mod, sep="\n", file=basename(f))
            combined_df[,var] <<- y_mod     # updates the main dataset
            y_mod <- data.frame("l"=unlist(y_mod))
            names(y_mod) <- var
            setwd(wd)
            rm(list=setdiff(ls(), "y_mod")); gc()   # this is to release memory
            return(y_mod)
        } else {
            cat(y, sep="\n", file=basename(f))
            y <- data.frame("l"=unlist(y))
            y[y == NA] <- -99    # sets -99 as NA
            names(y) <- var
            setwd(wd)
            rm(list=setdiff(ls(), "y")); gc()   # this is to release memory
            return(y)
        }
    }
}


main <- function(dir, idt) {
    # makes the date a
    # global variable
    idt <<- idt

    # number of nn that
    # each station will
    # get to use in modelling
    k <<- 31

    # read the fork and
    # calculates all nn
    forks <<- fork.load(file.path(dir, "tmp.txt"))  # changed "pcp.txt" to "tmp.txt"
    forks[,1] <<- rownames(forks)
    write.csv(forks, 'forks.csv', row.names=TRUE)
    closest <<- nn2(data=forks[,3:4], k=k)[[1]]
    closest <<- data.frame(closest)
    write.csv(closest, 'closest.csv', row.names=TRUE)

    # this makes a list of all fullpath
    # filenames in the directory 'dir',
    # while excluding forks
    files <- c(list.files(path=dir, pattern='^t.*\\.txt$', full.names=TRUE))    # changed '^p.*\\.txt$' to '^t.*\\.txt$'
    files <- files[!files %in% file.path(dir, "tmp.txt")]   # changed "pcp.txt" to "tmp.txt"
    names(files) <- sapply(files, function(x) {sub('\\..*$','', basename(x))})

    # this loads all data from txt files
    # into a list of dataframe objects
    print("Sorting stations by missing values")
    files <- files.order(files)

    # this reads multiple files
    # into a single dataframe
    # in parallel
    print("Loading dataset")
    combined_df <<- dataset.load(files)
    combined_df <<- data.frame(combined_df)

    # this is to release memory
    rm(list=setdiff(ls(), "files")); gc()

    for(i in 1:3) {
        # this makes de iter
        # accessable from within
        # the model function
        iter <<- i

        # this prints the number
        # of the iteration
        print(paste(c("Processing iteration", iter), collapse = " "))

        # this creates the
        # folder that will
        # hold the interation
        # results
        dir.create(file.path(getwd(), iter), showWarnings = FALSE)

        # this process cubist in
        # for each file
        # combined_df <<- lapply(files, model.run)
        # combined_df <<- do.call("cbind", combined_df)
        # combined_df[combined_df == -99] <<- NA    # sets -99 as NA
        # combined_df <<- data.frame(combined_df)
        for (f in files) {
            model.run(f)
        }
        
        # this is to release memory
        gc()
    }

    return(NULL)
}


dir <- "/mnt/d/SUPer/Clima/EST_1961_a_2021_3_ComRAD/3_final"     # Linux version
# dir <- "D:\SUPer\Clima\EST_1961_a_2021_3.ComRAD\3_final"    # Windows version
main(dir, idt=as.Date("1961-01-01"))
