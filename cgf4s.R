# code by Rodrigo Miranda (rodrigo.qmiranda@gmail.com)
# and Josicleda Galvincio (josicleda@gmail.com)

library("future")
library("RANN")
library("dplyr")
library("Cubist")

make_dataset <- function(files, cl) {
    # this loads all data from txt files
    # into a list of dataframe objects
    files_df <- parallel::parLapply(cl, files, fun=function(x) {read.table(file=x, header=FALSE, col.names=c( sub('\\..*$','', basename(x)) ))})

    # this merges all dataframes
    # into a single one
    combined_df <- do.call("cbind", files_df)
    combined_df <- combined_df[2:nrow(combined_df),]    # removes the first line which is just the starting date (YYYYMMMDD)
    combined_df[combined_df == -99] <- NA    # sets -99 as NA

    # this is to release memory
    rm(list=setdiff(ls(), "combined_df")); gc()

    return(combined_df)
}


read_fork <- function(fork) {
    f <- read.table(file=fork, header=TRUE, sep=',')

    return(f)
}


model <- function(f) {
    # this gets the basename of the file
    # without extension
    var <- sub('\\..*$','', basename(f))

    # this subsets the fork
    # dataset based on the
    # name of the station
    sub_forks <- subset(forks, Name==var)
    id <- sub_forks[1,1]+1  # this transforms the real id into an index

    # this subssets the closest
    # dataset based on the id in
    # the column X1
    sub_closest <- closest[id,,drop=F]  # this ',,' is important to keep
    ids <- sub_closest[1,1:k]-1  # this transforms the indices into real ids (as a vector)

    # this creates the vector
    # with the names of all
    # closest station to the target
    sub_forks <- subset(forks, ID %in% ids)
    nam <- as.vector(sub_forks$Name)

    # this pre-processes the
    # combined_df into a new dataset
    # with only the closest stations
    n_combined_df <- dplyr::select(combined_df, nam)

    # this creates the x matrix
    # and the y vector used by
    # the model
    x <- dplyr::select(n_combined_df, -one_of(var))
    y <- read.table(file=f, header=FALSE)
    y <- y[2:nrow(y),]    # removes the first line which is just the starting date (YYYYMMMDD)
    y[y == -99] <- NA    # sets -99 as NA

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
        m <- cubist(x=x_train, y=y_train, control=ctrl)

        # this makes the final prediction
        y_mod = predict(m, x)

        # this is to export various
        # types of results
        cat(m$output, file=sub('.txt', '_out.txt', basename(f)))
    })

    if (exists("y_mod")) {
        cat(y_mod, sep="\n", file=basename(f))
        y_mod <- data.frame("l"=unlist(y_mod))
        names(y_mod) <- var
        return(y_mod)
    } else {
        cat(y, sep="\n", file=basename(f))
        y <- data.frame("l"=unlist(y))
        names(y) <- var
        return(y)
    }
}


main <- function(dir) {
    # number of nn that
    # each station will
    # get to use in modelling
    k <<- 11

    # read the fork and
    # calculates all nn
    forks <<- read_fork(file.path(dir, "pcp.txt"))
    closest <<- nn2(data=forks[,3:4], k=k)[[1]]
    closest <<- data.frame(closest)
    # write.csv(closest, 'closest.csv', row.names=TRUE)

    # this makes a list of all fullpath
    # filenames in the directory 'dir',
    # while excluding forks
    files <- c(list.files(path=dir, pattern='^p.*\\.txt$', full.names=TRUE))
    files <- files[!files %in% file.path(dir, "pcp.txt")]
    names(files) <- sapply(files, function(x) {sub('\\..*$','', basename(x))})

    # this gets the number of cores
    cores <- availableWorkers()

    # this reads multiple files
    # into a single dataframe
    # in parallel
    cl <- parallel::makePSOCKcluster(cores, type="PSOCK")   # this initializes all processes
    combined_df <<- make_dataset(files, cl)
    snow::stopCluster(cl)   # this closes all processes

    # this is to release memory
    rm(list=setdiff(ls(), "files")); gc()

    print(c(nrow(combined_df), length(combined_df)))


    for(i in 1:2) {
        # this prints the number
        # of the iteration
        print(paste(c("Processing iteration", i), collapse = " "))

        # this process cubist in
        # parallel for each file
        combined_df <<- lapply(files, model)
        combined_df <<- do.call("cbind", combined_df)

        print(c(nrow(combined_df), length(combined_df)))
    }

    return(NULL)
}


dir <- "/mnt/d/SUPer/Clima/EST_1961_a_2016.ComRAD"     # Linux version
# dir <- "D:/SUPer/Clima/EST_1961_a_2016.ComRAD"    # Windows version
main(dir)