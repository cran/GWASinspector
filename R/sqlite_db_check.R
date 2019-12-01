#' Check the data in the reference database
#'
#' This function displays the summary of the database, including how many tables are in the database file, number of data rows for each data table and the first row of each table.
#'
#' @details This function only checks databases in sqlite format.
#' @param inspector An instance of \linkS4class{Inspector} class. Check \code{\link{setup.inspector}} for more details.
#' @return This function returns a data table including the summary of the specified database. This is necessary to check the consistency and validity of an unknown or new database file.
#' @examples
#' config.file <- get.config(tempdir())
#' inspector <- setup.inspector(config.file , validate = FALSE)
#' # use sample database embedded in the package
#' inspector@@supplementaryFiles$allele_ref_std <- system.file("extdata",
#'                                                             "sample_db.sqlite",
#'                                                              package = "GWASinspector")
#' sqlite.db.check(inspector)
#'
#' @note
#' First column include the names of the tables in database
#'
#' Second column is the number of rows in each table
#'
#' Next columns are the first row of each table
#'
sqlite.db.check <- function(inspector)
{

  if(missing(inspector))
    stop('Function arguments are not set.',call. = FALSE)

  if (!is(inspector, "Inspector"))
    stop("Object must be of class Inspector.", call. = FALSE)


  # check if file is in SQLite format
  path <- inspector@supplementaryFiles$allele_ref_std

  if(tools::file_ext(path) != 'sqlite')
    stop('File is not sqlite type database!', call. = FALSE)

  if (!file.exists(path))
    stop('File not found!', call. = FALSE)

  # message(paste('Reading database file:','"' ,basename(path),'"'))

  check<-tryCatch(
    {
      mydb<-RSQLite::dbConnect(RSQLite::SQLite(),path)
      tables<-RSQLite::dbListTables(mydb)
      # pb <- txtProgressBar(min = 0, max = length(tables), style = 3,width = 65)
      # counter <- 1


      ### ==========================
      rowcount <- data.table()
      rowTable <- data.table()

      for(i in tables){
        query <- sprintf('select count(*) from %s',i)
        rowcount <- rbind(rowcount,
                          as.data.table(RSQLite::dbGetQuery(mydb, query)))

        rowTable <- rbind(rowTable,
                          as.data.table(RSQLite::dbGetQuery(mydb, sprintf('select * from %s LIMIT 1',i))), fill = TRUE)

        # setTxtProgressBar(pb, counter)
        # counter <- counter + 1
      }

      colnames(rowcount) <- 'number of rows'
      RSQLite::dbDisconnect(mydb)
      #close(pb)

      ### ==========================


      tables <- as.data.table(tables)
      colnames(tables) <- 'table name'
      ### ==========================
      number <- as.data.table(1:nrow(tables))
      colnames(number) <- '#'

      ### ==========================
      # message('\n- first column is tables of the database')
      # message( '- second column is the number of rows in that table')
      # message( '- next columns are the first row of each table')

      tbl <- cbind(number,tables,rowcount,rowTable)
     # print(knitr::kable(as.data.table(tbl),format = 'rst',align = 'l'))
      return(as.data.table(tbl))

    },error=function(err)
    {
      stop(paste('Error reading database file:',err$message), call. = FALSE)
    })

}
